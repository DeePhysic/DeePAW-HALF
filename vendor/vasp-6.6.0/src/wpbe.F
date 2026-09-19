#include "symbol.inc"

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!=====================================================================
!Functional wPBE short-range for calculation in separated interactions LR-SR
!with modified kernel Short-Range erf(\mu*r)/r in non-spin polarized case
! ref1: "Hybrid Functionals based on a screened Coulomb potential" PRB 2003, vol 118 p 8207
! ref2: "Generalized gradient approximation to the angle- and system-averaged exchange hole"
! JCP 1998, vol 109 p 3313
!Author: Iann Gerber
!Date : 01-04
!=====================================================================
      MODULE wpbe
      USE prec
      USE constant
      USE ldalib
      USE setexm_struct_def
      IMPLICIT NONE
      PRIVATE :: LOGRHO,S
! Table for the range separated enhancement factor
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: FS_RS(:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: SSPLINES_RHO_S(:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: SSPLINES_S_RHO(:,:)
! Dimensions of the table
      INTEGER, PARAMETER, PRIVATE :: NRHO=2000
      REAL(q), PARAMETER, PRIVATE :: INV_NRHOMIN1=1._q/(NRHO-1)
      INTEGER, PARAMETER, PRIVATE :: NS=2000
      REAL(q), PARAMETER, PRIVATE :: INV_NSMIN1=1._q/(NS-1)
      INTEGER, PRIVATE, SAVE :: LOGRHO0

! Constants for enforcement of local Lieb-Oxford bound
! Value of the maximum of s without the rescaling
! define the value of MAXS in the table construction

      REAL(q), PARAMETER,PRIVATE :: STRANS=8.3_q
      REAL(q), PARAMETER,PRIVATE :: SMAX=8.5728844_q
      REAL(q), PARAMETER,PRIVATE :: SCONST=18.79622316_q

      REAL(q), PARAMETER, PRIVATE :: MAXLOGRHO= 22.5_q
      REAL(q), PARAMETER, PRIVATE :: MINLOGRHO=-33.0_q
!      REAL(q), PARAMETER, PRIVATE :: MAXS=(STRANS+3.0_q)
      REAL(q), PARAMETER, PRIVATE :: MAXS=(9.3_q)
      REAL(q), PARAMETER, PRIVATE :: MINS= 0.0_q
      REAL(q), PARAMETER, PRIVATE :: LOGRHOSTEP=(MAXLOGRHO-MINLOGRHO)/(1._q*(NRHO-1))
      REAL(q), PARAMETER, PRIVATE :: INV_LOGRHOSTEP=1._q/LOGRHOSTEP
      REAL(q), PARAMETER, PRIVATE :: SSTEP=MAXS/(1._q*(NS-1))
      REAL(q), PARAMETER, PRIVATE :: INV_SSTEP=1._q/SSTEP
GPACC !$ACC DECLARE CREATE(LOGRHO0,SSPLINES_RHO_S,SSPLINES_S_RHO,FS_RS)
GPOMP !$OMP DECLARE TARGET(LOGRHO0,SSPLINES_RHO_S,SSPLINES_S_RHO,FS_RS)

      REAL(q), PRIVATE, SAVE :: LDASCREEN_WPBE_TABLE=-1
      LOGICAL, PRIVATE, SAVE :: INITIALIZE_WPBE_TABLE=.TRUE.

      CONTAINS

!************************ SUBROUTINE INIT_WPBE *************************
!
! (Re)initialize the wpbe interpolation tables
!
!***********************************************************************

      SUBROUTINE INIT_WPBE(XC)

      USE setexm_struct_def, ONLY : xc_info,XC_DATA_WPBE_TABLE
      USE xc_name
      TYPE (xc_info) XC

! Early exit if possible
      IF (((.NOT.((XC%NXC==1).AND.(XC%ID(1)==ID_XC_PBE))) .AND. &
        &  (.NOT.((XC%NXC==1).AND.(XC%ID(1)==ID_XC_PBESOL)))) .OR. &
        & (XC%LDASCREEN<=0) .OR. XC%LUSE_THOMAS_FERMI) RETURN

! Possibly reinitialize the table
      IF (.NOT.INITIALIZE_WPBE_TABLE) THEN
         IF ((LDASCREEN_WPBE_TABLE/=XC%LDASCREEN) .OR. &
           & (XC_DATA_WPBE_TABLE%NXC/=1) .OR. &
           & (XC_DATA_WPBE_TABLE%ID(1)/=XC%ID(1))) INITIALIZE_WPBE_TABLE=.TRUE.
      ENDIF

! Deallocate the interpolation tables if allocated
      IF (INITIALIZE_WPBE_TABLE.AND.ALLOCATED(FS_RS)) THEN
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
         DEALLOCATE(FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
      ENDIF
! (Re)initialize the wpbe tables
      CALL INIT_WPBE_TABLE(XC)
      RETURN
      END SUBROUTINE INIT_WPBE


!************************ SUBROUTINE INIT_WPBE_TABLE *******************
!
!***********************************************************************

      SUBROUTINE INIT_WPBE_TABLE(XC)

      USE setexm_struct_def, ONLY : xc_info,XC_DATA_WPBE_TABLE
      USE xc_name
      USE tutor, ONLY : vtutor

      IMPLICIT NONE
      TYPE (xc_info) XC
      INTEGER IGRID,JGRID
      REAL(q) RHO
      REAL(q) FXWPBE_SR
! For testing purposes
      REAL(q) DFSR_DR,DFSR_DS,LGRHO,SV

! early exit if possible
      IF (.NOT.INITIALIZE_WPBE_TABLE) RETURN

!$OMP CRITICAL (omp_lck_wpbe_table)
!$    IF (INITIALIZE_WPBE_TABLE) THEN
! Allocate tables
         ALLOCATE(FS_RS(NRHO,NS),SSPLINES_RHO_S(NRHO,NS),SSPLINES_S_RHO(NRHO,NS))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
#ifdef CRAYHIP
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
#endif
! Fill the table containing the short range part of the
! range separated enhancement factor
         DO IGRID=1,NRHO
            IF ((LOGRHO(IGRID)>0._q).AND.((LOGRHO(IGRID)-LOGRHOSTEP)<0._q)) THEN
               LOGRHO0=IGRID
            ENDIF
            RHO=EXP(LOGRHO(IGRID))
            DO JGRID=1,NS
!gK: converted LDASCREEN to  atomic units
!LS chose between PBE (with EP exchange hole) and PBEsol (HJS exchange hole)
               IF ((XC%NXC==1).AND.(XC%ID(1)==ID_XC_PBE)) THEN
                  CALL EXCHWPBE_R(XC%LDASCREEN*AUTOA,RHO,S(JGRID),FXWPBE_SR)
                  FS_RS(IGRID,JGRID)=FXWPBE_SR
               ELSE IF ((XC%NXC==1).AND.(XC%ID(1)==ID_XC_PBESOL)) THEN
                  CALL EXCHWPBEsol_R(XC%LDASCREEN*AUTOA,RHO,S(JGRID),FXWPBE_SR)
                  FS_RS(IGRID,JGRID)=FXWPBE_SR
               ELSE
                  CALL vtutor%error("INIT_WPBE_TABLE: A screening can be applied only when the semilocal component of "//&
                  "the functional is LDA (GGA=CA or PZ), PBE (GGA=PE) or PBEsol (GGA=PS).")
               ENDIF
            ENDDO
         ENDDO
! Calculate the 2nd derivative arrays
         CALL WPBE_SPLIE2
! Store the parameters for which the table was set up
         XC_DATA_WPBE_TABLE=XC
         LDASCREEN_WPBE_TABLE=XC%LDASCREEN
! and flag the tables as initialized
         INITIALIZE_WPBE_TABLE = .FALSE.
!$    ENDIF
!$OMP END CRITICAL (omp_lck_wpbe_table)
GPACC !$ACC UPDATE DEVICE(LOGRHO0,FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
GPOMP !$OMP TARGET UPDATE TO(LOGRHO0,FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)

      RETURN
      END SUBROUTINE INIT_WPBE_TABLE


      FUNCTION LOGRHO(I)
GPACC !$ACC ROUTINE SEQ
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(LOGRHO)
      INTEGER I
      REAL(q) LOGRHO
      LOGRHO=MINLOGRHO+(I-1)*(MAXLOGRHO-MINLOGRHO)/(NRHO-1)
      END FUNCTION LOGRHO

      FUNCTION S(I)
GPACC !$ACC ROUTINE SEQ
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(S)
      INTEGER I
      REAL(q) S
      S=MINS+(I-1)*(MAXS-MINS)/(NS-1)
! Rescaling the s values to ensure the Lieb-Oxford bound for s>8.3
!      IF (S>STRANS) S=SMAX-(SCONST/S**2)
      END FUNCTION S


      SUBROUTINE WPBE_SPLIE2
      IMPLICIT NONE
      INTEGER IGRID,JGRID
      REAL(q), ALLOCATABLE :: XTMP(:),YTMP(:),Y2TMP(:),WORK(:)

! Calculate SSPLINES_RHO_S: second derivative of FS_RS wrt. S at constant RHO
      ALLOCATE(XTMP(NS),YTMP(NS),Y2TMP(NS),WORK(NS))
      DO IGRID=1,NRHO
         DO JGRID=1,NS
            YTMP(JGRID)=FS_RS(IGRID,JGRID)
            XTMP(JGRID)=S(JGRID)
         ENDDO
         CALL WPBE_SPLIN2(XTMP,YTMP,WORK,NS,1E30_q,1E30_q,Y2TMP)
         DO JGRID=1,NS
            SSPLINES_RHO_S(IGRID,JGRID)=Y2TMP(JGRID)
         ENDDO
      ENDDO
      DEALLOCATE(XTMP,YTMP,Y2TMP,WORK)
! Calculate SSPLINES_S_RHO: second derivative of FS_RS wrt. RHO at constant S
      ALLOCATE(XTMP(NRHO),YTMP(NRHO),Y2TMP(NRHO),WORK(NRHO))
      DO IGRID=1,NS
         DO JGRID=1,NRHO
            YTMP(JGRID)=FS_RS(JGRID,IGRID)
            XTMP(JGRID)=LOGRHO(JGRID)
         ENDDO
         CALL WPBE_SPLIN2(XTMP,YTMP,WORK,NRHO,1E30_q,1E30_q,Y2TMP)
         DO JGRID=1,NRHO
            SSPLINES_S_RHO(IGRID,JGRID)=Y2TMP(JGRID)
         ENDDO
      ENDDO
      DEALLOCATE(XTMP,YTMP,Y2TMP,WORK)
      RETURN
      END SUBROUTINE WPBE_SPLIE2


      SUBROUTINE WPBE_SPLIN2(X,Y,U,N,YP1,YPN,Y2)
GPACC !$ACC ROUTINE SEQ
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(WPBE_SPLIN2)
      INTEGER N
      REAL(q) X(N),Y(N),U(N),Y2(N)
      REAL(q) YP1,YPN

      INTEGER I,K
      REAL(q)  P,QN,SIG,UN

! First point
      IF (YP1.GT.0.99E30_q) THEN
         Y2(1)=0; U(1)=0
      ELSE
!        Y2(1) =-0.5_q
         Y2(1) = 0.5_q
         U(1)  = (3._q/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
! Decomposition loop for the tridiagonal alg
      DO I=2,N-1
         SIG   = (X(I)-X(I-1))/(X(I+1)-X(I-1))
         P     = SIG*Y2(I-1)+2._q
         Y2(I) = (SIG-1._q)/P
         U(I)  = (6._q*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
        &        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
! Last point
      IF (YPN.GT.0.99E30_q) THEN
         QN=0; UN=0
      ELSE
         QN = 0.5_q
         UN = (3._q/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1._q)
! Backsubstitution loop of the tridiagonal alg
      DO K=N-1,1,-1
         Y2(K) = Y2(K)*Y2(K+1)+U(K)
      ENDDO

      RETURN
      END SUBROUTINE WPBE_SPLIN2


      SUBROUTINE WPBE_SPLINE(X1,X2,Y)
GPACC !$ACC ROUTINE SEQ
      USE string, ONLY : str
      USE tutor, ONLY : vtutor
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(WPBE_SPLINE)

      REAL(q) X1,X2,Y
      REAL(q) D1X1,D1X2

      INTEGER J,K,JLO,JHI,KLO,KHI,JTMP,KTMP,KMIN

      REAL(q) YYTMPX2(6),Y2X2(6),LOGRHOTMP(6)
      REAL(q) YYTMPX1(6),Y2X1(6),STMP(6),WORK(6)
      REAL(q) H,A,B
      REAL(q) D1HI1,D1LO2,DIFFX1,DIFFX2
      REAL(q) D1A,D1B,H13,X_XJ_3_H,XJ1_X_3_H
      REAL(q) YX2,FX
GPACC !$ACC ROUTINE (WPBE_SPLIN2) SEQ
GPACC !$ACC ROUTINE (LOGRHO) SEQ
! Rescaled the S parameter to lie in the bounds:
      IF (X2>STRANS) X2=SMAX-(SCONST/X2**2)
#ifndef OMP_OFFLOAD
! Check whether X1 and X2 are within bounds
      IF (X1<LOGRHO(1).OR.X1>LOGRHO(NRHO)) THEN
         WRITE(*,*) 'WPBE_SPLINE: ERROR(1), LOG(RHO) out of bounds: ', &
        & X1,LOGRHO(1),LOGRHO(NRHO)
      ENDIF
      IF (X2<S(1).OR.X2>S(NS)) THEN
         WRITE(*,*) 'WPBE_SPLINE: ERROR(1), S out of bounds: ', &
        & X2,S(1),S(NS),EXP(X1)
      ENDIF
#endif
! Find bracketing indices at the rho axis
      IF (X1.GE.0._q) THEN
         JLO=INT(X1/LOGRHOSTEP)+LOGRHO0
      ELSE
         JLO=INT(ABS((X1-MINLOGRHO)/LOGRHOSTEP))+1
      ENDIF
      JHI=JLO+1

      IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
         IF ((LOGRHO(JLO).GT.X1)) THEN
            JLO=JLO-1
         ELSE
            JLO=JLO+1
         ENDIF
         JHI=JLO+1
      ENDIF

      IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
__CPU    CALL vtutor%error('WPBE_SPLINE: ERROR(2), LOG(RHO) out of bounds: ' // &
__CPU       str(X1) // ' ' // str(LOGRHO(JLO)) // ' ' // str(LOGRHO(JHI)) // &
__CPU       ' ' // str(JLO) // ' ' // str(JHI) // ' ' // str(LOGRHO0))
      ENDIF

! Find bracketing indices at the S axis
      KLO=INT((X2/SSTEP))+1
      KHI=KLO+1

      IF (S(KLO).GT.X2) THEN
         KHI=KLO
         KLO=KLO-1
      ENDIF

      IF ((S(KLO).GT.X2).OR.(S(KHI)-X2.LT.-1.0E-8_q)) THEN
__CPU    CALL vtutor%error('WPBE_SPLINE: ERROR(2), S out of bounds: ' // str(X2) // ' ' // &
__CPU       str(S(KLO)) // ' ' // str(S(KHI)) // ' ' // str(KLO) // ' ' // str(KHI))
         IF (S(KLO).GT.X2) THEN
            KLO=NS-4
            KHI=NS-3
         ENDIF
      ENDIF

#if defined(RETAIN_DEAD_CODE)
! Interpolate S on a section of the rows
      DO J=1,6
         ! Calculate function value
         H=S(KHI)-S(KLO)
         A=(S(KHI)-X2)/H
         B=(X2-S(KLO))/H

         JTMP=JLO-3+J
         YYTMPX1(J)=A*FS_RS(JTMP,KLO)+B*FS_RS(JTMP,KHI)+ &
        &  ((A*A*A-A)*SSPLINES_RHO_S(JTMP,KLO)+ &
        &   (B*B*B-B)*SSPLINES_RHO_S(JTMP,KHI))*(H*H)/6._q
      ENDDO
! Calculate a spline from these 6 function values
      DIFFX1=LOGRHOSTEP
      D1LO2 =((YYTMPX1(2)-YYTMPX1(1))/DIFFX1)
      D1HI1 =((YYTMPX1(6)-YYTMPX1(5))/DIFFX1)

      DO J=1,6
         LOGRHOTMP(J)=LOGRHO(JLO-3+J)
      ENDDO
      CALL WPBE_SPLIN2(LOGRHOTMP,YYTMPX1,WORK,6,D1LO2,D1HI1,Y2X1)
#endif

! Now do the same thing along the other dimension
      ! Take care of X2 close to zero
      IF (KLO.LT.3) THEN
         KMIN=2
         YYTMPX2(1)=1._q
         IF (KLO.LT.2) THEN
            KMIN=3
         ENDIF
      ELSE
         KMIN=1
      ENDIF

      DO K=KMIN,6
         ! Calculate function value
         H=LOGRHO(JHI)-LOGRHO(JLO)
         A=(LOGRHO(JHI)-X1)/H
         B=(X1-LOGRHO(JLO))/H

         KTMP=KLO-3+K

         YYTMPX2(K)=A*FS_RS(JLO,KTMP)+B*FS_RS(JHI,KTMP)+ &
        &  ((A*A*A-A)*SSPLINES_S_RHO(KTMP,JLO)+ &
        &   (B*B*B-B)*SSPLINES_S_RHO(KTMP,JHI))*(H*H)/6._q
      ENDDO

      IF (KMIN.GT.1) THEN
         IF (KMIN.EQ.2) THEN
            YYTMPX2(1)=YYTMPX2(3)
         ENDIF
         IF (KMIN.EQ.3) THEN
            YYTMPX2(2)=YYTMPX2(4)
            YYTMPX2(1)=YYTMPX2(5)
         ENDIF
      ENDIF

! Calculate a spline from these 6 function values
      DIFFX2=SSTEP
      D1LO2 =((YYTMPX2(2)-YYTMPX2(1))/DIFFX2)
      D1HI1 =((YYTMPX2(6)-YYTMPX2(5))/DIFFX2)

      IF (KMIN.EQ.1) THEN
         DO K=1,6
            STMP(K)=S(KLO-3+K)
         ENDDO
         CALL WPBE_SPLIN2(STMP,YYTMPX2,WORK,6,D1LO2,D1HI1,Y2X2)
      ELSE
         IF (KLO.EQ.2) THEN
            STMP(1)=S(KLO-1)-SSTEP
         ENDIF
         IF (KLO.EQ.1) THEN
            STMP(1)=S(KLO)-2._q*SSTEP
            STMP(2)=S(KLO)-SSTEP
         ENDIF
         DO K=KMIN,6
            STMP(K)=S(KLO-3+K)
         ENDDO
         CALL WPBE_SPLIN2(STMP,YYTMPX2,WORK,6,D1LO2,D1HI1,Y2X2)
      ENDIF
#if defined(RETAIN_DEAD_CODE)
! Interpolate the first spline for X1 (LOGRHO)
      H=LOGRHOSTEP
! Derivative in X1 direction
      D1A = -1._q/H
      D1B = -D1A
      H13 = (H*H*H)**(-1._q/3._q)

      XJ1_X_3_H=(LOGRHO(JHI)-X1)*(LOGRHO(JHI)-X1)*3._q*H13
      X_XJ_3_H =(X1-LOGRHO(JLO))*(X1-LOGRHO(JLO))*3._q*H13
! dF/dRHO
!      D1X1=(D1A*YYTMPX1(3)+D1B*YYTMPX1(4)+ &
!     &      (H*H/6._q)*((-XJ1_X_3_H-D1A)*Y2X1(3)+ &
!     &      (X_XJ_3_H-D1B)*Y2X1(4)))

! Interpolate the second spline for X2 (S)
#endif
      H=SSTEP
      A=(S(KHI)-X2)/H
      B=(X2-S(KLO))/H
! F(rho,s)
      YX2=A*YYTMPX2(3)+B*YYTMPX2(4)+ &
     &   ((A*A*A-A)*Y2X2(3)+(B*B*B-B)*Y2X2(4))*(H*H)/6._q
#if defined(RETAIN_DEAD_CODE)
! Derivative in X2 direction
      D1A = -1._q/H
      D1B = -D1A
      H13 = (H*H*H)**(-1._q/3._q)

      XJ1_X_3_H=(S(KHI)-X2)*(S(KHI)-X2)*3._q*H13
      X_XJ_3_H =(X2-S(KLO))*(X2-S(KLO))*3._q*H13
! dF/dS
!      D1X2=(D1A*YYTMPX2(3)+D1B*YYTMPX2(4)+ &
!     &      (H*H/6._q)*((-XJ1_X_3_H-D1A)*Y2X2(3)+ &
!     &      (X_XJ_3_H-D1B)*Y2X2(4)))
#endif
      Y=YX2

!     IF (Y.LT.0._q) THEN
!        WRITE(*,*) 'WPBE_SPLINE: Y<1',X1,X2,Y
!     ENDIF
      RETURN
      END SUBROUTINE WPBE_SPLINE

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE WPBE_SPLIN2_N6_SIMD(X, Y, YP1, YPN, Y2, MASK0)
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: X(*), Y(*)
      TYPE(SIMD_REAL64), INTENT(IN)      :: YP1, YPN
      TYPE(SIMD_REAL64), INTENT(OUT)     :: Y2(*)
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: U(6), XN_SUB_XNM1(6), YN_SUB_YNM1(6)
      REAL(q)                            :: P, INV_P, QN, SIG, UN
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5
      INTEGER                            :: I, II, K

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4, SIG, P, INV_P, QN, UN)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         DO K = 2, 6
            XN_SUB_XNM1(K)%X(II) = X(K)%X(II) - X(K - 1)%X(II)
            YN_SUB_YNM1(K)%X(II) = Y(K)%X(II) - Y(K - 1)%X(II)
         ENDDO
! First point
         IF (MASK0%X(II)) THEN
            Y2(1)%X(II) = 0.0_q
         ENDIF
         U(1)%X(II) = 0.0_q
         IF (YP1%X(II) .LE. 0.99E30_q) THEN
            IF (MASK0%X(II)) THEN
               Y2(1)%X(II) = 0.5_q
            ENDIF
            TEMP_1 = 1.0_q / XN_SUB_XNM1(2)%X(II)
            U(1)%X(II) = (3.0_q * TEMP_1) * (YN_SUB_YNM1(2)%X(II) * TEMP_1 - YP1%X(II))
         ENDIF
! Decomposition loop for the tridiagonal alg
         TEMP_2 = YN_SUB_YNM1(2)%X(II) / XN_SUB_XNM1(2)%X(II)
         DO I = 2, 5
            TEMP_1 = TEMP_2
            TEMP_3 = 1.0_q / XN_SUB_XNM1(I + 1)%X(II)
            TEMP_2 = YN_SUB_YNM1(I + 1)%X(II) * TEMP_3
            TEMP_4 = 1.0_q / (X(I + 1)%X(II) - X(I - 1)%X(II))
            SIG = XN_SUB_XNM1(I)%X(II) * TEMP_4
            P = SIG * Y2(I - 1)%X(II) + 2.0_q
            INV_P = 1.0_q / P
            IF (MASK0%X(II)) THEN
               Y2(I)%X(II) = (SIG - 1.0_q) * INV_P
            ENDIF
            U(I)%X(II) = (6.0_q * (TEMP_2 - TEMP_1) * TEMP_4 - SIG * U(I - 1)%X(II)) * INV_P
         ENDDO
! Last point
         QN = 0.0_q
         UN = 0.0_q
         IF (YPN%X(II).LE.0.99E30_q) THEN
            QN = 0.5_q
            UN = (3.0_q * TEMP_3) * (YPN%X(II) - TEMP_2)
         ENDIF
         Y2(6)%X(II) = (UN - QN * U(5)%X(II)) / (QN * Y2(5)%X(II) + 1.0_q)
! Backsubstitution loop of the tridiagonal alg
         DO K = 5, 1, -1
            IF (MASK0%X(II)) THEN
               Y2(K)%X(II) = Y2(K)%X(II) * Y2(K + 1)%X(II) + U(K)%X(II)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE WPBE_SPLIN2_N6_SIMD


      SUBROUTINE WPBE_SPLINE_SIMD(X1_IN, X2_IN, Y, MASK0)
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: X1_IN, X2_IN
      TYPE(SIMD_REAL64), INTENT(OUT)     :: Y
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: X1, X2
      TYPE(SIMD_REAL64)                  :: YYTMPX2(6), Y2X2(6), STMP(6)
      TYPE(SIMD_REAL64)                  :: D1HI1, D1LO2, S_KLO, S_KHI, LOGRHO_JLO, LOGRHO_JHI
      TYPE(SIMD_REAL64)                  :: TEMP_FS_RS_JLO(6), TEMP_FS_RS_JHI(6)
      TYPE(SIMD_REAL64)                  :: TEMP_SSPLINES_S_RHO_JLO(6), TEMP_SSPLINES_S_RHO_JHI(6)
      TYPE(SIMD_INT64)                   :: JLO, JHI, KLO, KHI
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q)                            :: D1X1, D1X2
      REAL(q)                            :: H, INV_H, A, B
      INTEGER                            :: II, J, K, KTMP
      LOGICAL                            :: PRED1, PRED2
      INTEGER, PARAMETER                 :: INT1 = 1, INT3 = 3
      REAL(q), PARAMETER                 :: P_1 = 0.166666666666666657414808_q ! 1._q/6._q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(PRED1, PRED2)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Rescaled the S parameter to lie in the bounds:
         X1%X(II) = X1_IN%X(II)
         X2%X(II) = X2_IN%X(II)
         IF (X2%X(II) .GT. STRANS) X2%X(II) = SMAX - SCONST / (X2%X(II) * X2%X(II))
! Check whether X1 and X2 are within bounds
!        IF (X1.LT.LOGRHO(1).OR.X1.GT.LOGRHO(NRHO)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(1), LOG(RHO) out of bounds: ', &
!           &  X1,LOGRHO(1),LOGRHO(NRHO)
!        ENDIF
!        IF (X2.LT.S(1).OR.X2.GT.S(NS)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(1), S out of bounds: ', &
!           &  X2,S(1),S(NS),EXP(X1)
!        ENDIF
! Find bracketing indices at the rho axis
         IF (X1%X(II) .GE. 0.0_q) THEN
            JLO%X(II) = INT(X1%X(II) * INV_LOGRHOSTEP) + LOGRHO0
         ELSE
            JLO%X(II) = INT(ABS((X1%X(II) - MINLOGRHO) * INV_LOGRHOSTEP)) + 1
         ENDIF
         JHI%X(II) = JLO%X(II) + 1
         LOGRHO_JLO%X(II) = MINLOGRHO + (JLO%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         LOGRHO_JHI%X(II) = MINLOGRHO + (JHI%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         PRED1 = LOGRHO_JLO%X(II) .GT. X1%X(II)
         PRED2 = LOGRHO_JHI%X(II) .LT. X1%X(II)
         IF (PRED1) JLO%X(II) = JLO%X(II) - 1
         IF (PRED2) JLO%X(II) = JLO%X(II) + 1
         IF (PRED1 .OR. PRED2) THEN
            JHI%X(II) = JLO%X(II) + 1
            LOGRHO_JLO%X(II) = MINLOGRHO + (JLO%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
            LOGRHO_JHI%X(II) = MINLOGRHO + (JHI%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         ENDIF

!        IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(2), LOG(RHO) out of bounds: ', &
!           &  X1,LOGRHO(JLO),LOGRHO(JHI),JLO,JHI,LOGRHO0
!           STOP
!        ENDIF
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((LOGRHO_JLO%X(II) .GT. X1%X(II)) .OR. (LOGRHO_JHI%X(II) .LT. X1%X(II))) THEN
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF

! Find bracketing indices at the S axis
         KLO%X(II) = INT(X2%X(II) * INV_SSTEP) + 1
         KHI%X(II) = KLO%X(II) + 1
         S_KLO%X(II) = MINS + (KLO%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
         IF (S_KLO%X(II) .GT. X2%X(II)) THEN
            KHI%X(II) = KLO%X(II)
            KLO%X(II) = KLO%X(II) - 1
            S_KLO%X(II) = MINS + (KLO%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
         ENDIF
         S_KHI%X(II) = MINS + (KHI%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
!        IF ((S(KLO).GT.X2).OR.(S(KHI)-X2.LT.-1.0E-8_q)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(2), S out of bounds: ', &
!           &  X2,S(KLO),S(KHI),KLO,KHI
!           STOP
!           IF (S(KLO).GT.X2) THEN
!              KLO=NS-4
!              KHI=NS-3
!           ENDIF
!        ENDIF

         IF ((S_KLO%X(II) .GT. X2%X(II)) .OR. ((S_KHI%X(II) - X2%X(II)) .LT. -1.0E-8_q)) THEN
            MASK1%X(II) = .FALSE.
         ENDIF

         DO K = 1, 6
            STMP(K)%X(II) = MINS + (KLO%X(II) - 3 + K - 1) * (MAXS - MINS) * INV_NSMIN1
         ENDDO
         IF (KLO%X(II) .EQ. 2) THEN
            STMP(1)%X(II) = MINS - SSTEP
         ENDIF
         IF (KLO%X(II) .EQ. 1) THEN
            STMP(1)%X(II) = MINS - 2.0_q * SSTEP
            STMP(2)%X(II) = MINS - SSTEP
         ENDIF
      ENDDO

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK1%X(II)) THEN
            TRUE_FOR_ANY = .TRUE.
         ENDIF
      ENDDO
      
      IF (.NOT. TRUE_FOR_ANY) RETURN
      

      DO II = 0, SIMD_WIDTH_64BIT - 1
! Now do the same thing along the other dimension
! Take care of X2 close to zero
         DO K = 1, 6
            KTMP = MAX(INT1, KLO%X(II) - INT3 + K)
            IF (MASK1%X(II)) THEN
               TEMP_FS_RS_JLO(K)%X(II) = FS_RS(JLO%X(II), KTMP)
               TEMP_FS_RS_JHI(K)%X(II) = FS_RS(JHI%X(II), KTMP)
               TEMP_SSPLINES_S_RHO_JLO(K)%X(II) = SSPLINES_S_RHO(KTMP, JLO%X(II))
               TEMP_SSPLINES_S_RHO_JHI(K)%X(II) = SSPLINES_S_RHO(KTMP, JHI%X(II))
            ENDIF
         ENDDO
      ENDDO

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(H, INV_H, A, B)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         DO K = 1, 6
! Calculate function value
            H = LOGRHO_JHI%X(II) - LOGRHO_JLO%X(II)
            INV_H = 1.0_q / H
            A = (LOGRHO_JHI%X(II) - X1%X(II)) * INV_H
            B = (X1%X(II) - LOGRHO_JLO%X(II)) * INV_H
            YYTMPX2(K)%X(II) = A * TEMP_FS_RS_JLO(K)%X(II) + B * TEMP_FS_RS_JHI(K)%X(II) + &
            &  ((A * A * A - A) * TEMP_SSPLINES_S_RHO_JLO(K)%X(II) + &
            &   (B * B * B - B) * TEMP_SSPLINES_S_RHO_JHI(K)%X(II)) * (H * H) * P_1
         ENDDO
         IF (KLO%X(II) .EQ. 2) THEN
            YYTMPX2(1)%X(II) = YYTMPX2(3)%X(II)
         ENDIF
         IF (KLO%X(II) .EQ. 1) THEN
            YYTMPX2(2)%X(II) = YYTMPX2(4)%X(II)
            YYTMPX2(1)%X(II) = YYTMPX2(5)%X(II)
         ENDIF
! Calculate a spline from these 6 function values
         D1LO2%X(II) = ((YYTMPX2(2)%X(II) - YYTMPX2(1)%X(II)) * INV_SSTEP)
         D1HI1%X(II) = ((YYTMPX2(6)%X(II) - YYTMPX2(5)%X(II)) * INV_SSTEP)
      ENDDO

      CALL WPBE_SPLIN2_N6_SIMD(STMP, YYTMPX2, D1LO2, D1HI1, Y2X2, MASK1)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(A, B)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Interpolate the second spline for X2 (S)
         A = (S_KHI%X(II) - X2%X(II)) * INV_SSTEP
         B = (X2%X(II) - S_KLO%X(II)) * INV_SSTEP
! F(rho,s)
         IF (MASK1%X(II)) THEN
            Y%X(II) = A * YYTMPX2(3)%X(II) + B * YYTMPX2(4)%X(II) + &
            &  ((A * A * A - A) * Y2X2(3)%X(II) + (B * B * B - B) * Y2X2(4)%X(II)) * (SSTEP * SSTEP) * P_1
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE WPBE_SPLINE_SIMD
#endif


!=====================================================================
!Functional wPBE long-range for calculation in separated interactions LR-SR
!with modified kernel Short-Range erf(\mu*r)/r
! ref1: "Hybrid Functionals based on a screened Coulomb potential" PRB 2003, vol 118 p 8207
! ref2: "Generalized gradient approximation to the angle- and system-averaged exchange hole"
! JCP 1998, vol 109 p 3313
! ref 3 jcp 2004 vol 120 7274
!Author: Iann Gerber
!Date : 05-04
!=====================================================================
      SUBROUTINE EXCHWPBE_R(OMEGA,RHO,SW,FXWPBE_SR)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT OMEGA : SPLITTING RANGE
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxwPBE(s,w,rho)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF,ERRFC
!      USE ldalib
      IMPLICIT NONE
      REAL(q) OMEGA,RHO
      REAL(q) FXWPBE_SR
      REAL(q) F13,F43,F12,F14,F32,F34,F94,F98,F1516
      REAL(q) AX,UM,UK,UL,P0,FXPBE
      REAL(q) A,B,C,D,E
      REAL(q) HA1,HA2,HA3,HA4,HA5
      REAL(q) FC1,FC2,EA1,EA2,EA3,EA4,EA5,EA6,EA7,EA8,EB1,WCUT
      REAL(q) EGSCUT,EGA1,EGA2,EGA3
      REAL(q) EXPCUT,EXEI1,EXEI2,EXEI3,EXEI4
      REAL(q) PI2,SRPI,F89M,SREAL,STRA1S
      REAL(q) A2,A3,A4,A12,A32,A52,W,W2,W3,W4,W5,W6,W7,W8,XKF
      REAL(q) SW,S2,S3,S4,S5,S6
      REAL(q) HNUM,HDEN,H,HNU1S,HDE1S,H1S,F
      REAL(q) HSBW,HSBW2,HSBW3,HSBW4,HSBW6,HSBW12,HSBW32,HSBW52,HSBW72
      REAL(q) DHSB,DHSB2,DHSB3,DHSB4,DHSB5,DHSB12,DHSB32,DHSB52,DHSB72
      REAL(q) DHSB92,HA94,HA942,HA943,HA945,HA9412
      REAL(q) DHS,DHS2,DHS3,DHS4,DHS72,DHS92,DHSW,DHSW2,DHSW52,DHSW72
      REAL(q) GA,GB,EG,TM1,TM2,TM3,TM4,TM5,T10
      REAL(q) EXER,EXHA94,EIHA94,EXEI,T1,PN1,PN2
      REAL(q) F2,F3,F4,F5,F6,F7,F8,F9,T2T9
! Numerical factors
      PARAMETER(F13=1._q/3._q,F43=4._q/3._q,F12=0.5_q,F14=0.25_q)
      PARAMETER(F32=1.5_q,F34=0.75_q,F94=2.25_q,F98=1.125_q,F1516=0.9375_q)
! values for PBE enhancement factor calculation
      PARAMETER(UM=0.2195149727645171_q,UK=0.8040_q,UL=UM/UK)
! Constants  from the PBE hole
      PARAMETER(A=1.0161144_q,B=-0.37170836_q)
      PARAMETER(C=-0.077215461_q,D=0.57786348_q)
      PARAMETER(E=-0.051955731_q)
! Values for H(s)
      PARAMETER(HA1=0.00979681_q,HA2=0.0410834_q,HA3=0.187440_q)
      PARAMETER(HA4=0.00120824_q,HA5=0.0347188_q)
! Values for F(s)
      PARAMETER(FC1=6.4753871_q,FC2=0.47965830_q)
! Coefficients of the erfc(x) expansion (eb1 set later depending on wcut)
      PARAMETER(EA1=-1.128223946706117_q,EA2=1.452736265762971_q)
      PARAMETER(EA3=-1.243162299390327_q,EA4=0.971824836115601_q)
      PARAMETER(EA5=-0.568861079687373_q,EA6=0.246880514820192_q)
      PARAMETER(EA7=-0.065032363850763_q,EA8=0.008401793031216_q)
      PARAMETER(WCUT=14.0_q)
! Constants for polynomial expansion of EG for small s
      PARAMETER(EGSCUT=0.08_q,EGA1=-0.02628417880_q,EGA2=-0.07117647788_q)
      PARAMETER(EGA3=0.08534541323_q)
! Constants for large x in exp(x)*Ei(x)
      PARAMETER(EXPCUT=700_q,EXEI1=4.03640_q,EXEI2=1.15198_q)
      PARAMETER(EXEI3=5.03627_q,EXEI4=4.19160_q)

! General constants
      PI2=PI*PI
      SRPI=SQRT(PI)
      F89M=-8._q/9._q

!----------------------------------------------------------------------
! construct modified-PBE enhancement factor
!----------------------------------------------------------------------
!     INTERMEDIATE VARIABLES

! Calculate prelim variables
      XKF=(3._q*PI2*RHO)**F13
      A2=A*A ; A3=A2*A ; A4=A3*A ; A12=SQRT(A) ; A32=A12*A ; A52=A32*A
      W=OMEGA/XKF ; W2= W*W ; W3=W2*W ; W4=W2*W2 ; W5=W2*W3 ; W6=W3*W3
      W7=W6*W ; W8=W7*W

      S2=SW*SW
      S3=S2*SW
      S4=S2*S2
      S5=S4*SW
      S6=S4*S2

! Calculate H(s) and F(s) for the PBE hole
      HNUM=HA1*S2+HA2*S4
      HDEN=1.0_q+HA3*S4+HA4*S5+HA5*S6
      H=(HNUM)/(HDEN)
      HNU1S=2._q*HA1*SW+4._q*HA2*S3
      HDE1S=4.0_q*HA3*S3+5.0_q*HA4*S4+6.0_q*HA5*S5
      H1S=(HDEN*HNU1S-HNUM*HDE1S)/(HDEN*HDEN)
      F=FC1*H+FC2

! Set exponent of the Gaussian in the approximation of the erfc function
      IF (W<WCUT) THEN
         EB1= 1.455915450052607_q
      ELSE
         EB1=2.0_q
      ENDIF

! Calculate intermediate variables
      HSBW=S2*H+EB1*W2 ; HSBW2=HSBW*HSBW ; HSBW3=HSBW2*HSBW ;
      HSBW4=HSBW2*HSBW2 ; HSBW6=HSBW3*HSBW3
      HSBW12=SQRT(HSBW) ; HSBW32=HSBW12*HSBW ; HSBW52=HSBW32*HSBW
      HSBW72=HSBW52*HSBW

      DHSB=D+S2*H+EB1*W2 ; DHSB2=DHSB*DHSB ; DHSB3=DHSB2*DHSB
      DHSB4=DHSB2*DHSB2 ; DHSB5=DHSB4*DHSB
      DHSB12=SQRT(DHSB) ; DHSB32=DHSB12*DHSB ; DHSB52=DHSB32*DHSB
      DHSB72=DHSB52*DHSB ; DHSB92=DHSB72*DHSB

      HA94=F94*HSBW/A ; HA942=HA94*HA94 ; HA943=HA942*HA94
      HA945=HA943*HA942 ; HA9412=SQRT(HA94)

      DHS=D+S2*H ; DHS2=DHS*DHS ; DHS3=DHS2*DHS ; DHS4=DHS2*DHS2
      DHS72=DHS3*SQRT(DHS) ; DHS92=DHS72*DHS

      DHSW=DHS +W2 ; DHSW2=DHSW*DHSW ; DHSW52=SQRT(DHSW)*DHSW2
      DHSW72=DHSW52*DHSW

! Calculate G(s) using expansion for small s if necessary
      IF (SW>EGSCUT) THEN
         GA=SRPI*(15._q*E+6.0_q*C*(1.0_q+F*S2)*DHS + 4.0_q*B*(DHS2)&
        &   +8.0_q*A*(DHS3)) * (1.0_q/(16._q*DHS72))         &
        &   -F34*PI*SQRT(A)*EXP(F94*H*S2/A)* (1._q-ERRF(F32*SW*SQRT(H/A)))
         GB=F1516*SRPI*S2/DHS72
         EG=-(F34*PI+GA)/GB
      ELSE
         EG=EGA1+EGA2*S2+EGA3*S4
      ENDIF
! calculate the terms needed in any case

      TM2=(DHS2*B + DHS*C +2._q*E +DHS*S2*C*F +2._q*S2*EG )/2._q/DHS3
      TM3=-W*(4._q*DHSW2*B +6._q*DHSW*C + 15._q*E + 6.0_q*DHSW*S2*C*F + &
     &        15._q*S2*EG)/8._q/DHS/DHSW52
      TM4=-W3*(DHSW*C + 5._q*E + DHSW*S2*C*F + 5.0_q*S2*EG)/2._q/DHS2/DHSW52
      TM5=-W5*(E+S2*EG)/DHS3/DHSW52

! Calculate t10 unless that would generate a division by zero

      IF ((SW>0.0_q).OR.(W>0.0_q)) THEN
         T10=F12*A*LOG(HSBW/DHSB)
      ENDIF
! Calculate exp(x)*f(x) depending on the size of x

      IF (HA94<EXPCUT) THEN
         EXER=PI*EXP(HA94)*(ERRFC(HA9412))
         EXHA94=EXP(HA94)
         EIHA94=-EXPINT(1,HA94)
         EXEI=EXHA94*EIHA94
      ELSE
         EXER=PI*(1._q/(SRPI*HA9412)-1._q/(2._q*SQRT(PI*HA943))+ &
        &         3._q/(4._q*SQRT(PI*HA945)))
         EXEI=-(1._q/HA94)*(HA942+EXEI1*HA94+EXEI2)/            &
        &      (HA942+EXEI3*HA94+EXEI4)
      ENDIF
      IF (W==0.0_q) THEN
!Fall back to the PBE hole expression
         T1=-F12*A*EXEI
         IF (SW>0.0_q) THEN
            TM1=T1+T10
            FXWPBE_SR=F89M*(TM1+TM2)
          ELSE
            FXWPBE_SR=1._q
         ENDIF
      ELSE IF(W>WCUT) THEN
! Use simple gaussian approximation for large w
         TM1=-F12*A*(EXEI+LOG(DHSB)-LOG(HSBW))
         FXWPBE_SR=F89M*(TM1+TM2+TM3+TM4+TM5)
      ELSE
! For everything else use the full blown expression
!
! First calculate the polynomials for the first term
         PN1=-F32*EA1*A12*W + 27._q*EA3*W3/(8._q*A12)-243._q*EA5*W5/ &
        &     (32._q*A32) + 2187._q*EA7*W7/(128._q*A52)
         PN2=-A + F94*EA2*W2 - 81._q*EA4*W4/(16.0_q*A) + &
        &    729._q*EA6*W6/(64._q*A2) - 6561._q*EA8*W8/(256._q*A3)

! The first term is
         T1=F12*(PN1*EXER+PN2*EXEI)
! The factors for the main polynomials in w
         F2= F12*EA1*SRPI*A/DHSB12
         F3= F12*EA2*A/DHSB
         F4= EA3*SRPI*(-F98/HSBW12+F14*A/DHSB32)
         F5= EA4*(1._q/128._q)*(-144._q*(1._q/HSBW)+64._q*(1._q/DHSB2)*A)
         F6= EA5*(3._q*SRPI*(3._q*DHSB52*(9.0_q*HSBW-2._q*A)+4.0_q  &
        &            *HSBW32*A2))/(32._q*DHSB52*HSBW32*A)
         F7= EA6*(((32._q*A)/DHSB3 + (-36._q+(81._q*S2*H)/A)/HSBW2))&
        &      /32._q
         F8= EA7*(-3._q*SRPI*(-40._q*HSBW52*A3+9.0_q*DHSB72*(27._q  &
        &      *HSBW2-6.0_q*HSBW*A+4._q*A2)))/(128._q*DHSB72*HSBW52*A2)
         F9= (324._q*EA6*EB1*DHSB4*HSBW*A + EA8*(384._q*HSBW3*A3    &
        &        +DHSB4*(-729._q*HSBW2+324._q*HSBW*A-288._q*A2)))      &
        &        /(128._q*DHSB4*HSBW3*A2)

         T2T9= F2*W+F3*W2+F4*W3+F5*W4+F6*W5+F7*W6+F8*W7+f9*W8

! The final value of the first term for 0<omega<wcut is
         TM1= T1+ T2T9 +T10
         FXWPBE_SR=F89M*(TM1+TM2+TM3+TM4+TM5)
      ENDIF
!      write (*,*) RHO,S,FXWPBE_SR
      RETURN
      END SUBROUTINE EXCHWPBE_R

!LS
!
!=====================================================================
!Range-separated enhancementfactor for hybrid functionals
! ref1: "GGA model exchange holes for range-separated hybrids",
!        J.Chem.Phys. 128, 194105(2008)
! not yet implemented: stable version for \nu --> infinity
! J. Chem Theo. Comput. 5, 754 (2009)
!
!Date : 05-2009 ; Laurids Schimka
!=====================================================================
      SUBROUTINE EXCHWPBEsol_R(OMEGA,RHO,SW,FXWPBE_SR)
      USE prec
      USE constant
!      USE ldalib
      IMPLICIT NONE
      REAL(q) OMEGA,RHO
      REAL(q) FXWPBE_SR
      REAL(q) F13,F43,F12,F14,F32,F34,F94,F98,F1516,F54
      REAL(q) ABAR,B,C,D,E
      REAL(q) HA2,HA3,HA4,HA5,HA6,HA7
      REAL(q) HB1,HB2,HB3,HB4,HB5,HB6,HB7,HB8,HB9
      REAL(q) PI2,SRPI,F89M,F25M,F415M,F65M,F45M,F125M
      REAL(q) F49M,F32M,F158M,F38M
      REAL(q) XKF
      REAL(q) SW,S2,S3,S4,S5,S6,S7,S8,S9
      REAL(q) nu,eta,lambda,zeta,chi
      REAL(q) lambda2,lambda3,lambda72
      REAL(q) chi2,chi3,chi5
      REAL(q) zeta12,eta12,nu2
      REAL(q) HNUM,HDEN,H,FBAR,EGBAR
      REAL(q) T1,T2,T3,T4,T5,T6



     ! Numerical factors
      PARAMETER(F13=1._q/3._q,F43=4._q/3._q,F12=0.5_q,F14=0.25_q)
      PARAMETER(F32=1.5_q,F34=0.75_q,F94=2.25_q,F98=1.125_q,F1516=0.9375_q)
      PARAMETER(F54=1.25_q)
     ! Constants  from the PBE(sol) hole
      PARAMETER(ABAR=0.757211_q,B=-0.106364_q,C=-0.118649_q,D=0.609650_q,E=-0.0477963_q)
! this parameters are for the PBE exchange hole
     ! Values for H(s)-PBE; values for PBEsol,B88 and B97x are also available
!      PARAMETER(HA2=0.0159941_q,HA3=0.0852995_q,HA4=-0.160368_q,HA5=0.152645_q)
!      PARAMETER(HA6=-0.0971263_q,HA7=0.0422061_q)
!      PARAMETER(HB1=5.33319_q,HB2=-12.4780_q,HB3=11.0988_q,HB4=-5.11013_q,HB5=1.71468_q)
!      PARAMETER(HB6=-0.610380_q,HB7=0.307555_q,HB8=-0.0770547_q,HB9=0.0334840_q)

     ! Parameters for PBEsol
      PARAMETER(HA2=0.0047333_q,HA3=0.0403304_q,HA4=-0.0574615_q,HA5=0.0435395_q)
      PARAMETER(HA6=-0.0216251_q,HA7=0.0063721_q)
      PARAMETER(HB1=8.52056_q,HB2=-13.9885_q,HB3=9.28583_q,HB4=-3.27287_q,HB5=0.843499_q)
      PARAMETER(HB6=-0.235543_q,HB7=0.0847074_q,HB8=-0.0171561_q,HB9=0.0050552_q)

     ! General constants
      PI2=PI*PI
      SRPI=SQRT(PI)
      F89M=-8._q/9._q;F25M=-2._q/5._q;F415M=-4._q/15._q;F65M=-6._q/5._q
      F45M=-4._q/5._q;F125M=-12._q/5._q;F49M=-4._q/9._q;F32M=-3._q/2._q
      F158M=-15._q/8._q;F38M=-3._q/8._q

!----------------------------------------------------------------------
! construct modified-PBE enhancement factor
!----------------------------------------------------------------------
!     INTERMEDIATE VARIABLES

! Calculate prelim variables
      XKF=(3._q*PI2*RHO)**F13  !(k_f)
      S2=SW*SW;S3=S2*SW;S4=S2*S2;S5=S4*SW;S6=S4*S2;S7=S6*SW;S8=S4*S4;S9=S8*SW
! Calculate H(s)
      HNUM=HA2*S2+HA3*S3+HA4*S4+HA5*S5+HA6*S6+HA7*S7
      HDEN=1._q+HB1*SW+HB2*S2+HB3*S3+HB4*S4+HB5*S5+HB6*S6+HB7*S7+HB8*S8+HB9*S9
      H=(HNUM)/(HDEN)
! Calculate Fbar(s)
      FBAR=1._q-1._q/(27._q*C)*S2/(1._q+S2/4._q)-1._q/(2._q*C)*S2*H
! Calculate chi,lambda,nu,zeta,eta
      nu=OMEGA/XKF;nu2=nu*nu
      zeta=S2*H;zeta12=SQRT(zeta)
      eta=ABAR+S2*H;eta12=SQRT(eta)
      lambda=D+S2*H;lambda2=lambda*lambda;lambda3=lambda2*lambda
      lambda72=lambda3*SQRT(lambda)
      chi=nu/(SQRT(lambda+nu2));chi2=chi*chi;chi3=chi2*chi;chi5=chi2*chi3
! Calculate EGbar(s)
      EGBAR=F25M*C*FBAR*lambda+F415M*B*lambda2+F65M*ABAR*lambda3+F45M*SRPI*lambda72+ &
     &F125M*lambda72*(zeta12-eta12)
! Calculate terms of FXWPBE_SR
      T1=ABAR+F49M*B/lambda*(1._q-chi)
      T2=F49M*C*FBAR/lambda2*(1._q+F32M*chi+F12*chi3)
      T3=F89M*EGBAR/lambda3*(1._q+F158M*chi+F54*chi3+F38M*chi5)
      T4=2._q*nu*(SQRT(zeta+nu2)-SQRT(eta+nu2))
      T5=2._q*zeta*LOG((nu+SQRT(zeta+nu2))/(nu+SQRT(lambda+nu2)))
      T6=-2._q*eta*LOG((nu+SQRT(eta+nu2))/(nu+SQRT(lambda+nu2)))
! Calculate FXWPBE_SR
      FXWPBE_SR=T1+T2+T3+T4+T5+T6
!      WRITE(77,'(5G16.7)') OMEGA,RHO,nu,SW,FXWPBE_SR
      RETURN
      END SUBROUTINE EXCHWPBEsol_R
!End LS

      FUNCTION EXPINT(N,X)
!
!       ============================================
!       Purpose: Compute exponential integral EN(x)
!       ============================================
!
      USE prec
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      REAL(q) EXPINT
      REAL(q) X,EULER,EPS,FPMIN
      REAL(q) B,C,D,H,A,DEL,FACT,PSI
      INTEGER N,MAXIT,I,II,NM1

      PARAMETER (EPS=1E-9_q,FPMIN=1E-30_q)
      PARAMETER (MAXIT=100)
      PARAMETER (EULER=0.577215664901532860606512_q)

      NM1=N-1
      IF ((N<0.OR.X<0.0_q).OR.(X==0._q.AND.(N==0.0_q.OR.N==1))) THEN
         CALL vtutor%stopCode()
      ELSE IF (N==0) THEN
         EXPINT=EXP(-X)/X
      ELSE IF (X==0.0_q) THEN
         EXPINT=1._q/NM1
      ELSE IF (X>1.0_q) THEN
         B=X+N
         C=1._q/FPMIN
         D=1.0_q/B
         H=D
         DO I=1,MAXIT
            A=-I*(NM1+I)
            B=B+2._q
            D=1.0_q/(A*D+B)
            C=B+A/C
            DEL=C*D
            H=H*DEL
            IF (ABS(DEL-1._q)< EPS) THEN
               EXPINT=H*EXP(-X)
               RETURN
            ENDIF
         ENDDO
         CALL vtutor%stopCode()
      ELSE
         IF (NM1/=0) THEN
            EXPINT=1.0_q/NM1
         ELSE
            EXPINT=-LOG(X)-EULER
         ENDIF
         FACT=1._q
         DO I=1,MAXIT
            FACT=-FACT*X/I
            IF (I/=NM1) THEN
               DEL=-FACT/(I-NM1)
            ELSE
               PSI=-EULER
               DO II=1,NM1
                  PSI=PSI+1._q/II
               ENDDO
               DEL=FACT*(LOG(X)+PSI)
            ENDIF
            EXPINT=EXPINT+DEL
            IF (ABS(DEL)<ABS(EXPINT)*EPS) RETURN
         ENDDO
         CALL vtutor%stopCode()
      ENDIF
      RETURN
      END FUNCTION EXPINT

      END MODULE wpbe
