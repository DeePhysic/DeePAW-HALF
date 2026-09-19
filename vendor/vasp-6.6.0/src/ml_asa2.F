#include "symbol.inc"
!*********************************************************************
! RCS:  $Id: asa.F,v 1.7 2003/06/27 13:22:13 kresse Exp kresse $
!
!  this modul contains code to calculate Clebsch-Gordan coefficients
!  and code to evaluate the integrals of three spherical (real or
!  complex) harmonics
!
!  The part which calculates the Clebsch-Gordan coefficients
!  was kindly supplied by Helmut Nowotny
!  and is taken from an ASW (augmented spherical waves) code
!  written by Williams, Kuebler and Gelat.
!  port to f90, a lot of extensions and code clean up was
!  done by gK (I think I have also written the version for the real
!  spherical harmonics, but I am not shure ;-).
!
!*********************************************************************

  MODULE ml_asa2
    USE ML_FF_PREC

  CONTAINS

    SUBROUTINE ASA2_ALLOCATE(LMAXCG,ASA)
       USE ML_FF_STRUCT , ONLY : ASA_STUFF
       IMPLICIT NONE
       TYPE (ASA_STUFF)    :: ASA
       INTEGER, INTENT(IN) :: LMAXCG
       ALLOCATE(ASA%YLM3(1:ASA%NCG))
       ALLOCATE(ASA%YLM3I(1:ASA%NCG))
       ALLOCATE(ASA%JL(1:ASA%NCG))
       ALLOCATE(ASA%JS(1:ASA%NCG))
       ALLOCATE(ASA%INDCG(1:ASA%NCG))
       ALLOCATE(ASA%YLM3LOOKUP_TABLE(0:LMAXCG,0:LMAXCG))
    END SUBROUTINE ASA2_ALLOCATE

    SUBROUTINE ASA2_INITIALIZE(ASA)
       USE ML_FF_STRUCT , ONLY : ASA_STUFF
       IMPLICIT NONE
       TYPE (ASA_STUFF)    :: ASA
       ASA%YLM3=0.0_q
       ASA%YLM3I=0.0_q
       ASA%JL=0
       ASA%JS=0
       ASA%INDCG=0
       ASA%YLM3LOOKUP_TABLE=0
    END SUBROUTINE ASA2_INITIALIZE

    SUBROUTINE ASA2_DEALLOCATE(ASA)
       USE ML_FF_STRUCT , ONLY : ASA_STUFF
       IMPLICIT NONE
       TYPE (ASA_STUFF)    :: ASA
       DEALLOCATE(ASA%YLM3)
       DEALLOCATE(ASA%YLM3I)
       DEALLOCATE(ASA%JL)
       DEALLOCATE(ASA%JS)
       DEALLOCATE(ASA%INDCG)
       DEALLOCATE(ASA%YLM3LOOKUP_TABLE)
    END SUBROUTINE ASA2_DEALLOCATE

!**************** FUNCTION CLEBGO ************************************
!
! caculate Clebsch-Gordan-coeff. <J1 J2 M1 M2 I J3 M3>
! using racah-formel
! FAC is a user supplied array containing factorials
!
!*********************************************************************

    FUNCTION CLEBGO(FAC,J1,J2,J3,M1,M2,M3)

      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) FAC(40)
      REAL(q) CLEBGO

      IF(M3/=M1+M2) GO TO 2
      K1=J1+J2-J3+1
      K2=J3+J1-J2+1
      K3=J3+J2-J1+1
      K4=J1+J2+J3+2
      T= (2*J3+1)*FAC(K1)*FAC(K2)*FAC(K3)/FAC(K4)
      K1=J1+M1+1
      K2=J1-M1+1
      K3=J2+M2+1
      K4=J2-M2+1
      K5=J3+M3+1
      K6=J3-M3+1
      T=SQRT(T*FAC(K1)*FAC(K2)*FAC(K3)*FAC(K4)*FAC(K5)*FAC(K6))
      N1=MAX0(J2-J3-M1,J1-J3+M2,0)+1
      N2=MIN0(J1+J2-J3,J1-M1,J2+M2)+1
      IF(N1>N2) GO TO 2
      T1=0.0_q
      DO M=N1,N2
         N=M-1
         K1=J1+J2-J3-N+1
         K2=J1-M1-N+1
         K3=J2+M2-N+1
         K4=J3-J2+M1+N+1
         K5=J3-J1-M2+N+1
         T1=T1+ (1+4*(N/2)-2*N)/(FAC(M)*FAC(K1)*FAC(K2)*FAC(K3) &
              &  *FAC(K4)*FAC(K5))
      ENDDO
      CLEBGO=T*T1
      RETURN
! coefficient is zero, drop back
 2    CONTINUE
      CLEBGO=0.0_q
      RETURN

    END FUNCTION CLEBGO

!**************** FUNCTION CLEBG0 ************************************
!
! calculate Clebsch-Gordan-coeff. <L1 L2 0 0 I L3 0>
! using racah-formel
! FAC is a user supplied array containing factorials
!
!*********************************************************************


    FUNCTION CLEBG0(FAC,L1,L2,L3)
      USE ML_FF_PREC
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER X,P
      REAL(q) FAC(40)
      REAL(q) CLEBG0

      LT=L1+L2+L3
      P=LT/2
      IF(2*P/=LT) GO TO 1
      CLEBG0= SQRT( REAL(2*L3+1,KIND=q)/(LT+1))
      CLEBG0=CLEBG0*FAC(P+1)/SQRT(FAC(2*P+1))
      X=P-L1
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L2
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L3
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      IF(X>2*(X/2)) CLEBG0=-CLEBG0
      RETURN
! coefficient is zero, drop back
 1    CONTINUE
      CLEBG0=0.0_q
      RETURN
    END FUNCTION CLEBG0

!************************* YLM3ST   **********************************
!
! calculate the integral of the product of three real spherical
! harmonics
! i.e    Y_lm Y_l'm' Y_LM
!
! LMAX     max value for l and lp (maximum L is given by triagular rule
!             | l- lp | < L < | l + lp |
!
!*********************************************************************


    SUBROUTINE YLM3ST(LMAX, NCG,YLM3,YLM3I,JL,JS,INDCG,LMAXCG,YLM3LOOKUP_TABLE)
      USE ML_FF_PREC
      IMPLICIT REAL(q) (A-H,O-Z)
      ! Input and output variables
      INTEGER         :: NCG
      REAL(q)         :: YLM3(:)
      REAL(q)         :: YLM3I(:)
      INTEGER         :: JL(:)
      INTEGER         :: JS(:)
      INTEGER         :: INDCG(:)
      INTEGER         :: LMAXCG
      INTEGER         :: YLM3LOOKUP_TABLE(0:,0:)
!
      INTEGER S1,S2,S3,T1,T2,T3
      REAL(q) FAC(40)
      REAL(q), PARAMETER :: SRPI =1.772453850905516027_q
!---------------------------------------------------------------------
! function to evaluate (-1)^I
!---------------------------------------------------------------------
      IF (LMAXCG>0) RETURN
      LMAXCG=LMAX
!---------------------------------------------------------------------
! set up table for factorials
!---------------------------------------------------------------------
      IMAX=30
      FAC(1)=1._q
      DO I=1,IMAX
         FAC(I+1)= I*FAC(I)
      ENDDO

      IC=0

      LMIND=0
!---------------------------------------------------------------------
! loop over l,m         m =-l,+l
! loop over lp,mp       mp=-lp,+lp
!---------------------------------------------------------------------
      DO L1=0,LMAX
      DO L2=0,LMAX
      K2=(2*L1+1)*(2*L2+1)

      DO M1=-L1,L1
      DO M2=-L2,L2

         LMIND=LMIND+1
         INDCG(LMIND)=IC+1

         N1=IABS(M1)
         S1=0
         IF(M1<0) S1=1
         T1=0
         IF(M1==0) T1=1

         N2=IABS(M2)
         S2=0
         IF(M2<0) S2=1
         T2=0
         IF(M2==0) T2=1

!---------------------------------------------------------------------
! for integrals of 3 real spherical harmonics
! two M values are possibly nonzero
!---------------------------------------------------------------------
         IF(M1*M2<0) THEN
            M3=-N1-N2
            M3P=-IABS(N1-N2)
            IF(M3P==0) THEN
               NM3=1
            ELSE
               NM3=2
            ENDIF
         ELSE IF (M1*M2==0) THEN
            M3=M1+M2
            M3P=0     ! Dummy initialization, not used for this case
            NM3=1
         ELSE
            M3=N1+N2
            M3P=IABS(N1-N2)
            NM3=2
         ENDIF

 5       N3=IABS(M3)
         S3=0
         IF(M3<0) S3=1
         T3=0
         IF(M3==0) T3=1

!---------------------------------------------------------------------
! loop over L given by triangular rule
!---------------------------------------------------------------------
         Q1= 1/2._q*SQRT( REAL(K2,KIND=q))*FS(N3+(S1+S2+S3)/2)
         Q2= 1/(SQRT(2._q)**(1+T1+T2+T3))

         DO L3=ABS(L1-L2),L1+L2, 2

            IF(N3>L3) CYCLE
            T=0._q
            IF(N1+N2==-N3) T=T+CLEBG0(FAC(1),L1,L2,L3)
            IF(N1+N2==N3 ) &
     &           T=T+CLEBGO(FAC(1),L1,L2,L3, N1, N2, N3)*FS(N3+S3)
            IF(N1-N2==-N3) &
     &           T=T+CLEBGO(FAC(1),L1,L2,L3, N1,-N2,-N3)*FS(N2+S2)
            IF(N1-N2==N3 ) &
     &           T=T+CLEBGO(FAC(1),L1,L2,L3,-N1, N2,-N3)*FS(N1+S1)
            IC=IC+1

            IF (IC>NCG)  THEN
!               CALL vtutor%error("ERROR: in YLM3ST IC larger than NCG increase NCG")
            ENDIF

            T0=CLEBG0(FAC(1),L1,L2,L3)

            YLM3(IC) = Q1*Q2*T*T0/(SRPI* SQRT( REAL(2*L3+1,KIND=q)))
            IF (T0==0) THEN
               YLM3I(IC)=0
            ELSE
               YLM3I(IC)= T*Q2/Q1/T0*(SRPI* SQRT( REAL(2*L3+1,KIND=q)))
            ENDIF
!           WRITE(*,'(6I4,E14.7)')L3*(L3+1)+M3+1,0,L1,L2,M1,M2,YLM3(IC)

            JL(IC)=L3
            JS(IC)=L3*(L3+1)+M3+1
         ENDDO
! if there is a second M value calculate coefficients for this M
         NM3=NM3-1
         M3=M3P
         IF(NM3>0) GO TO 5

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      INDCG(LMIND+1)=IC+1

      LMIND=0

      DO L1=0,LMAXCG
      DO L2=0,LMAXCG

         YLM3LOOKUP_TABLE(L1,L2)=LMIND
         LMIND=LMIND+(2*L1+1)*(2*L2+1)
      ENDDO
      ENDDO

      RETURN
      
    CONTAINS

      FUNCTION FS(I)
        INTEGER FS,I
        FS=1-2*MOD(I+20,2)
      END FUNCTION FS

    END SUBROUTINE YLM3ST


!************************* YLM3LOOKUP ********************************
!
! function to look up a the startpoint in the array
! YLM3 for two quantumnumbers l lp
!
!*********************************************************************

    SUBROUTINE YLM3LOOKUP(L,LP,LMIND, LMAXCG,YLM3LOOKUP_TABLE)
      USE ML_FF_PREC

      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER             :: LMAXCG
      INTEGER             :: YLM3LOOKUP_TABLE(0:,0:)

      LMIND=YLM3LOOKUP_TABLE(L,LP)

      RETURN
    END SUBROUTINE YLM3LOOKUP


!************************* SETYLM_GRAD *******************************
!
! calculate spherical harmonics
! and the gradient of the spherical harmonics 
! for a set of grid points up to LMAX 
!
! YLM(:,1:LMAX)      = Y_lm(hat x)
! YLMD(:,1:LMAX,1:3) = d Y_lm(hat x)/d hat x_i
!
! hat x must lie on the unit sphere
!
! to obtain the final gradient one has to use
!
! d Y_lm(x)/ d x_i = sum_j d Y_lm(hat x)/ d hat x_j  d hat x_j/ d x_i
! d Y_lm(x)/ d x_i =
!   sum_j (delta_ij - hat x_i hat x_j)/|x| d Y_lm(hat x)/ d hat x_j 
!
! where hat x = x / ||x||
!
! written by Georg Kresse
!
!*********************************************************************

    SUBROUTINE SETYLM_GRAD(LYDIM,INDMAX,YLM,YLMD,X,Y,Z, YLM3I,JS,INDCG,LMAXCG,YLM3LOOKUP_TABLE)
      USE ML_FF_PREC
      USE ML_FF_CONSTANT 
      IMPLICIT NONE

      REAL(q)                :: YLM3I(:)
      INTEGER                :: JS(:)
      INTEGER                :: INDCG(:)
      INTEGER                :: LMAXCG
      INTEGER                :: YLM3LOOKUP_TABLE(0:,0:)  

      INTEGER LYDIM           ! maximum L
      INTEGER INDMAX          ! number of points (X,Y,Z)
      REAL(q) YLM(:,:)        ! spherical harmonics
      REAL(q) YLMD(:,:,:)     ! gradient of spherical harmonics
      REAL(q) X(:),Y(:),Z(:)  ! x,y and z coordinates
! local variables
      REAL(q) FAK
      INTEGER IND,LSET,LM,LP,LMINDX,ISTART,IEND,LNEW,L,M,MP,IC,ITMP
!-----------------------------------------------------------------------
! runtime check of workspace
!-----------------------------------------------------------------------
! gfortran compiler bug workaround suggested by jF
      ITMP=UBOUND(YLM,1)

      FAK=1/(2._q * SQRT(PI))
!-----------------------------------------------------------------------
! here is the code for L=0, hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <0) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,1)  =FAK
        YLMD(IND,1,:)=0
      ENDDO
!-----------------------------------------------------------------------
! here is the code for L=1, once again hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <1) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,2)  = (FAK*SQRT(3._q))*Y(IND)
        YLM(IND,3)  = (FAK*SQRT(3._q))*Z(IND)
        YLM(IND,4)  = (FAK*SQRT(3._q))*X(IND)
        ! gradient with respect to x
        YLMD(IND,2,1)= 0
        YLMD(IND,3,1)= 0
        YLMD(IND,4,1)= (FAK*SQRT(3._q))
        ! gradient with respect to y
        YLMD(IND,2,2)= (FAK*SQRT(3._q))
        YLMD(IND,3,2)= 0
        YLMD(IND,4,2)= 0
        ! gradient with respect to z
        YLMD(IND,2,3)= 0
        YLMD(IND,3,3)= (FAK*SQRT(3._q))
        YLMD(IND,4,3)= 0
      ENDDO
!-----------------------------------------------------------------------
! code for L=2,
!-----------------------------------------------------------------------
      IF (LYDIM <2) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,5)= (FAK*SQRT(15._q))  *X(IND)*Y(IND)
        YLM(IND,6)= (FAK*SQRT(15._q))  *Y(IND)*Z(IND)
        YLM(IND,7)= (FAK*SQRT(5._q)/2._q)*(3*Z(IND)*Z(IND)-1)
        YLM(IND,8)= (FAK*SQRT(15._q))  *X(IND)*Z(IND)
        YLM(IND,9)= (FAK*SQRT(15._q)/2._q)*(X(IND)*X(IND)-Y(IND)*Y(IND))
        ! gradient with respect to x
        YLMD(IND,5,1)= (FAK*SQRT(15._q))  *Y(IND)
        YLMD(IND,6,1)= 0
        YLMD(IND,7,1)= 0
        YLMD(IND,8,1)= (FAK*SQRT(15._q))  *Z(IND)
        YLMD(IND,9,1)= (FAK*SQRT(15._q)/2._q)*2*X(IND)
        ! gradient with respect to y
        YLMD(IND,5,2)= (FAK*SQRT(15._q))  *X(IND)
        YLMD(IND,6,2)= (FAK*SQRT(15._q))  *Z(IND)
        YLMD(IND,7,2)= 0
        YLMD(IND,8,2)= 0
        YLMD(IND,9,2)= (FAK*SQRT(15._q)/2._q)*(-2*Y(IND))
        ! gradient with respect to z
        YLMD(IND,5,3)= 0
        YLMD(IND,6,3)= (FAK*SQRT(15._q))  *Y(IND)
        YLMD(IND,7,3)= (FAK*SQRT(5._q)/2._q)*6*Z(IND)
        YLMD(IND,8,3)= (FAK*SQRT(15._q))  *X(IND)
        YLMD(IND,9,3)= 0
      ENDDO
!-----------------------------------------------------------------------
! initialize all componentes L>2 to zero
!-----------------------------------------------------------------------
      IF (LYDIM <3) GOTO 100
      LSET=2

      DO LM=(LSET+1)*(LSET+1)+1,(LYDIM+1)*(LYDIM+1)
      DO IND=1,INDMAX
        YLM(IND,LM) = 0
        YLMD(IND,LM,1) = 0
        YLMD(IND,LM,2) = 0
        YLMD(IND,LM,3) = 0
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
! for L>2 we use (some kind of) Clebsch-Gordan coefficients
! i.e. the inverse of the integral of three reel sperical harmonics
!      Y_LM = \sum_ll'mm'  C_ll'mm'(L,M) Y_lm Y_l'm'
!-----------------------------------------------------------------------
      LP=1
      DO L=LSET,LYDIM-1
         CALL YLM3LOOKUP(L,LP,LMINDX,LMAXCG,YLM3LOOKUP_TABLE)
         LNEW=L+LP
         DO M = 1, 2*L +1
         DO MP= 1, 2*LP+1
            LMINDX=LMINDX+1

            ISTART=INDCG(LMINDX)
            IEND  =INDCG(LMINDX+1)

            DO IC=ISTART,IEND-1
               LM=JS(IC)
               IF (LM > LNEW*LNEW       .AND. &
                   LM <= (LNEW+1)*(LNEW+1)) THEN
!DIR$ IVDEP
!OCL NOVREC
                  DO IND=1,INDMAX
                     YLM(IND,LM) = YLM(IND,LM)+ &
                         YLM3I(IC)*YLM(IND,L*L+M)*YLM(IND,LP*LP+MP)
                     ! gradient
                     YLMD(IND,LM,1) = YLMD(IND,LM,1)+ &
                         YLM3I(IC)*(YLMD(IND,L*L+M,1)*YLM(IND,LP*LP+MP)+YLM(IND,L*L+M)*YLMD(IND,LP*LP+MP,1))
                     YLMD(IND,LM,2) = YLMD(IND,LM,2)+ &
                         YLM3I(IC)*(YLMD(IND,L*L+M,2)*YLM(IND,LP*LP+MP)+YLM(IND,L*L+M)*YLMD(IND,LP*LP+MP,2))
                     YLMD(IND,LM,3) = YLMD(IND,LM,3)+ &
                         YLM3I(IC)*(YLMD(IND,L*L+M,3)*YLM(IND,LP*LP+MP)+YLM(IND,L*L+M)*YLMD(IND,LP*LP+MP,3))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         ENDDO
       ENDDO

 100  CONTINUE
        
    END SUBROUTINE SETYLM_GRAD


!************************* SETYLM_GRAD2 ******************************
!
! calculate spherical harmonics
! and the gradient of the spherical harmonics
! for a set of grid points up to LMAX
!
! YLM(:,1:LMAX)      = Y_lm(hat x)
! YLMD(:,1:LMAX,1:3) = 
!       sum_j (delta_ij - hat x_i hat x_j)  d Y_lm(hat x)/ d hat x_j
!
! to obtain the final gradient the YLMD must be divided by the |x|
! (see SETYLM_GRAD)
!
! hat x must be lying on the unit sphere, i.e. hat x = x / |x|
!
! written by Georg Kresse
!
!*********************************************************************

    SUBROUTINE SETYLM_GRAD2(LYDIM,INDMAX,YLM,YLMD,X,Y,Z,NORMHELP, YLM3I,JS,INDCG,LMAXCG,YLM3LOOKUP_TABLE)
      USE ML_FF_PREC
      USE ML_FF_CONSTANT 
      IMPLICIT NONE
      INTEGER LYDIM           ! maximum L
      INTEGER INDMAX          ! number of points (X,Y,Z)
      REAL(q) YLM(:,:)        ! spherical harmonics
      REAL(q) YLMD(:,:,:)     ! gradient of spherical harmonics
      REAL(q) X(:),Y(:),Z(:)  ! x,y and z coordinates
      REAL(q) NORMHELP(:)
      REAL(q)                :: YLM3I(:)
      INTEGER                :: JS(:)
      INTEGER                :: INDCG(:)
      INTEGER                :: LMAXCG
      INTEGER                :: YLM3LOOKUP_TABLE(0:,0:)
    ! local
      INTEGER LM, IND
      REAL(q) YLMX, YLMY, YLMZ

      CALL SETYLM_GRAD(LYDIM,INDMAX,YLM,YLMD,X,Y,Z, YLM3I,JS,INDCG,LMAXCG,YLM3LOOKUP_TABLE)
      DO LM=1,(LYDIM+1)*(LYDIM+1)
         DO IND=1,INDMAX
            YLMX=YLMD(IND,LM,1)
            YLMY=YLMD(IND,LM,2)
            YLMZ=YLMD(IND,LM,3)

            YLMD(IND,LM,1) = (YLMX-X(IND)*X(IND)*YLMX-X(IND)*Y(IND)*YLMY-X(IND)*Z(IND)*YLMZ)/NORMHELP(IND)
            YLMD(IND,LM,2) = (YLMY-Y(IND)*X(IND)*YLMX-Y(IND)*Y(IND)*YLMY-Y(IND)*Z(IND)*YLMZ)/NORMHELP(IND)
            YLMD(IND,LM,3) = (YLMZ-Z(IND)*X(IND)*YLMX-Z(IND)*Y(IND)*YLMY-Z(IND)*Z(IND)*YLMZ)/NORMHELP(IND)
         ENDDO
      ENDDO

    END SUBROUTINE SETYLM_GRAD2

  END MODULE ml_asa2
