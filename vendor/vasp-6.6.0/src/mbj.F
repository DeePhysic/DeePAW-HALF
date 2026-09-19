#include "symbol.inc"

  MODULE mbj
    USE prec
    USE xc_name

    ! MBJ related constants and variables
    LOGICAL, SAVE :: LMETA_SC_CMBJ=.FALSE.
    REAL(q), SAVE :: CMBJ=-1,CMBJA,CMBJB,CMBJE
    REAL(q),ALLOCATABLE, SAVE :: CMBJ_TYP(:)
    RGRID, ALLOCATABLE, SAVE, PRIVATE :: CMBJ_AUX(:)           ! grid quantity, usually stores C values for MBJ
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: CMBJ_AUG(:)          ! CMBJ averaged in each atomic sphere
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: WEIGHT_AUG_AUG(:,:)  ! weights for averaging over aug. spheres
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: WEIGHT_AUG_PW(:,:)   ! weights for averaging over grid points
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: NORM_AUG_PW(:)       ! norm of weights
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: NORM_PW_AUG(:)       ! norm of weights
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: NORM_AUG_AUG(:)      ! norm of weights
    REAL(q),ALLOCATABLE, SAVE, PRIVATE :: SCA_VAC(:)           ! scaling factor for vacuum
    REAL(q), SAVE, PRIVATE :: VACRAT = 0.0_q                   ! ratio of the vacuum over the whole volume, calculated selfconsistently
    REAL(q), SAVE, PRIVATE :: GRHO_OVER_RHO_PW = 0.0_q
    REAL(q), SAVE :: GRHO_OVER_RHO_ONE_CENTER = 0.0_q
    REAL(q),ALLOCATABLE, SAVE :: GRHO_OVER_RHO_AUG(:)          ! for averages in each atomic sphere
    REAL(q), SAVE :: SMBJ = 2.0_q                              ! spread of the Gaussian in Angstrom
    REAL(q), SAVE :: RSMBJ = 7.0_q                             ! density threshold for PAW contribution against vacuum in bohr
    REAL(q), SAVE :: MBJ_SMALL = 1.E-8_q                       ! "SMALL" value where to cut RHO

GPACC !$ACC DECLARE CREATE(CMBJ,CMBJ_TYP,CMBJ_AUX,LMETA_SC_CMBJ,GRHO_OVER_RHO_PW,MBJ_SMALL)
GPOMP !$OMP DECLARE TARGET(CMBJ,CMBJ_TYP,CMBJ_AUX,LMETA_SC_CMBJ,GRHO_OVER_RHO_PW,MBJ_SMALL)


    CONTAINS

!***********************************************************************
!
!***********************************************************************

      SUBROUTINE GET_CMBJ_PW(XC,I,CVMBJ)
GPACC !$ACC ROUTINE SEQ
      USE setexm_struct_def, only : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET (GET_CMBJ_PW)
      REAL(q) CVMBJ
      INTEGER I
      TYPE (xc_info) XC
      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) THEN
         IF(CMBJ_TYP(1)/=-1 .OR. XC%ID(1)==ID_XC_LMBJ) CVMBJ=CMBJ_AUX(I)
      ENDIF
      RETURN
      END SUBROUTINE GET_CMBJ_PW


!***********************************************************************
!
!***********************************************************************

      SUBROUTINE GET_CMBJ_RAD(XC,N,CVMBJ)
      USE setexm_struct_def, only : xc_info
      IMPLICIT NONE
      REAL(q) CVMBJ
      INTEGER N
      TYPE (xc_info) XC
      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) THEN
         IF (CMBJ_TYP(1)/=-1) CVMBJ=CMBJ_TYP(N)
         IF (XC%ID(1)==ID_XC_LMBJ) CVMBJ=CMBJ_AUG(N)
      ENDIF
      RETURN
      END SUBROUTINE GET_CMBJ_RAD


!***********************************************************************
!
!***********************************************************************
      FUNCTION LscMBJ()
GPACC !$ACC ROUTINE SEQ
      IMPLICIT NONE
      LOGICAL LscMBJ
GPOMP !$OMP DECLARE TARGET (LscMBJ)
      LscMBJ=LMETA_SC_CMBJ
      END FUNCTION LscMBJ


!***********************************************************************
!
!***********************************************************************
      SUBROUTINE SUM_GRHO_OVER_RHO_PW(RHO,GRHO)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET (SUM_GRHO_OVER_RHO_PW)
      REAL(q) RHO,GRHO
__CPU !$OMP ATOMIC
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
      GRHO_OVER_RHO_PW=GRHO_OVER_RHO_PW+GRHO/RHO
GPOMP !$OMP END ATOMIC
      RETURN
      END SUBROUTINE SUM_GRHO_OVER_RHO_PW

!***********************************************************************
!
!***********************************************************************

      SUBROUTINE CREATE_CMBJ_AUX(GRIDC,T_INFO,LATT_CUR,XC)

__GPU USE moffload

      USE poscar
      USE mgrid
      USE lattice
      USE setexm_struct_def, only : xc_info
      IMPLICIT NONE
      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (latt)        LATT_CUR
      TYPE (xc_info)     XC

      INTEGER NG, NC, N1, N2, N3, NX, NY, NZ, NI, NJ, NI_FOUND, I
      INTEGER I1,I2,I3,M1,M2,M3
      REAL(q) X, Y, Z, X1, X2, X3, D, D2, DMIN, TSMBJ2, WTMP
      REAL(q) NORM_PW_AUG_TMP,NORM_AUG_PW_TMP,NORM_AUG_AUG_TMP,WEIGHT_AUG_PW_TMP,WEIGHT_AUG_AUG_TMP

      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) THEN

         IF ((CMBJ_TYP(1)==-1).AND.(XC%ID(1)==ID_XC_MBJ)) RETURN
         IF (ALLOCATED(CMBJ_AUX)) THEN
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: CMBJ_AUX)
            DEALLOCATE(CMBJ_AUX)
         ENDIF

         IF (ALLOCATED(NORM_PW_AUG)) DEALLOCATE(NORM_PW_AUG)
         IF (ALLOCATED(NORM_AUG_PW)) DEALLOCATE(NORM_AUG_PW)
         IF (ALLOCATED(NORM_AUG_AUG)) DEALLOCATE(NORM_AUG_AUG)
         IF (ALLOCATED(WEIGHT_AUG_PW)) DEALLOCATE(WEIGHT_AUG_PW)
         IF (ALLOCATED(WEIGHT_AUG_AUG)) DEALLOCATE(WEIGHT_AUG_AUG)
         IF (ALLOCATED(SCA_VAC)) DEALLOCATE(SCA_VAC)
         ALLOCATE(CMBJ_AUX(GRIDC%RL%NP))
! For OpenMP Offloading, allocatable arrays must be brought to the device,
! even if they are declared as OMP TARGET.
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:CMBJ_AUX)
#ifdef CRAYHIP
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CMBJ_AUX)
#endif

         ALLOCATE(NORM_PW_AUG(GRIDC%RL%NP))
         ALLOCATE(NORM_AUG_PW(T_INFO%NIONS))
         ALLOCATE(NORM_AUG_AUG(T_INFO%NIONS))
         ALLOCATE(WEIGHT_AUG_PW(GRIDC%RL%NP,T_INFO%NIONS))
         ALLOCATE(WEIGHT_AUG_AUG(T_INFO%NIONS,T_INFO%NIONS))
         ALLOCATE(SCA_VAC(GRIDC%RL%NP))

         IF (XC%ID(1)==ID_XC_MBJ) THEN
         ! loop over all grid points
         DO NC=1,GRIDC%RL%NCOL
         N2= GRIDC%RL%I2(NC)
         N3= GRIDC%RL%I3(NC)

         DO N1=1,GRIDC%RL%NROW
            NG=(NC-1)*GRIDC%RL%NROW+N1
            IF (GRIDC%RL%NFAST==3) THEN
               NX=N2
               NY=N3
               NZ=N1
            ELSE
               NX=N1
               NY=N2
               NZ=N3
            ENDIF

            ! now search the ion that is closest to this grid point
            DMIN=1000
            NI_FOUND=-1

            DO NI=1,T_INFO%NIONS
               X1=MOD((REAL(NX,q)-1)*(1._q/GRIDC%NGX)-T_INFO%POSION(1,NI)+10.5_q,1._q)-0.5_q
               X2=MOD((REAL(NY,q)-1)*(1._q/GRIDC%NGY)-T_INFO%POSION(2,NI)+10.5_q,1._q)-0.5_q
               X3=MOD((REAL(NZ,q)-1)*(1._q/GRIDC%NGZ)-T_INFO%POSION(3,NI)+10.5_q,1._q)-0.5_q

               X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
               Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
               Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

               D=SQRT(X*X+Y*Y+Z*Z)
            !   WRITE(*,'(7F14.7)') X1, X2, X3, X, Y, Z, D

               IF (D<DMIN) THEN
                  DMIN=D
                  NI_FOUND=NI
               ENDIF
            ENDDO
            IF (NI_FOUND<0) THEN
               WRITE(0,*) 'internal error in CREATE_CMBJ_AUX: NI_FOUND not correct', NI_FOUND
            ENDIF
            CMBJ_AUX(NG)=CMBJ_TYP(T_INFO%ITYP(NI_FOUND))
         !   WRITE(*,*) 'found',NG, KINEDEN%AUX(NG)
         ENDDO
         ENDDO

         ELSEIF (XC%ID(1)==ID_XC_LMBJ) THEN
            ! TR: local SC MBJ version
            ! find the weights for averaging over ions and PW grid points

#ifdef OFFLOADING
            CALL OFFLOAD_ON_PUSH(.TRUE.)
            CALL COPYIN_TYPED_VAR(GRIDC)
#endif
GPACC !$ACC ENTER DATA COPYIN(T_INFO,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:T_INFO,LATT_CUR) IF(OFFLOAD_ON)

            TSMBJ2 = 2._q*SMBJ*SMBJ
            M1=NINT(9._q*SMBJ*LATT_CUR%BNORM(1))
            M2=NINT(9._q*SMBJ*LATT_CUR%BNORM(2))
            M3=NINT(9._q*SMBJ*LATT_CUR%BNORM(3))
GPACC !$ACC ENTER DATA COPYIN(TSMBJ2,M1,M2,M3) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:TSMBJ2,M1,M2,M3) IF(OFFLOAD_ON)

GPACC !$ACC ENTER DATA CREATE(NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) IF(OFFLOAD_ON)
            
GPACC !$ACC KERNELS PRESENT(NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            NORM_PW_AUG = 0._q
            NORM_AUG_PW = 0._q
            WEIGHT_AUG_PW = 0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(TSMBJ2,M1,M2,M3) &
GPACC !$ACC PRIVATE(N2,N3,NG,NX,NY,NZ,NORM_PW_AUG_TMP,NORM_AUG_PW_TMP,WEIGHT_AUG_PW_TMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:TSMBJ2,M1,M2,M3) &
GPOMP !$OMP PRIVATE(N2,N3,NG,NX,NY,NZ,NORM_PW_AUG_TMP,NORM_AUG_PW_TMP,WEIGHT_AUG_PW_TMP) IF(OFFLOAD_ON)
            ! loop over all ions
            DO NI=1,T_INFO%NIONS
               DO NC=1,GRIDC%RL%NCOL
__CPU             N2= GRIDC%RL%I2(NC)
__CPU             N3= GRIDC%RL%I3(NC)
       
               DO N1=1,GRIDC%RL%NROW
__GPU             N2= GRIDC%RL%I2(NC)
__GPU             N3= GRIDC%RL%I3(NC)
                  NG=(NC-1)*GRIDC%RL%NROW+N1
                  IF (GRIDC%RL%NFAST==3) THEN
                     NX=N2
                     NY=N3
                     NZ=N1
                  ELSE
                     NX=N1
                     NY=N2
                     NZ=N3
                  ENDIF

                  NORM_PW_AUG_TMP = 0._q
                  NORM_AUG_PW_TMP = 0._q
                  WEIGHT_AUG_PW_TMP = 0._q
                  ! lattice summation (within the range of the Gaussian)
GPACC !$ACC LOOP COLLAPSE(3) REDUCTION(+:NORM_PW_AUG_TMP,NORM_AUG_PW_TMP,WEIGHT_AUG_PW_TMP) &
GPACC !$ACC PRIVATE(X1,X2,X3,X,Y,Z,D2,WTMP)
GPOMP !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:NORM_PW_AUG_TMP,NORM_AUG_PW_TMP,WEIGHT_AUG_PW_TMP) &
GPOMP !$OMP PRIVATE(X1,X2,X3,X,Y,Z,D2,WTMP)
                  DO I1=-M1,M1
                   DO I2=-M2,M2
                    DO I3=-M3,M3
                     X1=I1+MODULO((REAL(NX,q)-1)*(1._q/GRIDC%NGX)-T_INFO%POSION(1,NI)+10.5_q,1._q)-0.5_q
                     X2=I2+MODULO((REAL(NY,q)-1)*(1._q/GRIDC%NGY)-T_INFO%POSION(2,NI)+10.5_q,1._q)-0.5_q
                     X3=I3+MODULO((REAL(NZ,q)-1)*(1._q/GRIDC%NGZ)-T_INFO%POSION(3,NI)+10.5_q,1._q)-0.5_q
       
                     X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                     Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                     Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
  
                     D2=X*X+Y*Y+Z*Z
            
                     WTMP = EXP(-D2/TSMBJ2)

                     NORM_PW_AUG_TMP = NORM_PW_AUG_TMP + WTMP
                     NORM_AUG_PW_TMP = NORM_AUG_PW_TMP + WTMP
                     WEIGHT_AUG_PW_TMP = WEIGHT_AUG_PW_TMP + WTMP

                    ENDDO
                   ENDDO
                  ENDDO

GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                  NORM_PW_AUG(NG) = NORM_PW_AUG(NG) + NORM_PW_AUG_TMP
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                  NORM_AUG_PW(NI) = NORM_AUG_PW(NI) + NORM_AUG_PW_TMP

                  WEIGHT_AUG_PW(NG,NI) = WEIGHT_AUG_PW(NG,NI) + WEIGHT_AUG_PW_TMP

               ENDDO
               ENDDO
            ENDDO
GPACC !$ACC EXIT DATA COPYOUT(NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:NORM_PW_AUG,NORM_AUG_PW,WEIGHT_AUG_PW) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
            CALLMPI( M_sum_d(GRIDC%COMM, NORM_AUG_PW(1), T_INFO%NIONS) )

GPACC !$ACC ENTER DATA CREATE(NORM_AUG_AUG,WEIGHT_AUG_AUG) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:NORM_AUG_AUG,WEIGHT_AUG_AUG) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(NORM_AUG_AUG,WEIGHT_AUG_AUG) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:NORM_AUG_AUG,WEIGHT_AUG_AUG) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            NORM_AUG_AUG = 0._q
            WEIGHT_AUG_AUG = 0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(TSMBJ2,M1,M2,M3) &
GPACC !$ACC PRIVATE(NORM_AUG_AUG_TMP,WEIGHT_AUG_AUG_TMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:TSMBJ2,M1,M2,M3) &
GPOMP !$OMP PRIVATE(NORM_AUG_AUG_TMP,WEIGHT_AUG_AUG_TMP) IF(OFFLOAD_ON)
            DO NI=1,T_INFO%NIONS
               DO NJ=1,T_INFO%NIONS

               NORM_AUG_AUG_TMP = 0._q
               WEIGHT_AUG_AUG_TMP = 0._q

GPACC !$ACC LOOP COLLAPSE(3) REDUCTION(+:NORM_AUG_AUG_TMP,WEIGHT_AUG_AUG_TMP) &
GPACC !$ACC PRIVATE(X1,X2,X3,X,Y,Z,D2,WTMP)
GPOMP !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:NORM_AUG_AUG_TMP,WEIGHT_AUG_AUG_TMP) &
GPOMP !$OMP PRIVATE(X1,X2,X3,X,Y,Z,D2,WTMP)
                DO I1=-M1,M1
                 DO I2=-M2,M2
                  DO I3=-M3,M3
                   X1=I1+MODULO(T_INFO%POSION(1,NI)-T_INFO%POSION(1,NJ)+10.5_q,1._q)-0.5_q
                   X2=I2+MODULO(T_INFO%POSION(2,NI)-T_INFO%POSION(2,NJ)+10.5_q,1._q)-0.5_q
                   X3=I3+MODULO(T_INFO%POSION(3,NI)-T_INFO%POSION(3,NJ)+10.5_q,1._q)-0.5_q
       
                   X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                   Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                   Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
       
                   D2=X*X+Y*Y+Z*Z

                   WTMP = EXP(-D2/TSMBJ2)

                   NORM_AUG_AUG_TMP = NORM_AUG_AUG_TMP + WTMP
                   WEIGHT_AUG_AUG_TMP = WEIGHT_AUG_AUG_TMP + WTMP

                  ENDDO
                 ENDDO
                ENDDO

GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                NORM_AUG_AUG(NI) = NORM_AUG_AUG(NI) + NORM_AUG_AUG_TMP

                WEIGHT_AUG_AUG(NI,NJ) = WEIGHT_AUG_AUG(NI,NJ) + WEIGHT_AUG_AUG_TMP

               ENDDO
            ENDDO
GPACC !$ACC EXIT DATA COPYOUT(NORM_AUG_AUG,WEIGHT_AUG_AUG,TSMBJ2,M1,M2,M3) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:NORM_AUG_AUG,WEIGHT_AUG_AUG,TSMBJ2,M1,M2,M3) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

            DO I=1,GRIDC%RL%NP
            ! an initial value maybe given on INCAR (as a "better starting point" ...)
!               CMBJ_AUX(I) = 1.0_q
               CMBJ_AUX(I) = CMBJ
            ENDDO

GPACC !$ACC EXIT DATA DELETE(T_INFO,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:T_INFO,LATT_CUR) IF(OFFLOAD_ON)

#ifdef OFFLOADING
            CALL DELETE_TYPED_VAR(GRIDC)
            CALL OFFLOAD_ON_POP
#endif

         ENDIF
      ENDIF
GPACC !$ACC UPDATE DEVICE(CMBJ_AUX)
GPOMP !$OMP TARGET UPDATE TO(CMBJ_AUX)

      RETURN
      END SUBROUTINE CREATE_CMBJ_AUX

!***********************************************************************
! TR: output of CMBJ_AUX
!***********************************************************************

      SUBROUTINE PRINT_CMBJ_AUX(GRIDC,LATT_CUR,XC)
      USE poscar
      USE mgrid
      USE lattice
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      TYPE (xc_info)     XC
      INTEGER NG, NC, N1, N2, N3, NX, NY, NZ 
      REAL(q) X, Y, Z, X1, X2, X3 
      REAL(q), ALLOCATABLE :: CMBJ_AUX_SER(:,:,:)

      IF (XC%ID(1)==ID_XC_LMBJ) THEN
         ! accumulate the complete CMBJ_AUX into CMBJ_AUX_SER
         ALLOCATE(CMBJ_AUX_SER(GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ))
         CALL MRG_GRID_RL(GRIDC,CMBJ_AUX_SER,CMBJ_AUX)

         OPEN(99,FILE='CMBJCAR')
         ! loop over all grid points
         DO NX=1,GRIDC%NGX
         DO NY=1,GRIDC%NGY
         DO NZ=1,GRIDC%NGZ

            ! only every second point plotted
            IF (MOD(NX,2)==0.AND.MOD(NY,2)==0.AND.MOD(NZ,2)==0) THEN
               X1=(REAL(NX,q)-1)*(1._q/GRIDC%NGX)
               X2=(REAL(NY,q)-1)*(1._q/GRIDC%NGY)
               X3=(REAL(NZ,q)-1)*(1._q/GRIDC%NGZ)

               X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
               Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
               Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
            
               WRITE(99,'(4F14.7)') X, Y, Z, CMBJ_AUX_SER(NX,NY,NZ)
            ENDIF
               
         ENDDO
         ENDDO
         ENDDO
         CLOSE(99)

         DEALLOCATE(CMBJ_AUX_SER)
      ENDIF

      RETURN
      END SUBROUTINE PRINT_CMBJ_AUX

!***********************************************************************
! TR: Create CMBJ_AUG and GRHO_OVER_RHO_AUG
!***********************************************************************

      SUBROUTINE CREATE_CMBJ_AUG(T_INFO,XC)
      USE poscar
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
      TYPE (type_info)   T_INFO
      TYPE (xc_info)     XC
      INTEGER NI

      IF (XC%ID(1)==ID_XC_LMBJ) THEN
         IF (ALLOCATED(CMBJ_AUG)) DEALLOCATE(CMBJ_AUG)
         ALLOCATE(CMBJ_AUG(T_INFO%NIONS))
         IF (ALLOCATED(GRHO_OVER_RHO_AUG)) DEALLOCATE(GRHO_OVER_RHO_AUG)
         ALLOCATE(GRHO_OVER_RHO_AUG(T_INFO%NIONS))
         DO NI=1,T_INFO%NIONS
         ! an initial value maybe given on INCAR (as a "better starting point" ...)
!            CMBJ_AUG(NI) = 1.0_q
            CMBJ_AUG(NI) = CMBJ
            GRHO_OVER_RHO_AUG(NI) = 0.0_q
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CREATE_CMBJ_AUG
 

!***********************************************************************
!
!***********************************************************************
      SUBROUTINE UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,XC,IU6)

__GPU USE moffload

      USE poscar
      USE mgrid
      USE lattice
      USE vaspxml
      USE constant
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (latt)        LATT_CUR
      TYPE (xc_info)     XC

      INTEGER IU6
      ! local variables
      REAL(q) GRHO_OVER_RHO
      REAL(q) CMBJ_AUX_AVG ! average of CMBJ_AUX over all grid points
      INTEGER I,J, N, NC, N1, N2, N3, NX, NY, NZ, NGX, NGY, NGZ
      INTEGER NT
      REAL(q) X, Y, Z
      REAL(q) CMBJ_AVG
      REAL(q) VOL_AUG(T_INFO%NIONS) ! volume of the aug. sphere for PAW-contribution to CMBJ
      REAL(q), ALLOCATABLE :: CMBJ_AUX_H(:)
      REAL(q), ALLOCATABLE :: GRHO_OVER_RHO_H(:)
      ! write(*,*) "starting UPDATE_CMBJ, OFFLOAD_ON=",OFFLOAD_ON
 
      IF ((XC%ID(1)==ID_XC_MBJ).AND.LscMBJ()) THEN
         ! write current CMBJ to OUTCAR and vasprun.xml
         IF (IU6>=0) THEN
GPACC !$ACC UPDATE SELF(CMBJ) IF(OFFLOAD_ON) WAIT
GPOMP !$OMP TARGET UPDATE FROM(CMBJ) IF(OFFLOAD_ON)
            WRITE(IU6,'(4X,A,F10.4/)') 'CMBJ =',CMBJ
            CALL XML_TAG_REAL("CMBJ", CMBJ)
         ENDIF
         ! communicate \nabla rho / rho
GPACC !$ACC UPDATE SELF(GRHO_OVER_RHO_PW) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(GRHO_OVER_RHO_PW) IF(OFFLOAD_ON)
__GPU    CALL OFFLOAD_ON_PUSH(.FALSE.)
         CALLMPI( M_sum_d(GRIDC%COMM, GRHO_OVER_RHO_PW, 1) )
         ! calculate CMBJ
         GRHO_OVER_RHO=GRHO_OVER_RHO_PW/GRIDC%NPLWV+GRHO_OVER_RHO_ONE_CENTER/LATT_CUR%OMEGA
         CMBJ=CMBJA+CMBJB*GRHO_OVER_RHO**CMBJE
         ! just to be sure
         IF (CMBJ<1._q) CMBJ=1._q
         IF (CMBJ>2._q) CMBJ=2._q
         ! initialize for the next step
         GRHO_OVER_RHO_PW=0._q; GRHO_OVER_RHO_ONE_CENTER=0._q
__GPU    CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(GRHO_OVER_RHO_PW,CMBJ) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(GRHO_OVER_RHO_PW,CMBJ) IF(OFFLOAD_ON)
      ENDIF
     
      ! TR: for local MBJ 
      !     calculate the averages of the augmentation spheres
      !     and add the missing parts (PW and AUG) to respective CMBJ
      IF (XC%ID(1)==ID_XC_LMBJ) THEN

         ! merge averages inside augmentation spheres
         CALLMPI( M_sum_d(GRIDC%COMM, GRHO_OVER_RHO_AUG(1), T_INFO%NIONS) )
         CALLMPI( M_sum_d(GRIDC%COMM_KINTER, GRHO_OVER_RHO_AUG(1), T_INFO%NIONS) )

         ! smear AUG contributions
         ALLOCATE(GRHO_OVER_RHO_H(T_INFO%NIONS))

         DO I=1,T_INFO%NIONS
            GRHO_OVER_RHO_H(I) = 0._q
            DO J=1,T_INFO%NIONS
               GRHO_OVER_RHO_H(I) = GRHO_OVER_RHO_H(I) + GRHO_OVER_RHO_AUG(J)*WEIGHT_AUG_AUG(I,J)/NORM_AUG_AUG(I)
            ENDDO
         ENDDO

         ! add the AUG contribution to PW
         ALLOCATE(CMBJ_AUX_H(GRIDC%RL%NP))

         ! volume of the aug. sphere for AUG contribution
         VOL_AUG = ((1.0_q-VACRAT)*LATT_CUR%OMEGA/T_INFO%NIONS)

         ! To allow the use of CMBJE with LMBJ, the arrays CMBJ_AUX, CMBJ_AUX_H
         ! and CMBJ_AUX_AVG are temporarily storing grad(rho)/rho instead of c.
         DO I=1,GRIDC%RL%NP
            CMBJ_AUX(I) = (CMBJ_AUX(I) - CMBJA)/(SIGN(1.0_q,CMBJB)*MAX(ABS(CMBJB),1.E-10_q))
            CMBJ_AUX_H(I) = CMBJ_AUX(I)
            DO J=1,T_INFO%NIONS
               ! Add AUG to PW 
               CMBJ_AUX_H(I) = CMBJ_AUX_H(I) + SCA_VAC(I)*GRHO_OVER_RHO_AUG(J)/VOL_AUG(J)*WEIGHT_AUG_PW(I,J)/NORM_PW_AUG(I)
            ENDDO
         ENDDO

         DO I=1,T_INFO%NIONS
            GRHO_OVER_RHO_AUG(I) = GRHO_OVER_RHO_H(I)
         ENDDO
         DEALLOCATE(GRHO_OVER_RHO_H)

         ! average of CMBJ_AUX and add it to AUG
         DO I=1,T_INFO%NIONS
            CMBJ_AUX_AVG = 0._q
            DO J=1,GRIDC%RL%NP
               CMBJ_AUX_AVG = CMBJ_AUX_AVG + CMBJ_AUX(J)*WEIGHT_AUG_PW(J,I)/NORM_AUG_PW(I)
            ENDDO
            CALLMPI( M_sum_d(GRIDC%COMM, CMBJ_AUX_AVG, 1) )
            CMBJ_AUG(I) = CMBJA + CMBJB*(CMBJ_AUX_AVG+GRHO_OVER_RHO_AUG(I)/VOL_AUG(I))**CMBJE
        ENDDO

         ! set CMBJ_AUX to the full value
         DO I=1,GRIDC%RL%NP
            CMBJ_AUX(I) = CMBJA + CMBJB*CMBJ_AUX_H(I)**CMBJE
         ENDDO
         DEALLOCATE(CMBJ_AUX_H)

      ! find the average CMBJ
         CMBJ_AVG = 0.0_q

         DO I=1,GRIDC%RL%NP
            CMBJ_AVG = CMBJ_AVG + CMBJ_AUX(I)
         ENDDO 

         CALLMPI( M_sum_d(GRIDC%COMM, CMBJ_AVG, 1) )

         CMBJ_AVG = CMBJ_AVG / GRIDC%NPLWV
         IF (IU6>=0) THEN
            WRITE(IU6,'(4X,A,F10.4/)') 'CMBJ_AVG_PW =',CMBJ_AVG
            CALL XML_TAG_REAL("CMBJ_AVG_PW", CMBJ_AVG)
            WRITE(IU6,'(4X,A,F10.4/)') 'VACRAT =',VACRAT
            CALL XML_TAG_REAL("VACRAT", VACRAT)
            DO I=1,T_INFO%NIONS
               WRITE(IU6,'(4X,I3,X,A,F10.4/)') I,'CMBJ_AUG(I) =',CMBJ_AUG(I)
               CALL XML_TAG_REAL("CMBJ_AUG(I)", CMBJ_AUG(I))
            ENDDO
         ENDIF

         GRHO_OVER_RHO_AUG = 0._q

GPACC !$ACC UPDATE DEVICE(CMBJ_AUX) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(CMBJ_AUX) IF(OFFLOAD_ON)

      ENDIF

      RETURN
      END SUBROUTINE UPDATE_CMBJ

!***********************************************************************
! TR: update CMBJ_AUX, for local SC CMBJ
!***********************************************************************

      SUBROUTINE UPDATE_CMBJ_AUX(GRIDC,LATT_CUR,XC,RHO,DENCOR,DRHO1,DRHO2,DRHO3,NCDIJ)

__GPU USE moffload_struct_def

      USE poscar
      USE prec
      USE mgrid
      USE lattice
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      TYPE (xc_info)     XC
      INTEGER NCDIJ
      RGRID RHO(DIMREAL(GRIDC%MPLWV),NCDIJ)
      RGRID DENCOR(GRIDC%RL%NP)
      REAL(q) DRHO1(GRIDC%RL%NP,NCDIJ),DRHO2(GRIDC%RL%NP,NCDIJ),DRHO3(GRIDC%RL%NP,NCDIJ)
      
      COMPLEX(q) GRHO_OVER_RHO_LOC(GRIDC%MPLWV), GAUSSIAN(GRIDC%MPLWV), RES(GRIDC%MPLWV)
      RGRID GRHO_OVER_RHO_LOC_R(GRIDC%RL%NP), RES_R(GRIDC%RL%NP)

      INTEGER NC, N1, N2, N3
      INTEGER NX, NY, NZ, NG, IS
      INTEGER NGX, NGY, NGZ
      REAL(q) X, Y, Z 
      REAL(q) X1, X2, X3 
      REAL(q) GX, GY, GZ
      INTEGER I
      REAL(q) RHOUP, RHODW, MAG_NORM
      REAL(q) ABSNABUP,ABSNABDW

      REAL(q), PARAMETER :: THR_HALF = 3.0_q/2.0_q
      REAL(q) :: PREF 
      REAL(q) :: CHGMIN=1E-10_q
      REAL(q) RSQUARE,SMBJ2H,RSMBJ3,TPI2
      REAL(q) RINPL
      REAL(q) TMBJ

      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) THEN
         IF (CMBJ_TYP(1)/=-1) RETURN
         IF (XC%ID(1)==ID_XC_LMBJ) THEN
         IF (ALLOCATED(CMBJ_AUX)) THEN
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CMBJ_AUX)
            DEALLOCATE(CMBJ_AUX)
         ENDIF
         ! write(*,*) 'allocating CMBJ_AUX'
         ALLOCATE(CMBJ_AUX(GRIDC%RL%NP))
! Allocatable arrays need to be mapped to the device,
! even if they are declared as TARGET in OpenMP offloading.
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:CMBJ_AUX)
#ifdef CRAYHIP
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CMBJ_AUX)
#endif
        
         RSMBJ3=RSMBJ*RSMBJ*RSMBJ*AUTOA3
         TMBJ = (8.0_q/3.0_q)*PI*RSMBJ3
         SMBJ2H=SMBJ*SMBJ*0.5_q
         TPI2=TPI*TPI

         ! loop over all grid points
         NGX = GRIDC%NGX
         NGY = GRIDC%NGY
         NGZ = GRIDC%NGZ

         ! calculate the VACRAT parameter
         VACRAT = 0.0_q
GPACC !$ACC PARALLEL LOOP PRESENT(RHO,DENCOR,LATT_CUR) PRIVATE(RHOUP,RHODW) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:RHO,DENCOR,LATT_CUR) PRIVATE(RHOUP,RHODW) IF(OFFLOAD_ON)
         DO I=1,GRIDC%RL%NP
            IF ((NCDIJ==1).OR.(NCDIJ==4)) THEN
               RHOUP= MAX(REAL((RHO(I,1)+DENCOR(I))*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               RHODW= RHOUP
            ELSEIF (NCDIJ==2) THEN
               RHOUP= MAX(REAL((RHO(I,1)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               RHODW= MAX(REAL((RHO(I,2)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            ENDIF
            IF ((RHOUP+RHODW)*0.5_q < 1._q/TMBJ) THEN
               VACRAT = VACRAT + 1.0_q
            ENDIF
         ENDDO
         CALLMPI( M_sum_d(GRIDC%COMM, VACRAT, 1) )
         VACRAT = VACRAT / GRIDC%NPLWV

GPACC !$ACC ENTER DATA CREATE(GRHO_OVER_RHO_LOC_R,GRHO_OVER_RHO_LOC,RES,RES_R) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:GRHO_OVER_RHO_LOC_R,GRHO_OVER_RHO_LOC,RES,RES_R) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL PRESENT(RHO,DENCOR,DRHO1,DRHO2,DRHO3) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(RHOUP,RHODW,ABSNABUP,ABSNABDW,MAG_NORM)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:RHO,DENCOR,DRHO1,DRHO2,DRHO3) &
GPOMP !$OMP PRIVATE(RHOUP,RHODW,ABSNABUP,ABSNABDW,MAG_NORM) IF(OFFLOAD_ON)

         DO I=1,GRIDC%RL%NP
            ! First: Calculate GRHO_OVER_RHO_LOC

            IF (NCDIJ==1) THEN
               ! nonspin-polarized case
               RHOUP= MAX(REAL((RHO(I,1)+DENCOR(I))*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               RHODW= RHOUP

               ABSNABUP= 0.5_q*SQRT(DRHO1(I,1)*DRHO1(I,1)+DRHO2(I,1)*DRHO2(I,1) &
                       +DRHO3(I,1)*DRHO3(I,1))
               ABSNABDW= ABSNABUP
            ELSEIF (NCDIJ==2) THEN
               ! Spin-polarized case
               RHOUP= MAX(REAL((RHO(I,1)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               RHODW= MAX(REAL((RHO(I,2)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
  
               ABSNABUP=SQRT(DRHO1(I,1)*DRHO1(I,1)+DRHO2(I,1)*DRHO2(I,1) &
                          +DRHO3(I,1)*DRHO3(I,1))
 
               ABSNABDW=SQRT(DRHO1(I,2)*DRHO1(I,2)+DRHO2(I,2)*DRHO2(I,2) &
                          +DRHO3(I,2)*DRHO3(I,2))
            ELSEIF (NCDIJ==4) THEN
               ! noncollinear case
               !     
               ! | m |
               !
               MAG_NORM=SQRT(ABS( &
               &    RHO(I,2)*RHO(I,2)+RHO(I,3)*RHO(I,3)+RHO(I,4)*RHO(I,4)))
               !       
               ! RHOUP: \rho_up   = ( \rho_tot + | m | )/2
               ! RHODW: \rho_down = ( \rho_tot - | m | )/2
               !
               RHOUP= MAX(REAL((RHO(I,1)+DENCOR(I)+MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               RHODW= MAX(REAL((RHO(I,1)+DENCOR(I)-MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
               ! \nabla \rho already calculated in upper routine
               ABSNABUP=SQRT(DRHO1(I,1)*DRHO1(I,1)+DRHO2(I,1)*DRHO2(I,1) &
                          +DRHO3(I,1)*DRHO3(I,1))

               ABSNABDW=SQRT(DRHO1(I,2)*DRHO1(I,2)+DRHO2(I,2)*DRHO2(I,2) &
                          +DRHO3(I,2)*DRHO3(I,2))
            ENDIF

            IF (VACRAT>0.0_q) THEN
               ! f = erf(RHO*TMBJ)
               GRHO_OVER_RHO_LOC_R(I) = ((1._q-CMBJA)/CMBJB + (AUTOA*ABSNABUP/RHOUP - (1._q-CMBJA)/CMBJB)*ERF(RHOUP*TMBJ))*0.5_q + &
              &                         ((1._q-CMBJA)/CMBJB + (AUTOA*ABSNABDW/RHODW - (1._q-CMBJA)/CMBJB)*ERF(RHODW*TMBJ))*0.5_q 
               ! factor for vacuum
               SCA_VAC(I) = ERF((RHOUP+RHODW)*0.5_q*TMBJ)
            ELSE
               ! f = 1
               GRHO_OVER_RHO_LOC_R(I) = (AUTOA*ABSNABUP/RHOUP + AUTOA*ABSNABDW/RHODW)*0.5_q
               SCA_VAC(I) = 1.0_q
            ENDIF

         ENDDO
GPACC !$ACC END PARALLEL

         ! FFT to reciprocal space
         RINPL=1._q/GRIDC%NPLWV
         ! GRHO_OVER_RHO_LOC         
         CALL RL_ADD(GRHO_OVER_RHO_LOC_R,RINPL,GRHO_OVER_RHO_LOC_R,0._q,GRHO_OVER_RHO_LOC,GRIDC)
         CALL FFT3D(GRHO_OVER_RHO_LOC(1),GRIDC,-1)
         CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, GRHO_OVER_RHO_LOC(1))

GPACC !$ACC UPDATE SELF(GRHO_OVER_RHO_LOC) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET UPDATE FROM(GRHO_OVER_RHO_LOC) IF(OFFLOAD_ON)

         ! Multiply in RC
         DO I=1,GRIDC%RC%NP
            ! Calculate Gaussian
            N1= MOD((I-1),GRIDC%RC%NROW) +1
            NC= (I-1)/GRIDC%RC%NROW+1
            N2= GRIDC%RC%I2(NC)
            N3= GRIDC%RC%I3(NC)
            GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
            GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
            GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

            RSQUARE=(GX*GX+GY*GY+GZ*GZ)*TPI2
            GAUSSIAN(I) = EXP(-SMBJ2H*RSQUARE)
            RES(I) = GRHO_OVER_RHO_LOC(I)*GAUSSIAN(I)
         ENDDO

GPACC !$ACC UPDATE DEVICE(RES) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(RES) IF(OFFLOAD_ON)

         ! FFT RES back to real space
         CALL SETUNB(RES(1),GRIDC)
         CALL FFT3D(RES(1),GRIDC,1)
         CALL RL_ADD(RES,1.0_q,RES,0._q,RES_R,GRIDC)

GPACC !$ACC UPDATE SELF(RES_R) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET UPDATE FROM(RES_R) IF(OFFLOAD_ON)

         ! Calculate CMBJ_AUX
         DO I=1,GRIDC%RL%NP
            CMBJ_AUX(I) = CMBJA + CMBJB*RES_R(I) 
         ENDDO

GPACC !$ACC EXIT DATA DELETE(GRHO_OVER_RHO_LOC_R,GRHO_OVER_RHO_LOC,RES,RES_R) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:GRHO_OVER_RHO_LOC_R,GRHO_OVER_RHO_LOC,RES,RES_R) IF(OFFLOAD_ON)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE UPDATE_CMBJ_AUX


!***********************************************************************
!
!***********************************************************************
      SUBROUTINE DESTROY_CMBJ_AUX()
      IMPLICIT NONE
      IF (ALLOCATED(CMBJ_AUX)) DEALLOCATE(CMBJ_AUX)
      RETURN
      END SUBROUTINE DESTROY_CMBJ_AUX

  END MODULE mbj
