#include "symbol.inc"

  MODULE tau_mu
    USE prec

    TYPE tau_handle
       COMPLEX(q), POINTER, CONTIGUOUS :: TAU(:,:)  => NULL() ! kinetic energy density
       COMPLEX(q), POINTER, CONTIGUOUS :: TAUL(:,:) => NULL() ! kinetic energy density (previous step)
       COMPLEX(q), POINTER, CONTIGUOUS :: TAUW(:,:) => NULL() ! Weisaecker kinetic energy density
       RGRID,      POINTER, CONTIGUOUS :: TAUC(:)   => NULL() ! kinetic energy density of the core electrons
    END TYPE tau_handle

   CONTAINS


      SUBROUTINE MGGA_READER(INCAR_F,XC_DATA,IO)
        USE incar_reader, ONLY: incar_file, process_incar
        USE base, ONLY: in_struct
        USE setexm_struct_def, only : xc_info
        USE tutor, ONLY: vtutor
        IMPLICIT NONE

        TYPE(incar_file), INTENT(INOUT)  :: INCAR_F
        TYPE(xc_info),    INTENT(INOUT)  :: XC_DATA
        TYPE (in_struct), INTENT(INOUT)  :: IO

        ! write TAUCAR or tau to vaspwave.h5
        ! default: if MGGA T, else F
        IO%WRT_DENSITY%LTAU = XC_DATA%LTAU
        CALL PROCESS_INCAR(INCAR_F, "LTAU", IO%WRT_DENSITY%LTAU)
      END SUBROUTINE MGGA_READER

!***********************************************************************
!
!***********************************************************************

      SUBROUTINE ALLOCATE_MU(LMU,MU,MUTOT,GRID,GRIDC,WDES)
      USE mgrid
      USE ini
      USE wave
      IMPLICIT NONE
      TYPE (grid_3d) GRID, GRIDC
      TYPE (wavedes) WDES
      COMPLEX(q),POINTER :: MUTOT(:,:)  ! d e_xc / d tau
      RGRID , POINTER    :: MU(:,:)     ! soft version of the above
      LOGICAL :: LMU

      IF (LMU) THEN
         ALLOCATE(MU(DIMREAL(GRID%MPLWV),WDES%NCDIJ),MUTOT(GRIDC%MPLWV,WDES%NCDIJ))
         CALL REGISTER_ALLOCATE(16._q*(SIZE(MUTOT)+SIZE(MU)), "grid")
         MU=0; MUTOT=0
      ELSE
         NULLIFY(MUTOT, MU)
      ENDIF

      RETURN
      END SUBROUTINE ALLOCATE_MU


!***********************************************************************
!
!***********************************************************************

      SUBROUTINE DEALLOCATE_MU(LMU,MU,MUTOT)
      USE mgrid
      USE ini
      IMPLICIT NONE
      COMPLEX(q),POINTER :: MUTOT(:,:)  ! d e_xc / d tau
      RGRID , POINTER    :: MU(:,:)     ! soft version of the above
      LOGICAL :: LMU

      IF (LMU) THEN
         CALL DEREGISTER_ALLOCATE(16._q*(SIZE(MUTOT)+SIZE(MU)), "grid")
         DEALLOCATE(MU,MUTOT)
      ENDIF
      NULLIFY(MUTOT, MU)

      RETURN
      END SUBROUTINE DEALLOCATE_MU


!!!***********************************************************************
!!!
!!! generate the tau_handle
!!!
!!!***********************************************************************
!!
!!      SUBROUTINE GENERATE_TAU_HANDLE(LTAU, KINEDEN, GRIDC, NCDIJ)
!!      USE mgrid
!!      IMPLICIT NONE
!!      TYPE(tau_handle), POINTER :: KINEDEN
!!      TYPE(grid_3d) GRIDC
!!      INTEGER NCDIJ
!!      LOGICAL LTAU
!!
!!      IF (.NOT.ASSOCIATED(KINEDEN).AND.LTAU) THEN
!!         ALLOCATE(KINEDEN)
!!         ALLOCATE(KINEDEN%TAU(GRIDC%MPLWV,NCDIJ),KINEDEN%TAUC(GRIDC%RL%NP))
!!         KINEDEN%TAU=0; KINEDEN%TAUC=0
!!         IF (LMIX_TAU()) THEN
!!            ALLOCATE(KINEDEN%TAUL(GRIDC%MPLWV,NCDIJ))
!!            KINEDEN%TAUL=0
!!         ELSE
!!            NULLIFY(KINEDEN%TAUL)
!!         ENDIF
!!      ENDIF
!!
!!      RETURN
!!      END SUBROUTINE GENERATE_TAU_HANDLE


!!!***********************************************************************
!!!
!!! deallocate the tau_handle
!!!
!!!***********************************************************************
!!
!!      SUBROUTINE DEALLOCATE_TAU_HANDLE(KINEDEN)
!!      IMPLICIT NONE
!!      TYPE(tau_handle), POINTER :: KINEDEN
!!
!!      IF (ASSOCIATED(KINEDEN)) THEN
!!         IF (ASSOCIATED(KINEDEN%TAU)) DEALLOCATE(KINEDEN%TAU)
!!         IF (ASSOCIATED(KINEDEN%TAUL)) DEALLOCATE(KINEDEN%TAUL)
!!         IF (ASSOCIATED(KINEDEN%TAUC)) DEALLOCATE(KINEDEN%TAUC)
!!         DEALLOCATE(KINEDEN)
!!      ENDIF
!!
!!      RETURN
!!      END SUBROUTINE DEALLOCATE_TAU_HANDLE


!***********************************************************************
!
! generate the tau_handle
!
!***********************************************************************

      SUBROUTINE ALLOCATE_TAU(LTAU, TAU, TAUC, TAUL, GRIDC, NCDIJ)
      USE mgrid
      IMPLICIT NONE
      COMPLEX(q), POINTER :: TAU(:,:),TAUL(:,:)
      RGRID, POINTER      :: TAUC(:)
      TYPE(grid_3d) GRIDC
      LOGICAL :: LTAU
      INTEGER NCDIJ

      IF (LTAU) THEN
         ALLOCATE(TAU(GRIDC%MPLWV,NCDIJ),TAUC(GRIDC%RL%NP))
         TAU=0; TAUC=0
         IF (LMIX_TAU()) THEN
            ALLOCATE(TAUL(GRIDC%MPLWV,NCDIJ))
            TAUL=0
         ELSE
            NULLIFY(TAUL)
         ENDIF
      ELSE
         NULLIFY(TAU,TAUC,TAUL)
      ENDIF

      RETURN
      END SUBROUTINE ALLOCATE_TAU


!***********************************************************************
!
! deallocate the tau_handle
!
!***********************************************************************

      SUBROUTINE DEALLOCATE_TAU(LTAU,TAU,TAUC,TAUL)
      IMPLICIT NONE
      COMPLEX(q), POINTER :: TAU(:,:),TAUL(:,:)
      RGRID, POINTER      :: TAUC(:)
      LOGICAL :: LTAU

      IF (LTAU) THEN
         DEALLOCATE(TAU,TAUC)
         IF (LMIX_TAU()) DEALLOCATE(TAUL)
      ENDIF
      NULLIFY(TAU,TAUC,TAUL)

      RETURN
      END SUBROUTINE DEALLOCATE_TAU


!***********************************************************************
!
!***********************************************************************
      FUNCTION LMIX_TAU()
      USE setexm_struct_def, ONLY : LMIXTAU
      IMPLICIT NONE
      LOGICAL LMIX_TAU
      LMIX_TAU=LMIXTAU
      END FUNCTION LMIX_TAU


!******************** SUBROUTINE SET_KINEDEN ***************************
!
! SET_KINEDEN calculates the kinetic energy density of the PW part of
! the wavefunctions (0.5*|grad psi|**2).
!
! The actual work is done in TAU_PW
!
!***********************************************************************

      SUBROUTINE SET_KINEDEN(&
     &   GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM,LTAU,NIOND,W,WDES,KINEDEN &
     &)

__GPU USE moffload

      USE prec
      USE lattice
      USE mgrid
      USE msymmetry
      USE base
      USE wave
      USE mpimy
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC,GRID,GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (symmetry)    SYMM
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      TYPE (tau_handle)  KINEDEN

      INTEGER NIOND
      LOGICAL LTAU

      ! quick return if possible
      IF (.NOT.LTAU) THEN
         RETURN
      ENDIF
      PROFILING_START('set_kineden')

      ! just a failsafe
      IF (.NOT.ASSOCIATED(KINEDEN%TAU)) THEN
         CALL vtutor%bug("SET_KINEDEN: internal ERROR, KIN_EDEN handle not allocated", __FILE__, __LINE__)
      ENDIF

      CALL TAU_PW(&
     &   GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM,NIOND,W,WDES,KINEDEN%TAU)

GPACC !$ACC UPDATE SELF(KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON)
      PROFILING_STOP('set_kineden')

      RETURN
      END SUBROUTINE SET_KINEDEN


!************************* SUBROUTINE TAUPAR ***************************
!
!  subroutine to calculate the partial kinetic energy density of the
!  core electrons from overlapping atoms in real space and store
!  it in an real array
!  NOTE: the DENS array is not necessarily suffiently large
!  to allow for 3d-FFT, hence a work array is temporarily allocated
!
!***********************************************************************

      SUBROUTINE TAUPAR(GRIDC,T_INFO,B,OMEGA,P,CSTRF,DENS)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      RGRID        DENS(GRIDC%RL%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3),OMEGA
      ! work arrays
      COMPLEX(q),ALLOCATABLE::    CWORK1(:),CWORK2(:)

      PROFILING_START('taupar')
      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: CWORK1,CWORK2) IF(OFFLOAD_ON)
      CALL TAUATO(.TRUE.,GRIDC,T_INFO,B,P,CSTRF,CWORK1,CWORK2)
      CALL RL_ADD(CWORK1,1.0_q/OMEGA,CWORK1,0.0_q,DENS,GRIDC)

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: CWORK1,CWORK2) IF(OFFLOAD_ON)
      DEALLOCATE(CWORK1,CWORK2)
      PROFILING_STOP('taupar')

      RETURN
      END SUBROUTINE TAUPAR


!************************* SUBROUTINE TAUATO ***************************
!
!  This routine calculates the partial kinetic energy density of the
!  core electrons corresponding to overlapping atoms
!
!  LFOUR
!   .FALSE. set up density in reciprocal space
!   .TRUE.  set up density in real space
!
!***********************************************************************

      SUBROUTINE TAUATO(LFOUR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CHDER)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC

      INTEGER :: N, N1, N2, N3, NADDR, NC, NT

      COMPLEX(q)   CHTOT(GRIDC%RC%NP),CHDER(GRIDC%RC%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

      REAL(q)      B(3,3), ARG, ARGSC, G, GX, GY, GZ, PSGMA2, REM, T0, T1, T2, T3, V1, V2, V3, V4

      LOGICAL LFOUR
      ! local variables
      REAL(q), POINTER :: PRHO(:)

      PROFILING_START('tauato')
GPACC !$ACC KERNELS PRESENT(CHTOT,CHDER) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHTOT,CHDER) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CHTOT=0
      CHDER=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
       PRHO=>P(NT)%PSPTAU
       IF (.NOT.ASSOCIATED(PRHO)) CYCLE
!=======================================================================
! calculate the scaling factor ARGSC that converts the magnitude of a
! reciprocal lattice vector to the correponding position in the
! pseudopotential arrays
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

GPACC !$ACC PARALLEL LOOP PRESENT(P,GRIDC,CHTOT,CSTRF,CHDER) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(N1,N2,N3,NC,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:P,GRIDC,CHTOT,CSTRF,CHDER) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3) IF(OFFLOAD_ON) 
      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX*GX+GY*GY+GZ*GZ)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          CHTOT(N)=CHTOT(N)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N,NT)
          CHDER(N)=CHDER(N)+(T1+REM*(2*T2+REM*3*T3))*ARGSC*CSTRF(N,NT)
        ELSE IF (G==0) THEN
          CHTOT(N)=CHTOT(N)+PRHO(1)*CSTRF(N,NT)
          CHDER(N)=0
        ENDIF
      ENDDO
      ENDDO typ
!=======================================================================
! set the charge-density of unbalanced lattice-vectors to 0
! and transform the charge-density to real space
!=======================================================================
      CALL SETUNB(CHTOT(1),GRIDC)
      CALL SETUNB(CHDER(1),GRIDC)

      IF (LFOUR) THEN
        CALL FFT3D(CHTOT,GRIDC,1)
      ENDIF

      PROFILING_STOP('tauato')

      RETURN
      END SUBROUTINE TAUATO

!
!  small subroutine which allocates CHDER calls the subroutine
!  above, and then returns
!

      SUBROUTINE TAUATO_WORK(LFOUR,GRIDC,T_INFO,B,P,CSTRF,CHTOT)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%RC%NP)
      COMPLEX(q), ALLOCATABLE :: CHDER(:)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR
      ! local variables
      REAL(q), POINTER :: PRHO(:)

      ALLOCATE(CHDER(GRIDC%RC%NP))
      CALL TAUATO(LFOUR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CHDER)
      DEALLOCATE(CHDER)

      END SUBROUTINE TAUATO_WORK


!************************ SUBROUTINE FORTAU ****************************
!
! This subroutine calculates the correction to the Hellmann-Feynman
! forces due to the (possible) presence of a partial core kinetic energy
! density
!
!***********************************************************************

      SUBROUTINE FORTAU(GRIDC,P,T_INFO,LATT_CUR,MUTOT,TAUFOR)

__GPU USE moffload_struct_def

      USE prec

      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      INTEGER :: N,N1,N1P,N2,N3,NADDR,NC,NG,NGP,NI,NIADD,NIS,NT

      REAL(q)      TAUFOR(3,T_INFO%NIONS)
      ! work arrays
      REAL(q), POINTER :: PRHO(:)
      REAL(q), ALLOCATABLE :: WORK(:)

      REAL(q) ARG,ARGSC,FACTM,FOR,FOR1,FOR2,FOR3,G,G1,G2,G3,GX,GY,GZ,PSGMA2,REM,SCALE,T0,T1,T2,T3,V1,V2,V3,V4

      COMPLEX(q)   MUTOT(GRIDC%MPLWV)

      COMPLEX(q) CE,CEXPF,CX

      PROFILING_START('fortau')

      ALLOCATE(WORK(GRIDC%RC%NP))
GPACC !$ACC ENTER DATA CREATE(WORK,TAUFOR) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: WORK,TAUFOR) IF(OFFLOAD_ON)

      SCALE=1._q/(RYTOEV*AUTOA2)

GPACC !$ACC KERNELS PRESENT(TAUFOR) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:TAUFOR) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      TAUFOR(1:3,1:T_INFO%NIONS)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! loop over all types of atoms
!=======================================================================
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
      NIADD=T_INFO%NITYP(NT)

      PRHO=>P(NT)%PSPTAU

      IF (ASSOCIATED(PRHO)) THEN
!=======================================================================
! iterpolate the pseudopotential on the grid of reciprocal
! lattice-vectors
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

GPACC !$ACC PARALLEL LOOP PRIVATE(N1,NC,N2,N3,FACTM,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3) &
GPACC !$ACC PRESENT(GRIDC,LATT_CUR,P,WORK) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,LATT_CUR,P,WORK) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3,FACTM,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3) IF(OFFLOAD_ON)
      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)

        SETFACT1
        SETFACT
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX*GX+GY*GY+GZ*GZ)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          WORK(N)=T0+REM*(T1+REM*(T2+REM*T3))
!       ELSE IF (G==0) THEN
!         WORK(N)=PRHO(1)
        ELSE
          WORK(N)=0
        ENDIF

      ENDDO

#ifdef vector
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(FOR1,FOR2,FOR3) &
#else
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(FOR1,FOR2,FOR3,CX,G1) &
#endif
GPACC !$ACC PRESENT(GRIDC,T_INFO,T_INFO%POSION,T_INFO%VCA,WORK,TAUFOR,MUTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:GRIDC,T_INFO,T_INFO%POSION,T_INFO%VCA,WORK,TAUFOR,MUTOT) &
GPOMP !$OMP PRIVATE(FOR1,FOR2,FOR3,CX,G1,NC) IF(OFFLOAD_ON)
      ion: DO NI=NIS,NIADD+NIS-1
!=======================================================================
! initialise the force on the ion to zero
!=======================================================================
         FOR1=0._q
         FOR2=0._q
         FOR3=0._q
!=======================================================================
! calculate the total force on the ions by summing over reciprocal
! lattice vectors
! first calculate phase factor:
! there are two version for calculating the phase factor
! on vector machines you might try the first version
! but usually the second version is much faster (less calls to CEXP)
!=======================================================================
#ifdef vector
GPACC !$ACC LOOP VECTOR REDUCTION(-:FOR1,FOR2,FOR3) PRIVATE(N1,NC,N2,N3,G1,G2,G3,FACTM,CEXPF,FOR)
         DO NG=1,GRIDC%RC%NP

           N1= MOD((NG-1),GRIDC%RC%NROW) +1
           NC= (NG-1)/GRIDC%RC%NROW+1
           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)

           G1=T_INFO%POSION(1,NI)*GRIDC%LPCTX(N1)
           G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
           G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
           SETFACT1
           SETFACT
           CEXPF=EXP(-CITPI*(G3+G2+G1))*T_INFO%VCA(NT)
#else
GPOMP !$OMP PARALLEL FIRSTPRIVATE(CX,G1) PRIVATE(NGP,N2,N3,G2,G3,CE,N1,NG,FACTM,CEXPF,FOR) &
GPOMP !$OMP SHARED(GRIDC,T_INFO,NI,NT,WORK,FOR1,FOR2,FOR3) IF(OFFLOAD_ON)
         CX =EXP(-CITPI*T_INFO%POSION(1,NI))
         G1 =T_INFO%POSION(1,NI)*xmin(GRIDC%NGX)

GPACC !$ACC LOOP VECTOR REDUCTION(-:FOR1,FOR2,FOR3) PRIVATE(NGP,N2,N3,G2,G3,CE,N1,NG,FACTM,CEXPF,FOR)
GPOMP !$OMP DO REDUCTION(+:FOR1,FOR2,FOR3)
         DO NC=1,GRIDC%RC%NCOL
           NGP=(NC-1)*GRIDC%RC%NROW+1

           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)
           G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
           G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
           CE=EXP(-CITPI*(G3+G2+G1))*T_INFO%VCA(NT)

GPACC !$ACC LOOP SEQ
           DO N1P=0,GRIDC%RC%NROW-1
           SETN1fromN1P
           NG=NGP+N1
           N1=N1+1
           SETFACT1
           SETFACT
           CEXPF=CE
           CE=CE*CX
#endif
!=======================================================================
! add the contribution to the force from the present reciprocal lattice
! vector  and multiply by i (ie take imaginary part)
!=======================================================================
           FOR=WORK(NG)*MULFACT AIMAG(CONJG(MUTOT(NG))*CEXPF)

           FOR1=FOR1-GRIDC%LPCTX_(N1)*FOR
           FOR2=FOR2-GRIDC%LPCTY_(N2)*FOR
           FOR3=FOR3-GRIDC%LPCTZ_(N3)*FOR
         ENDDO
#ifndef vector
         ENDDO
#endif
GPOMP !$OMP END PARALLEL
!=======================================================================
! multiply forces by 2*Pi
!=======================================================================
         TAUFOR(1,NI)=FOR1*TPI*SCALE
         TAUFOR(2,NI)=FOR2*TPI*SCALE
         TAUFOR(3,NI)=FOR3*TPI*SCALE

      ENDDO ion
      ENDIF
      NIS=NIS+NIADD
      ENDDO typ
GPACC !$ACC EXIT DATA COPYOUT(TAUFOR) DELETE(WORK) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:TAUFOR) MAP(RELEASE: WORK) IF(OFFLOAD_ON)
!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
      CALLMPI( M_sum_d(GRIDC%COMM,TAUFOR(1,1),T_INFO%NIONS*3))
      CALL DIRKAR(T_INFO%NIONS,TAUFOR,LATT_CUR%B)
      DEALLOCATE(WORK)

      PROFILING_STOP('fortau')

      RETURN
      END SUBROUTINE FORTAU
!***********************************************************************
!***********************************************************************
  END MODULE tau_mu
!***********************************************************************
!***********************************************************************

      SUBROUTINE MUxTAU(GRIDC,NCDIJ,LATT_CUR,DMUTOT,DKINEDEN,TAUC,DE)

__GPU USE moffload_struct_def

      USE prec
      USE constant
      USE lattice
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(grid_3d) GRIDC
      TYPE(latt) LATT_CUR
      INTEGER NCDIJ
      RGRID DMUTOT(DIMREAL(GRIDC%MPLWV),NCDIJ)
      RGRID DKINEDEN(DIMREAL(GRIDC%MPLWV),NCDIJ),TAUC(GRIDC%RL%NP)
      REAL(q) DE
      ! local variables
      INTEGER I
      REAL(q) RINPL

      PROFILING_START('mu_x_tau')

      DE=0
      RINPL=1._q/GRIDC%NPLWV

      IF (NCDIJ==1) THEN
GPACC !$ACC PARALLEL LOOP REDUCTION(+:DE) PRESENT(DMUTOT,DKINEDEN,TAUC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO REDUCTION(+:DE) MAP(PRESENT,TO:DMUTOT,DKINEDEN,TAUC) IF(OFFLOAD_ON)
         DO I=1,GRIDC%RL%NP
            DE=DE+DMUTOT(I,1)*CHGCONJG(DKINEDEN(I,1)+TAUC(I))
         ENDDO
      ELSEIF (NCDIJ==2) THEN
GPACC !$ACC PARALLEL LOOP REDUCTION(+:DE) PRESENT(DMUTOT,DKINEDEN,TAUC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO REDUCTION(+:DE) MAP(PRESENT,TO:DMUTOT,DKINEDEN,TAUC) IF(OFFLOAD_ON)
         DO I=1,GRIDC%RL%NP
            DE=DE+DMUTOT(I,1)*CHGCONJG(DKINEDEN(I,1)+TAUC(I)/2._q)+ &
           &       DMUTOT(I,2)*CHGCONJG(DKINEDEN(I,2)+TAUC(I)/2._q)
         ENDDO
      ELSEIF (NCDIJ==4) THEN
GPACC !$ACC PARALLEL LOOP REDUCTION(+:DE) PRESENT(DMUTOT,DKINEDEN,TAUC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO REDUCTION(+:DE) MAP(PRESENT,TO:DMUTOT,DKINEDEN,TAUC) IF(OFFLOAD_ON)
         DO I=1,GRIDC%RL%NP
            DE=DE+DMUTOT(I,1)*CHGCONJG(DKINEDEN(I,1)+TAUC(I)/2._q)+ &
           &       DMUTOT(I,2)*CHGCONJG(DKINEDEN(I,2))+DMUTOT(I,3)*CHGCONJG(DKINEDEN(I,3))+ &
           &        DMUTOT(I,4)*CHGCONJG((DKINEDEN(I,4)+TAUC(I)/2._q))
         ENDDO
      ELSE
         CALL vtutor%bug("MUxTAU: internal error, NCDIJ= " // str(NCDIJ), __FILE__, __LINE__)
      ENDIF

      DE=DE*RINPL*LATT_CUR%OMEGA/(RYTOEV*AUTOA2)

      CALLMPI( M_sum_d(GRIDC%COMM,DE,1))

      PROFILING_STOP('mu_x_tau')
      RETURN
      END SUBROUTINE MUxTAU

!=======================================================================
!
! SUBROUTINE TAU_PW
!
! This subroutine calculates the kinetic energy of the PW part of the
! wavefunctions (0.5*|grad psi|**2) on GRIDC
!
! INPUT: GRIDC,  LATT_CUR, SYMM, NIOND, W, WDES
! OUTPUT: TAU(GRIDC%MPLWC, WDES%NCDIJ)
!
! Robin Hirschl 20001221
!=======================================================================

    SUBROUTINE TAU_PW(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM,NIOND,W,WDES,TAU)

__GPU USE moffload

      USE prec
      USE lattice
      USE mgrid
      USE msymmetry
      USE base
      USE wave
      USE mpimy
      USE constant

      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC,GRID,GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (symmetry)    SYMM
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
! result
      RGRID TAU(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)  ! kinetic energy density
! work arrays
      RGRID,  ALLOCATABLE :: CKIN(:)
      COMPLEX(q), POINTER :: CFA(:),CFB(:),CFAR(:),CFBR(:)

      INTEGER :: ISPINOR,ISPINOR_
      INTEGER :: MPLWV,ISP,N,NK,NPL,IDIR,I,NIOND
      REAL(q) :: WEIGHT,G1,G2,G3
      COMPLEX(q) :: GC

#ifdef OFFLOADING
GPACC !$ACC WAIT IF(OFFLOAD_ON)
      CALL OFFLOAD_ON_PUSH(.TRUE.)

GPACC ACC_ASYNC_Q = ACC_ASYNC_ASYNC

      CALL COPYIN_TYPED_VAR(WDES)
      CALL COPYIN_TYPED_VAR(W)
      CALL COPYIN_TYPED_VAR(GRID)
      CALL COPYIN_TYPED_VAR(GRID_SOFT)
      CALL COPYIN_TYPED_VAR(GRIDC)
      CALL COPYIN_TYPED_VAR(SOFT_TO_C)

GPACC !$ACC ENTER DATA COPYIN(LATT_CUR) __IF_ASYNC__
GPACC !$ACC ENTER DATA CREATE(TAU) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:TAU) IF(OFFLOAD_ON)
#endif

      MPLWV=MAX(GRID%MPLWV,GRID_SOFT%MPLWV)

      ALLOCATE(CFA(WDES%NRPLWV),CFAR(MPLWV),CKIN(DIMREAL(MPLWV)))

      IF (WDES%LNONCOLLINEAR) THEN
! the non-collinear case needs work arrays for two separate spinor components
         ALLOCATE(CFB(WDES%NRPLWV),CFBR(MPLWV))
      ELSE
! in the collinear case there is no second spinor component
         CFB=>CFA ; CFBR=>CFAR
      ENDIF

GPACC !$ACC ENTER DATA CREATE(CFA,CFB,CFAR,CFBR,CKIN) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CFA,CFB,CFAR,CFBR,CKIN) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(TAU) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:TAU) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      TAU=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC UPDATE SELF(W%FERWE) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(W%FERWE) IF(OFFLOAD_ON)

spin: DO ISP=1,WDES%ISPIN

      DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

GPACC !$ACC KERNELS PRESENT(CFAR,CFBR,CKIN) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CFAR,CFBR,CKIN) IF(OFFLOAD_ON)
      CKIN=0 ; CFAR=0 ; IF (WDES%LNONCOLLINEAR) CFBR=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE kpoints
#endif
         NPL=WDES%NGVECTOR(NK)

         bands: DO N=1,WDES%NBANDS

            WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)/LATT_CUR%OMEGA/GRID_SOFT%NPLWV

            dirs: DO IDIR=1,3
GPACC !$ACC KERNELS PRESENT(CFA,CFB) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CFA,CFB) IF(OFFLOAD_ON)
               CFA=0 ; IF (WDES%LNONCOLLINEAR) CFB=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

! loop over plane wave coefficients
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(WDES,LATT_CUR,W,CFA) PRIVATE(G1,G2,G3,GC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES,LATT_CUR,W,CFA) &
GPOMP !$OMP PRIVATE(G1,G2,G3,GC) IF(OFFLOAD_ON)
               DO I=1,NPL
! get k-vector of respective k-point and coefficient
                  G1=WDES%IGX(I,NK)+WDES%VKPT(1,NK)
                  G2=WDES%IGY(I,NK)+WDES%VKPT(2,NK)
                  G3=WDES%IGZ(I,NK)+WDES%VKPT(3,NK)
                  GC=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))
                  CFA(I)=GC*W%CPTWFP(I+ISPINOR *NPL,N,NK,ISP)*CITPI
               ENDDO
! fourier trafo of gradient of wave-function
               CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CFAR(1),CFA(1),GRID)

! and for LNONCOLLINEAR=T repeat this procedure for the second spinor component
               IF (WDES%LNONCOLLINEAR) THEN
! loop over plane wave coefficients
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(WDES,LATT_CUR,W,CFB) PRIVATE(G1,G2,G3,GC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES,LATT_CUR,W,CFB) &
GPOMP !$OMP PRIVATE(G1,G2,G3,GC) IF(OFFLOAD_ON)
               DO I=1,NPL
! get k-vector of respective k-point and coefficient
                  G1=WDES%IGX(I,NK)+WDES%VKPT(1,NK)
                  G2=WDES%IGY(I,NK)+WDES%VKPT(2,NK)
                  G3=WDES%IGZ(I,NK)+WDES%VKPT(3,NK)
                  GC=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))
                  CFB(I)=GC*W%CPTWFP(I+ISPINOR_*NPL,N,NK,ISP)*CITPI
               ENDDO
! fourier trafo of gradient of wave-function
               CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CFBR(1),CFB(1),GRID)
               ENDIF
! update kinetic energy density (in CW2) and grad rho^2 (in CW4)
GPACC !$ACC PARALLEL LOOP PRESENT(CKIN,CFBR,CFAR,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CKIN,CFBR,CFAR,GRID) IF(OFFLOAD_ON)
               DO I=1,GRID%RL%NP
#ifdef realmode
                  CKIN(I)=CKIN(I)+HSQDTM*REAL(CONJG(CFBR(I))*CFAR(I),KIND=q)*WEIGHT
#else
                  CKIN(I)=CKIN(I)+HSQDTM*(CONJG(CFBR(I))*CFAR(I))*WEIGHT
#endif
               ENDDO

            ENDDO dirs

         ENDDO bands
      ENDDO kpoints

! merge results from nodes
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,  CKIN(1), GRID%RL%NP))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, CKIN(1), GRID%RL%NP))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER,  CKIN(1), GRID%RL%NP))
      CALLMPI( M_sum_z(WDES%COMM_KINTER, CKIN(1), GRID%RL%NP))
#endif

! ACC: this kernel caused the code to crash when executed asynchronously.
!      Fix: the operation is absorbed into the WEIGHT.
! rescaling
!!!$ACC KERNELS PRESENT(CKIN,GRID_SOFT) __IF_ASYNC__
!!      CKIN=CKIN/GRID_SOFT%NPLWV
!!!$ACC END KERNELS

! to rec space
      CALL FFT3D(CKIN(1),GRID_SOFT,-1)

! transition to finer grid
      CALL CPB_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CKIN(1),TAU(1,ISP+ISPINOR_+2*ISPINOR))
!     CALL SETUNB_COMPAT(TAU(1,ISP+ISPINOR_+2*ISPINOR),GRIDC)

      ENDDO ! ispinor_
      ENDDO ! ispinor
      ENDDO spin

! deallocate work arrays
GPACC !$ACC EXIT DATA DELETE(CFA,CFB,CFAR,CFBR,CKIN) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CFA,CFB,CFAR,CFBR,CKIN) IF(OFFLOAD_ON)
      DEALLOCATE(CFA,CFAR,CKIN) ; IF (WDES%LNONCOLLINEAR) DEALLOCATE(CFB,CFBR)
      NULLIFY(CFA,CFB,CFAR,CFBR)

      IF (WDES%LNONCOLLINEAR) THEN
         ! needs (total,mag) instead of up,dw
         CALL RC_FLIP(TAU,GRIDC,WDES%NCDIJ,.FALSE.)
         IF (SYMM%ISYM>0) THEN
            CALL RHOSYM(TAU(1,1),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,1)
            IF (.NOT.WDES%LSPIRAL) THEN
               CALL SYMFIELD(TAU(1,2),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
            ENDIF
         ENDIF
      ELSE
         IF (SYMM%ISYM>0) THEN
            ! symmetrization of result TAU(:,ISP)
            ! needs (total,mag) instead of up,dw
            CALL RC_FLIP(TAU, GRIDC,WDES%NCDIJ,.FALSE.)
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(TAU(1,ISP),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,ISP)
            ENDDO
            ! in the collinear case we go back to up,dw
            CALL RC_FLIP(TAU,GRIDC,WDES%NCDIJ,.TRUE.)
         ENDIF
      ENDIF

! Apply a spherical cutoff in reciprocal space
      DO ISP=1,WDES%NCDIJ
         CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, TAU(1,ISP))
      ENDDO

#ifdef OFFLOADING
      CALL DELETE_TYPED_VAR(W)
      CALL DELETE_TYPED_VAR(WDES)
      CALL DELETE_TYPED_VAR(GRID)
      CALL DELETE_TYPED_VAR(GRID_SOFT)
      CALL DELETE_TYPED_VAR(GRIDC)
      CALL DELETE_TYPED_VAR(SOFT_TO_C)

GPACC !$ACC EXIT DATA DELETE(LATT_CUR) __IF_ASYNC__
GPACC !$ACC EXIT DATA COPYOUT(TAU) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LATT_CUR) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:TAU) IF(OFFLOAD_ON)

GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC ACC_ASYNC_Q = ACC_ASYNC_SYNC

      CALL OFFLOAD_ON_POP
#endif

      RETURN
    END SUBROUTINE TAU_PW
