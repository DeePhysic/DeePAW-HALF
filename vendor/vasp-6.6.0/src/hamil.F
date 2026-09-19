#include "symbol.inc"
!***********************************************************************
!
!> this module implements most of the low and high level
!> routines to calculated the action of a local Hamiltonian
!> onto the wavefunctions
!
!***********************************************************************
MODULE hamil

__GPU USE moffload

  USE prec
  USE wave_high

  INTERFACE
     SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)
       USE prec
       USE wave
       REAL(q) FNORM
       TYPE (wavedes1)    WDES1
       COMPLEX(q) :: C
       REAL(q), OPTIONAL :: FMETRIC
       REAL(q), OPTIONAL :: METRIC(WDES1%NPL)
     END SUBROUTINE PW_NORM_WITH_METRIC
  END INTERFACE

  INTERFACE
     SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CW, CH)
       USE mgrid
       USE wave

       TYPE (wavedes1)    WDES1
       TYPE (grid_3d)     GRID
       COMPLEX(q) :: CVR
       LOGICAL    :: LADD
       REAL(q)    :: DATAKE
       REAL(q)    :: EVALUE
       COMPLEX(q) :: CW
       COMPLEX(q) :: CH
     END SUBROUTINE KINHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       RGRID   SV
       COMPLEX(q) :: CR
       COMPLEX(q) :: CVR
     END SUBROUTINE VHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       RGRID   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       GDEF   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE_CMPLX
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)
       USE prec
       USE mgrid
       USE wave

       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1

       GDEF   SV
       COMPLEX(q) :: CR,CVR
       REAL(q)    :: WEIGHT
     END SUBROUTINE VHAMIL_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       GDEF   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       GDEF   CHARGE
       GDEF :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_GDEF
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       GDEF   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
       USE prec
       IMPLICIT NONE
       COMPLEX(q)  CNL
       INTEGER LMMAXC, LMDIM
       OVERLAP CDIJ,CQIJ
       REAL(q) EVALUE
       GDEF CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL_
  END INTERFACE

CONTAINS

!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions st be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!
!***********************************************************************

  SUBROUTINE ECCP(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR, ISPINOR_
    INTEGER NPRO,NPRO_

    GDEF    CNL
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) !< local potential
    OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
GPACC !$ACC ROUTINE (ECCP_NL) VECTOR

    PROFILING_START('eccp')

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CE) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CE) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE) IF(OFFLOAD_ON)
    CE=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,W1,W2,WDES1,SV,CE) &
GPACC !$ACC PRIVATE(MM,MM_) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:GRID,W1,W2,WDES1,SV,CE) &
GPOMP !$OMP PRIVATE(MM,MM_) REDUCTION(+:CE) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2,NGVECTOR,LMDIM,CDIJ) REDUCTION(+:CE) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_,NPRO,NPRO_,NI,NT,CNL)

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CE=CE+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO

GPACC !$ACC KERNELS PRESENT(CE,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE,GRID) IF(OFFLOAD_ON)
    CE=CE/GRID%NPLWV
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! kinetic energy contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1,W2,WDES1,CE) &
GPACC !$ACC PRIVATE(MM) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W1,W2,WDES1,CE) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CE) IF(OFFLOAD_ON)

    DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CE=CE+W1%CW(MM)*CONJG(W2%CW(MM))*WDES1%DATAKE(M,ISPINOR+1)
       ENDDO
__CPU !$OMP END DO
    ENDDO
!=======================================================================
! non local contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CE) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NT,CNL) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:W1,W2,CDIJ,WDES1) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NT,CNL) REDUCTION(+:CE) IF(OFFLOAD_ON)
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             IF (WDES1%LMMAX(NT)==0) CYCLE
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             CNL=0

! currently it is better to inline manually the call to ECCP_NL
GPACC !$ACC LOOP VECTOR COLLAPSE(2) REDUCTION(+:CNL)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:CNL)
               DO L=1,WDES1%LMMAX(NT)
                  DO LP=1,WDES1%LMMAX(NT)
                     CNL=CNL+CDIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*W1%CPROJ(LP+NPRO_)*GCONJG(W2%CPROJ(L+NPRO))
                  ENDDO
               ENDDO
             CE=CE+CNL
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO spinor
__CPU !$OMP END PARALLEL

GPACC !$ACC KERNELS PRESENT(CE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE) IF(OFFLOAD_ON)
    CE=GREAL(CE)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

GPACC !$ACC UPDATE SELF(CE) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC EXIT DATA DELETE(CE) __IF_ASYNC__

GPOMP !$OMP TARGET UPDATE FROM(CE) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CE) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp')

  END SUBROUTINE ECCP

!************************* SUBROUTINE ECCP_VEC *************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions must be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!> @n
!> Furthermore the workarrays CWORK1 and CWORK2 acquire a second
!> dimension (size=3).  This additional workspace improves the
!> performance of ::kinhamil_vec under OpenMP parallelization.
!
!***********************************************************************

  SUBROUTINE ECCP_VEC(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,AVEC, CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)


    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR
    INTEGER NPRO,NPRO_
!$  INTEGER NPRO2,NPRO2_

    GDEF    CNL
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS)
    RGRID   AVEC(:,:)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CVR,CWORK1,CWORK2
#endif

    OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
GPACC !$ACC ROUTINE (ECCP_NL) VECTOR

    PROFILING_START('eccp_vec')

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CLOCAL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CLOCAL) IF(OFFLOAD_ON)
    CLOCAL=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,WDES1,SV,W1,W2,CLOCAL) &
GPACC !$ACC PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:GRID,WDES1,SV,W1,W2,CLOCAL) &
GPOMP !$OMP PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2) PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_) REDUCTION(+:CLOCAL)
    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO
__CPU !$OMP END PARALLEL

GPACC !$ACC KERNELS PRESENT(CLOCAL,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CLOCAL) IF(OFFLOAD_ON)
    CLOCAL=CLOCAL/GRID%NPLWV
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA COPYOUT(CLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CLOCAL) IF(OFFLOAD_ON)

!=======================================================================
! kinetic energy contribution
!=======================================================================
GPACC !$ACC ENTER DATA CREATE(CKIN,CVR,CWORK1,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CKIN,CVR,CWORK1,CWORK2) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVR) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CVR=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALL KINHAMIL_VEC( WDES1, GRID, CVR,  .FALSE., &
         WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
         CWORK1(1 __NOOFF_omp_arg(1)), CWORK2(1 __NOOFF_omp_arg(1)), 0.0_q, W1%CW(1), CWORK1(1 __NOOFF_omp_arg(1)))

GPACC !$ACC EXIT DATA DELETE(CVR,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CVR,CWORK2) IF(OFFLOAD_ON)

GPACC !$ACC ENTER DATA CREATE(CKIN,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CKIN,CTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CKIN,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CKIN,CTMP) IF(OFFLOAD_ON)
    CKIN=0 ; CTMP=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CKIN,CWORK1,W2) &
GPACC !$ACC PRIVATE(MM) REDUCTION(+:CKIN) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CKIN,CWORK1,W2) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CKIN) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,NGVECTOR,CWORK1,W1,W2,LMDIM,CDIJ) REDUCTION(+:CKIN,CTMP) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,NRPO,NPRO_,NI,NT,LMMAXC,CNL)
    DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CWORK1(MM __NOOFF_omp_arg(1))*CONJG(W2%CW(MM))
       ENDDO
__CPU !$OMP END DO
    ENDDO

GPACC !$ACC EXIT DATA COPYOUT(CKIN) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(CWORK1) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CKIN) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

!=======================================================================
! non local contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CTMP) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) __IF_ASYNC__

!! INTEL COMPILER (IFX) 2024.2.1 20240711 crashes when IF(OFFLOAD_ON) is used in the following directive
!! unless the compiler switch "-switch no-use-host-usm-for-implicit-reduction-map" is used.
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:W1,W2,CDIJ,WDES1,CTMP) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) IF(OFFLOAD_ON)

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CNL=0; CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
             CTMP=CTMP+CNL
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO spinor
__CPU !$OMP END PARALLEL

GPACC !$ACC EXIT DATA COPYOUT(CTMP) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CTMP) IF(OFFLOAD_ON)

    CE=GREAL(CLOCAL+CKIN+CTMP)
    CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

    PROFILING_STOP('eccp_vec')

  END SUBROUTINE ECCP_VEC


!************************* SUBROUTINE ECCP_TAU *************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions must be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!> @n
!> Furthermore the workarrays CWORK1 and CWORK2 acquire a second
!> dimension (size=3 and 6, respectively). This additional workspace
!> improves the performance of ::kinhamil_tau under OpenMP
!> parallelization.
!
!***********************************************************************

  SUBROUTINE ECCP_TAU(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,LATT_CUR,MU,CE)
    USE prec
    USE mpimy
    USE mgrid
    USE lattice
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID
    TYPE (latt)     :: LATT_CUR

    INTEGER NGVECTOR, ISPINOR
    INTEGER NPRO,NPRO_
!$  INTEGER NPRO2,NPRO2_

    GDEF    CNL
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS)
    RGRID   MU(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CH(WDES1%NRPLWV)
    COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CVR,CH,CWORK1,CWORK2
#endif

    OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
GPACC !$ACC ROUTINE (ECCP_NL) VECTOR

    PROFILING_START('eccp_tau')

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CE,CLOCAL,CKIN,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CE,CLOCAL,CKIN,CTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CLOCAL) IF(OFFLOAD_ON)
    CLOCAL=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,WDES1,SV,W1,W2,CLOCAL) &
GPACC !$ACC PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:GRID,WDES1,SV,W1,W2,CLOCAL) &
GPOMP !$OMP PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2) PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_) REDUCTION(+:CLOCAL)
    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO
__CPU !$OMP END PARALLEL

GPACC !$ACC KERNELS PRESENT(CLOCAL,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CLOCAL,GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CLOCAL=CLOCAL/GRID%NPLWV
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

!=======================================================================
! kinetic energy contribution
!=======================================================================
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__
GPACC !$ACC ENTER DATA CREATE(CVR,CH,CWORK1,CWORK2) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CVR,CH,CWORK1,CWORK2) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CVR,CH,CKIN) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVR,CH,CKIN) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CKIN=0; CVR=0; CH=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALL KINHAMIL_TAU( WDES1, GRID, CVR,  .FALSE., .TRUE., &
         WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,1), &
         CWORK1(1 __NOOFF_omp_arg(1)), CWORK2(1 __NOOFF_omp_arg(1)), 0.0_q, W1%CW(1), CH(1))

GPACC !$ACC EXIT DATA DELETE(LATT_CUR%B,LATT_CUR,CVR,CWORK1,CWORK2) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LATT_CUR,LATT_CUR%B,CVR,CWORK1,CWORK2) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CTMP=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES1,W2,CH,CKIN) &
GPACC !$ACC PRIVATE(MM) REDUCTION(+:CKIN) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:WDES1,W2,CH,CKIN) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CKIN) IF(OFFLOAD_ON)

#ifndef OFFLOADING
!$OMP PARALLEL &
!$OMP SHARED(WDES1,NGVECTOR,CH,W1,W2,LMDIM,CDIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,NRPO,NPRO_,NI,NT,LMMAXC,CNL) &
!$OMP REDUCTION(+:CKIN,CTMP)
#endif

    DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CH(MM)*CONJG(W2%CW(MM))
       ENDDO
__CPU !$OMP END DO
    ENDDO

GPACC !$ACC EXIT DATA DELETE(CH) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CH) IF(OFFLOAD_ON)
!=======================================================================
! non local contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CTMP) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:W1,W2,CDIJ,WDES1,CTMP) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) IF(OFFLOAD_ON)

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CNL=0; CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
             CTMP=CTMP+CNL
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO spinor
__CPU !$OMP END PARALLEL


GPACC !$ACC KERNELS PRESENT(CE,CLOCAL,CKIN,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE,CLOCAL,CKIN,CTMP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CE=GREAL(CLOCAL+CKIN+CTMP)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA DELETE(CLOCAL,CKIN,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CLOCAL,CKIN,CTMP) IF(OFFLOAD_ON)

    CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

GPACC !$ACC UPDATE SELF(CE) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC!$ACC EXIT DATA DELETE(CE) __IF_ASYNC__

GPOMP !$OMP TARGET UPDATE FROM(CE) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CE) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp_tau')

  END SUBROUTINE ECCP_TAU


!*********************************************************************
!
!> calculate the band structure energy
!
!*********************************************************************

  FUNCTION BANDSTRUCTURE_ENERGY(WDES, W) RESULT(E)
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W      ! wavefunction
    REAL(q)            E
    INTEGER            ISP, NK, NB

    E=0
    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          DO NB=1,WDES%NB_TOT
             E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) *WDES%WTKPT(NK)*W%FERTOT(NB,NK,ISP)
          ENDDO
       ENDDO
    ENDDO
    CALLMPI( M_sum_d(WDES%COMM_KINTER, E, 1))
  END FUNCTION BANDSTRUCTURE_ENERGY


!**********************************************************************
!
!> calculate the kinetic energy of each wavefunction
!
!**********************************************************************

  SUBROUTINE KINETIC_ENERGY(W )
    USE wave_high
    IMPLICIT NONE
    TYPE (wavespin), TARGET :: W
    !    local
    TYPE (wavedes1)    WDES1
    TYPE (wavedes), POINTER :: WDES

    INTEGER ISP, NK, NB, ISPINOR, MM, M
    REAL(q) :: CKIN


    WDES=>W%WDES

    spin:  DO ISP=1,WDES%ISPIN
       kpoint: DO NK=1,WDES%NKPTS

#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!PK MRG_AUX already zeros contributions from non-kpoint nodes, but be cautious and zero workspace anyway
             DO NB=1,WDES%NBANDS
                W%AUX(NB,NK,ISP)=0.
             END DO
          ELSE
#endif
          CALL SETWDES(WDES,WDES1,NK)

          DO NB=1,WDES%NBANDS
             CKIN=0
             DO ISPINOR=0,WDES1%NRSPINORS-1
                DO M=1,WDES1%NGVECTOR
                   MM=M+ISPINOR*WDES1%NGVECTOR
                   CKIN=CKIN+W%CW(MM,NB,NK,ISP)*CONJG(W%CW(MM,NB,NK,ISP))*WDES1%DATAKE(M,ISPINOR+1)
                ENDDO
             ENDDO
             W%AUX(NB,NK,ISP)=CKIN
          ENDDO
#ifdef MPI
          ENDIF
#endif
       END DO kpoint
    END DO spin

    CALL MRG_AUX(WDES,W)
  END SUBROUTINE KINETIC_ENERGY

!************************* SUBROUTINE SETUP_PRECOND ********************
!
!> subroutine to set up a preconditioning matrix
!> IALGO determines the type of preconditioning
!>
!>-   0,6,8 TAP preconditioning
!>-    9    Jacobi like preconditioning
!>-   else  no preconditioning
!>
!>- the Jacobi like preconditer required the eigenvalue and
!>- a complex shift
!
!***********************************************************************


  SUBROUTINE SETUP_PRECOND( W1, IALGO, IDUMP, PRECON, EVALUE, DE_ATT )
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    INTEGER IALGO         !< chosen algorithm
    INTEGER IDUMP         !< dump flag
    REAL(q) :: EVALUE     !< eigenvalue minus average local potential
    REAL(q) :: DE_ATT     !< complex shift
    REAL(q) :: PRECON(*)  !< return: preconditioning matrix
    ! local
    INTEGER NRSPINORS, NGVECTOR
    INTEGER M, MM, NK, ISPINOR
    REAL(q) :: EKIN, FAKT, X, X2
    COMPLEX(q) :: CPT

    PROFILING_START('setup_precond')

    NRSPINORS= W1%WDES1%NRSPINORS
    NGVECTOR = W1%WDES1%NGVECTOR

GPACC !$ACC ENTER DATA COPYIN(W1,W1%CW,W1%WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:W1,W1%CW,W1%WDES1) IF(OFFLOAD_ON)

    IF (IALGO==0 .OR. IALGO==8 .OR. IALGO==6) THEN
       EKIN=0

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR REDUCTION(+:EKIN) PRIVATE(MM,CPT) &
GPACC !$ACC PRESENT(W1) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W1) &
GPOMP !$OMP PRIVATE(MM,CPT) REDUCTION(+:EKIN) IF(OFFLOAD_ON)

       DO ISPINOR=0,NRSPINORS-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP PARALLEL DO PRIVATE(M,MM,CPT) SHARED(W1,ISPINOR) REDUCTION(+:EKIN)
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
             CPT=W1%CW(MM)
             EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * W1%WDES1%DATAKE(M,ISPINOR+1)
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CALLMPI( M_sum_d(W1%WDES1%COMM_INB, EKIN, 1))

       IF (EKIN<2.0_q) EKIN=2.0_q
       EKIN=EKIN*1.5_q
       IF (IDUMP==2)  WRITE(*,'(E9.2,"K")',ADVANCE='NO') EKIN

       FAKT=2._q/EKIN

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM,X,X2) &
GPACC !$ACC PRESENT(PRECON,W1) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,ALLOC:PRECON(1:NGVECTOR*NRSPINORS),W1) &
GPOMP !$OMP PRIVATE(MM,X,X2) IF(OFFLOAD_ON)

       DO ISPINOR=0,NRSPINORS-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP PARALLEL DO PRIVATE(M,MM,X,X2) SHARED(W1,ISPINOR,EKIN,FAKT,PRECON)
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
             X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
             X2= 27+X*(18+X*(12+8*X))
             PRECON(MM)=X2/(X2+16*X*X*X*X)*FAKT
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

    ELSE IF (IALGO==9) THEN

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM,X) &
GPACC !$ACC PRESENT(PRECON,W1) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:PRECON(1:NGVECTOR*NRSPINORS),W1) &
GPOMP !$OMP PRIVATE(MM,X) IF(OFFLOAD_ON)

       DO ISPINOR=0,NRSPINORS-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP PARALLEL DO PRIVATE(M,MM,X) SHARED(W1,EVALUE,DE_ATT,PRECON)
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
             X=MAX(W1%WDES1%DATAKE(M,ISPINOR+1)-EVALUE,0._q)
             PRECON(MM)= REAL( 1._q/(X+ CMPLX( 0._q , DE_ATT ,KIND=q) ) ,KIND=q) !new
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

    ELSE

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM) &
GPACC !$ACC PRESENT(PRECON,W1) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:PRECON(1:NGVECTOR*NRSPINORS),W1) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

       DO ISPINOR=0,NRSPINORS-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(W1,PRECON)
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
             PRECON(MM)=1.0_q
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

    ENDIF

GPACC !$ACC EXIT DATA DELETE(W1,W1%CW,W1%WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1,W1%CW,W1%WDES1) IF(OFFLOAD_ON)
    PROFILING_STOP('setup_precond')

  END SUBROUTINE SETUP_PRECOND


!************************* SUBROUTINE APPLY_PRECOND ********************
!
!> apply W2%CW <- W1%CW * PRECON * MUL
!
!***********************************************************************

  SUBROUTINE APPLY_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    PROFILING_START('apply_precond')

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1._q
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CW(1) , MUL_, W2%CW(1), 0.0_q, W2%CW(1))

    PROFILING_STOP('apply_precond')

  END SUBROUTINE APPLY_PRECOND

!************************* SUBROUTINE ADD_PRECOND **********************
!
!> apply W2%CW < W2%CW + W1%CW * PRECON * MUL
!
!***********************************************************************

  SUBROUTINE ADD_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    PROFILING_START('add_precond')

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1.0
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CW(1) , MUL_, W2%CW(1), 1.0_q, W2%CW(1))

    PROFILING_STOP('add_precond')

  END SUBROUTINE ADD_PRECOND


  SUBROUTINE SIMPLE_PRECOND( W1 )
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1

    REAL(q) :: EKIN, FAKT, X, X2
    INTEGER :: ISPINOR, M, MM

! old version
    EKIN=10
    FAKT=2._q/EKIN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1) PRIVATE(MM,X,X2) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W1) &
GPOMP !$OMP PRIVATE(MM,X,X2) IF(OFFLOAD_ON)
    DO ISPINOR=0,W1%WDES1%NRSPINORS-1
__CPU !$OMP PARALLEL DO &
__CPU !$OMP SHARED(W1,ISPINOR,EKIN,FAKT) &
__CPU !$OMP PRIVATE(M,MM,X,X2)
       DO M=1,W1%WDES1%NGVECTOR
          MM=M+ISPINOR*W1%WDES1%NGVECTOR
          X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
          X2=27+X*(18+X*(12+8*X))
          W1%CW(MM)= W1%CW(MM)*X2/(X2+16*X*X*X*X)*FAKT
       ENDDO
__CPU !$OMP END PARALLEL DO
    ENDDO

  END SUBROUTINE SIMPLE_PRECOND


!***********************************************************************
!
!> truncate high frequency components when LDELAY is set
!> and/or in case of spin-spiral calculations (LSPIRAL=.TRUE.)
!
!***********************************************************************

  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1( W1, LDELAY, ENINI)
    USE prec
    IMPLICIT NONE
    LOGICAL LDELAY     !< usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      !< cutoff at which the wavefunction is truncated
    TYPE (wavefun1)    W1
    ! local
    INTEGER ISPINOR

    PROFILING_START('truncate_high_frequency_w1')

    IF (LDELAY.OR.W1%WDES1%LSPIRAL) THEN
       DO ISPINOR=1,W1%WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(W1%WDES1%NGVECTOR, &
               W1%CW(1+(ISPINOR-1)*W1%WDES1%NGVECTOR), W1%WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF

    PROFILING_STOP('truncate_high_frequency_w1')

  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1


!***********************************************************************
!
!> determine norm of a wavefunction and norm with a metric
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC_W1(W1, FNORM, FMETRIC, METRIC)
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavefun1)    W1
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(W1%WDES1%NPL)

    PROFILING_START('pw_norm_with_metric')

    IF (PRESENT(FMETRIC) .AND. PRESENT(METRIC)) THEN
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CW(1), FNORM, FMETRIC, METRIC)
    ELSE
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CW(1), FNORM)
    ENDIF

    PROFILING_STOP('pw_norm_with_metric')
  END SUBROUTINE PW_NORM_WITH_METRIC_W1



!************************* SUBROUTINE HAMILT ***************************
!
!> this subroutine calculates the H acting onto a wavefuntion
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE Q |phi>
!> where Q is the nondiagonal part of the overlap matrix S
!***********************************************************************

  SUBROUTINE HAMILT( &
       &    W1, NONLR_S, NONL_S, EVALUE, CDIJ, CQIJ, SV, ISP, CH)
    USE mpimy
    USE mgrid

    USE nonl_high
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S

    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    INTEGER    ISP
    COMPLEX(q) CH(:)
    REAL(q)    EVALUE
    ! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1
#endif

    IF (NONLR_S%LREAL) THEN
       ! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1))
       ! non-local contribution in real-space

       CALL RACC(NONLR_S, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .FALSE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))
    ELSE
       ! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1))

       ! calculate the non local contribution in reciprocal space
       CALL VNLACC(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE, CH)
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .TRUE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))

    ENDIF
    RETURN
  END SUBROUTINE HAMILT


!************************* SUBROUTINE HAMILTMU *************************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE S |phi>
!> where S is the overlap matrix S
!
!***********************************************************************

  SUBROUTINE HAMILTMU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    REAL(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP,NSTRIP,NSTRIP_ACT
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1
#endif

    PROFILING_START('hamiltmu')

GPACC !$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ELSE
       ! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu')

    RETURN
  END SUBROUTINE HAMILTMU


!************************* SUBROUTINE HAMILTMU_VEC *********************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE S |phi>
!> where S is the overlap matrix S
!> this version includes a vector potential A
!
!> @details @ref openmp :
!> under OpenMP the workarrays CWORK2 and CWORK3 acquire a second
!> dimension (size=3). This additional workspace improves the
!> performance of ::kinhamil_vec under OpenMP parallelization.
!
!***********************************************************************

  SUBROUTINE HAMILTMU_VEC( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, AVEC, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    REAL(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    RGRID      AVEC(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
#ifdef OFFLOADING
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __omp_arg(3))
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1,CWORK2,CWORK3
#endif
    INTEGER NP

    PROFILING_START('hamiltmu_vec')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
#ifdef OFFLOADING
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             ! calculate the non local contribution in reciprocal space
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
#ifdef OFFLOADING
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_vec')

    RETURN
  END SUBROUTINE HAMILTMU_VEC


!************************* SUBROUTINE HAMILTMU_TAU *********************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!> for a kinetic energy functional
!>
!> @details @ref openmp :
!> under OpenMP the workarrays CWORK2 and CWORK3 acquire a second
!> dimension (size=3 and 6, respectively). This additional workspace
!> improves the performance of ::kinhamil_tau under OpenMP
!> parallelization.
!
!***********************************************************************

  SUBROUTINE HAMILTMU_TAU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, LATT_CUR, MU, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE lattice
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (latt)        LATT_CUR
    REAL(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    RGRID      MU(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
#ifdef OFFLOADING
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __omp_arg(6))
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1,CWORK2,CWORK3
#endif

    INTEGER NP

    PROFILING_START('hamiltmu_tau')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO


    IF (NONLR_S%LREAL) THEN
      ! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef OFFLOADING
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             ! calculate the non local contribution in reciprocal space
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef OFFLOADING
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_tau')

    RETURN
  END SUBROUTINE HAMILTMU_TAU


!************************* SUBROUTINE HAMILTMU_C ***********************
!
!> identical to previous subroutine but with a complex eigenvalue EVALUE
!
!***********************************************************************

  SUBROUTINE HAMILTMU_C( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    COMPLEX(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1
#endif

    PROFILING_START('hamiltmu_c')

!$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ELSE
       ! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1) IF(OFFLOAD_ON)
    ACC_ASYNC_Q = ACC_ASYNC_SYNC
#endif

    PROFILING_STOP('hamiltmu_c')

    RETURN
  END SUBROUTINE HAMILTMU_C


!************************* SUBROUTINE HAMILTMU_VEC_C *******************
!
!> identical to HAMILTMU_VEC but with a complex eigenvalue EVALUE
!
!***********************************************************************

  SUBROUTINE HAMILTMU_VEC_C( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, AVEC, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    COMPLEX(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    RGRID      AVEC(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
#ifdef ACC_OFFLOAD
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __omp_arg(3))
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1,CWORK2,CWORK3
#endif
    INTEGER NP

    PROFILING_START('hamiltmu_vec_c')

!$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_VEC_C( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
#ifdef ACC_OFFLOAD
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             ! calculate the non local contribution in reciprocal space
             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL_VEC_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
#ifdef ACC_OFFLOAD
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif

    PROFILING_STOP('hamiltmu_vec_c')

    RETURN
  END SUBROUTINE HAMILTMU_VEC_C


!************************* SUBROUTINE HAMILTMU_TAU_C *******************
!
!> identical to HAMILTMU_TAU but with a complex eigenvalue EVALUE
!
!***********************************************************************

  SUBROUTINE HAMILTMU_TAU_C( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, LATT_CUR, MU, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE lattice
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (latt)        LATT_CUR
    COMPLEX(q)    EVALUE(:)               ! eigenvalues
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    RGRID      SV(:,:)
    RGRID      MU(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
#ifdef ACC_OFFLOAD
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __omp_arg(6))
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1,CWORK2,CWORK3
#endif

    INTEGER NP

    PROFILING_START('hamiltmu_tau_c')

!$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
!$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO


    IF (NONLR_S%LREAL) THEN
      ! non-local contribution in real-space
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_TAU_C( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef ACC_OFFLOAD
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             ! calculate the non local contribution in reciprocal space
             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CW(:,NP))
             CALL KINHAMIL_TAU_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef ACC_OFFLOAD
                  CWORK2(1,NP), CWORK3(1,NP), &
#else
                  CWORK2(1 __omp_arg(1)),  CWORK3(1 __omp_arg(1)), &
#endif
                  EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif

    PROFILING_STOP('hamiltmu_tau_c')

    RETURN
  END SUBROUTINE HAMILTMU_TAU_C


!************************* SUBROUTINE HAMILT_LOCAL *********************
!
!> this subroutine calculates the local and kinetic energy part of H
!> acting onto a wavefuntion
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> LADD=.TRUE.
!> ~~~
!>  CH = CH + SV* W1 + T * W1
!> ~~~
!> LADD=.FALSE.
!>  CH = SV* W1 + T * W1
!> LKIN determines whether the kinetic energy contribution is added or not
!> the default is add the kinetic energy contribution
!
!***********************************************************************

  SUBROUTINE HAMILT_LOCAL(W1, SV, ISP, CH, LADD, LKIN)
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    RGRID   SV(:,:)
    INTEGER :: ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
    ! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1
#endif

    PROFILING_START('hamilt_local')

GPACC !$ACC DATA CREATE(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1))
    IF (.NOT. PRESENT(LKIN)) THEN
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))
    ELSE
       IF (LKIN) THEN
          CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
               W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))
       ELSE
          CALL FFTHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
                                     0.0_q, W1%CW(1), CH(1))
       ENDIF
    ENDIF

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC END DATA
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_local')

  END SUBROUTINE HAMILT_LOCAL


  SUBROUTINE HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, MU, ISP, CH, LADD, LKIN)
    USE mgrid
    USE lattice
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (latt)        LATT_CUR
    RGRID   SV(:,:)
    RGRID   MU(:,:)
    INTEGER ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
    ! local variables
    COMPLEX(q) :: CVR(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(__align_DIM(W1%WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK2(__align_DIM(W1%WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
    LOGICAL :: LKINETIC=.TRUE.
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CVR,CWORK1,CWORK2
#endif
    PROFILING_START('hamilt_local_tau')

    IF (PRESENT(LKIN)) LKINETIC=LKIN
GPACC !$ACC ENTER DATA CREATE(CVR,CWORK1,CWORK2) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CVR,CWORK1,CWORK2) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CVR(1))

    CALL KINHAMIL_TAU( W1%WDES1, W1%WDES1%GRID, CVR(1), LADD, LKINETIC, W1%WDES1%DATAKE(1,1), &
   &   W1%WDES1%IGX(1), W1%WDES1%IGY(1), W1%WDES1%IGZ(1), W1%WDES1%VKPT(1), LATT_CUR, MU(1,ISP), CWORK1, CWORK2, &
   &   0.0_q, W1%CW(1), CH(1))

GPACC !$ACC EXIT DATA DELETE(CVR,CWORK1,CWORK2,LATT_CUR%B,LATT_CUR) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CVR,CWORK1,CWORK2) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_local_tau')

  END SUBROUTINE HAMILT_LOCAL_TAU


!************************* SUBROUTINE ECCP_NL_ALL **********************
!
!> this subroutine calculates the non local contribution to the
!> expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
!> where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)

__GPU USE moffload_struct_def

    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    OVERLAP CDIJ(:,:,:,:)
    OVERLAP CQIJ(:,:,:,:)
    REAL(q) EVALUE
    ! return
    COMPLEX(q) CE
    ! local
    COMPLEX(q) CTMP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

!$ACC ROUTINE(ECCP_NL_) VECTOR

!$ACC KERNELS PRESENT(CE) __IF_ASYNC__
    CE=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP REDUCTION(+:CE) PRESENT(CE,CDIJ,CQIJ,W1%CPROJ,W2%CPROJ) &
!$ACC COLLAPSE(3) PRIVATE(CTMP,NPRO,NPRO_,NT,LMMAXC) __IF_ASYNC__
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CTMP=0
             CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                  CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                  W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CTMP)
             CE=CE+CTMP
          ENDDO
       ENDDO
    ENDDO spinor

    CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

    RETURN
  END SUBROUTINE ECCP_NL_ALL


!************************* SUBROUTINE ECCP_NL_ALL_ION ******************
!
!> this subroutine calculates the non local contribution to the
!> expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
!> where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL_ION(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    OVERLAP CDIJ(:,:,:,:)
    OVERLAP CQIJ(:,:,:,:)
    REAL(q) EVALUE
    ! return
    COMPLEX(q) CE(*)
    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                     CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                     W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CE(NI))
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    RETURN
  END SUBROUTINE ECCP_NL_ALL_ION

#ifndef oldhamil
!************************* SUBROUTINE ECCP_ ****************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions st be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!
!***********************************************************************

  SUBROUTINE ECCP_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1, W2, LMDIM, CDIJ, SV, ISP, CE)

    USE hamil_struct_def
    USE lattice
    USE mgrid

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (grid_3d) :: GRID
    TYPE (latt) :: LATT_CUR
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1, W2

    RGRID :: SV(:,:)
    OVERLAP :: CDIJ(:,:,:,:)

    INTEGER :: LMDIM
    INTEGER :: ISP

    COMPLEX(q) :: CE !< expectation value

    ! local variables
    COMPLEX(q) :: CH(W1%WDES1%NRPLWV)
    COMPLEX(q) :: CTMP
    GDEF :: CNL
    INTEGER :: NGVECTOR, ISPINOR, ISPINOR_, M, MM, MM_
    INTEGER :: NI, NT, L, LP, NPRO, NPRO_

GPACC !$ACC ROUTINE (ECCP_NL) VECTOR

    PROFILING_START('eccp_')

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CE) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CE) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE) IF(OFFLOAD_ON)
    CE=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,W1,W2,WDES1,SV,CE) &
GPACC !$ACC PRIVATE(MM,MM_) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:GRID,W1,W2,WDES1,SV,CE) &
GPOMP !$OMP PRIVATE(MM,MM_) REDUCTION(+:CE) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2) REDUCTION(+:CE) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CE=CE+SV(M,ISP+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO
__CPU !$OMP END PARALLEL

    IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) THEN

GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP) IF(OFFLOAD_ON)
       CTMP=0
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,W1,W2,WDES1,HAMILTONIAN%MUJVEC,CTMP) &
GPACC !$ACC PRIVATE(MM,MM_) REDUCTION(+:CTMP) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:GRID,W1,W2,WDES1,HAMILTONIAN%MUJVEC,CTMP) &
GPOMP !$OMP PRIVATE(MM,MM_) REDUCTION(+:CTMP) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,GRID,HAMILTONIAN,W1,W2) REDUCTION(+:CTMP) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
       DO ISPINOR =0,WDES1%NRSPINORS-1
          DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
             DO M=1,GRID%RL%NP
                MM =M+ISPINOR *GRID%MPLWV
                MM_=M+ISPINOR_*GRID%MPLWV
                CTMP=CTMP+HAMILTONIAN%MUJVEC(M,ISP+ISPINOR_+2*ISPINOR,4) *W1%CR(MM_)*CONJG(W2%CR(MM))
             ENDDO
__CPU !$OMP END DO
          ENDDO
       ENDDO
__CPU !$OMP END PARALLEL

GPACC !$ACC KERNELS PRESENT(CE,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE,CTMP) IF(OFFLOAD_ON)
       CE=CE-CTMP*(0._q,1._q)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
    ENDIF

GPACC !$ACC KERNELS PRESENT(CE,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE,GRID) IF(OFFLOAD_ON)
    CE=CE/GRID%NPLWV
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

!=======================================================================
! add contribution from vector potential (if required)
!=======================================================================
    IF (ASSOCIATED(HAMILTONIAN%AVEC)) &
       CALL ECCP_VEC_(HAMILTONIAN, GRID, WDES1, W1, W2, CE)

!!!$ACC KERNELS PRESENT(CE,GRID) __IF_ASYNC__
!!    CE=CE/GRID%NPLWV
!!!$ACC END KERNELS

!=======================================================================
! kinetic energy contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1,W2,WDES1,CE) &
GPACC !$ACC PRIVATE(MM) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W1,W2,WDES1,CE) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CE) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,W1,W2,NGVECTOR) PRIVATE(ISPINOR,M,MM) REDUCTION(+:CE)
    DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CE=CE+W1%CW(MM)*CONJG(W2%CW(MM))*WDES1%DATAKE(M,ISPINOR+1)
       ENDDO
__CPU !$OMP END DO
    ENDDO
__CPU !$OMP END PARALLEL

!=======================================================================
! add mGGA and spin-current contributions (if required)
!=======================================================================
    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
       IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) THEN
          ! both mGGA and spin-current
          CALL ECCP_JVEC_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1, W2, ISP, CE)
       ELSE
          ! only mGGA
          CALL ECCP_TAU_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1, W2, ISP, CE)
       ENDIF
    ENDIF

!=======================================================================
! non local contribution
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CE) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NT,CNL) REDUCTION(+:CE) __IF_ASYNC__

__CPU !$OMP PARALLEL SHARED(WDES1,CDIJ,W1,W2) REDUCTION(+:CE) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,NPRO,NPRO_,NI,NT,CNL,L,LP)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:W1,W2,CDIJ,WDES1) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NT,CNL) REDUCTION(+:CE) IF(OFFLOAD_ON)
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
__CPU !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             IF (WDES1%LMMAX(NT)==0) CYCLE
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             CNL=0

! currently it is better to inline manually the call to ECCP_NL
GPACC !$ACC LOOP VECTOR COLLAPSE(2) REDUCTION(+:CNL)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:CNL)
               DO L=1,WDES1%LMMAX(NT)
                  DO LP=1,WDES1%LMMAX(NT)
                     CNL=CNL+CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)*W1%CPROJ(LP+NPRO_)*GCONJG(W2%CPROJ(L+NPRO))
                  ENDDO
               ENDDO
             CE=CE+CNL
          ENDDO
__CPU !$OMP END DO
       ENDDO
    ENDDO spinor
__CPU !$OMP END PARALLEL

GPACC !$ACC KERNELS PRESENT(CE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CE) IF(OFFLOAD_ON)
    CE=GREAL(CE)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

GPACC !$ACC UPDATE SELF(CE) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC EXIT DATA DELETE(CE) __IF_ASYNC__

GPOMP !$OMP TARGET UPDATE FROM(CE) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CE) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp_')

  END SUBROUTINE ECCP_


!************************* SUBROUTINE HAMILTMU_ ************************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE S |phi>
!> where S is the overlap matrix S
!
!***********************************************************************

  SUBROUTINE HAMILTMU_( &
     HAMILTONIAN, WDES1, W1, NONLR_S, NONL_S, LATT_CUR, CDIJ, CQIJ, SV, EVALUE, ISP, WRESULT)

    USE hamil_struct_def
    USE lattice
    USE nonl_high

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct) :: NONL_S
    TYPE (latt) :: LATT_CUR
    TYPE (wavefuna) :: WRESULT
    REAL(q) :: EVALUE(:)
    OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
    RGRID :: SV(:,:)
    INTEGER :: ISP

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))

    INTEGER :: NP

    PROFILING_START('hamiltmu_')

GPACC !$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))

          IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) &
             CALL VMUJVEC(WDES1, WDES1%GRID, HAMILTONIAN%MUJVEC(1,1,4), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    ! add the action of the vector potential (if required)
    IF (ASSOCIATED(HAMILTONIAN%AVEC)) &
       CALL HAMILTMU_VEC_(HAMILTONIAN, WDES1, W1, CWORK1)

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ELSE
       ! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP), WRESULT%CW(:,NP))
             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ENDIF

    ! add mGGA and spin-current contributions (if required)
    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
       IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) THEN
          ! both mGGA and spin-current
          CALL HAMILTMU_JVEC_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, WRESULT%CW)
       ELSE
          ! only mGGA
          CALL HAMILTMU_TAU_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, WRESULT%CW)
       ENDIF
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_')

    RETURN
  END SUBROUTINE HAMILTMU_


  SUBROUTINE HAMILTMU_C_( &
     HAMILTONIAN, WDES1, W1, NONLR_S, NONL_S, LATT_CUR, CDIJ, CQIJ, SV, EVALUE, ISP, WRESULT)

    USE hamil_struct_def
    USE lattice
    USE nonl_high

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct) :: NONL_S
    TYPE (latt) :: LATT_CUR
    TYPE (wavefuna) :: WRESULT
    COMPLEX(q) :: EVALUE(:)
    OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
    RGRID :: SV(:,:)
    INTEGER :: ISP

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))

    INTEGER :: NP

    PROFILING_START('hamiltmu_c_')

GPACC !$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))

          IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) &
             CALL VMUJVEC(WDES1, WDES1%GRID, HAMILTONIAN%MUJVEC(1,1,4), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    ! add the action of the vector potential (if required)
    IF (ASSOCIATED(HAMILTONIAN%AVEC)) &
       CALL HAMILTMU_VEC_(HAMILTONIAN, WDES1, W1, CWORK1)

    IF (NONLR_S%LREAL) THEN
       ! non-local contribution in real-space
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL KINHAMIL_C(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ELSE
       ! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             SET_ACC_ASYNC_Q(NP)

             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP), WRESULT%CW(:,NP))
             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), WRESULT%CW(1,NP))
          ENDIF
       ENDDO
    ENDIF

    ! add mGGA and spin-current contributions (if required)
    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
       IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) THEN
          ! both mGGA and spin-current
          CALL HAMILTMU_JVEC_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, WRESULT%CW)
       ELSE
          ! only mGGA
          CALL HAMILTMU_TAU_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, WRESULT%CW)
       ENDIF
    ENDIF

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CWORK1) IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_c_')

    RETURN
  END SUBROUTINE HAMILTMU_C_


  SUBROUTINE HAMILT_VEC_(HAMILTONIAN, WDES1, W1, CVR)

    USE hamil_struct_def

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1

    COMPLEX(q) :: CVR(:)

    ! local variables
    COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))

    PROFILING_START('hamilt_vec_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2) IF(OFFLOAD_ON)

    CALL KINHAMIL_VEC_(WDES1, WDES1%GRID, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%AVEC(1,1), W1%CW(1), &
       CWORK1(1 __NOOFF_omp_arg(1)), CWORK2(1 __NOOFF_omp_arg(1)), &
       CVR(1))

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_vec_')

    RETURN
  END SUBROUTINE HAMILT_VEC_


  SUBROUTINE ECCP_VEC_(HAMILTONIAN, GRID, WDES1, W1, W2, CE)

    USE hamil_struct_def
    USE mgrid_struct_def

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (grid_3d) :: GRID
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1, W2

    COMPLEX(q) :: CE

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))

    INTEGER :: ISPINOR, M, MM, NGVECTOR

    PROFILING_START('eccp_vec_')

    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CWORK1=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALL KINHAMIL_VEC_(WDES1, WDES1%GRID, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%AVEC(1,1), W1%CW(1), &
       CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), &
       CWORK1(1))

    DO ISPINOR=0,WDES1%NRSPINORS-1
       CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), CWORK1(1+ISPINOR*WDES1%GRID%MPLWV), CWORK2(1+ISPINOR*WDES1%NGVECTOR __NOOFF_omp_arg(1)), GRID, .FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(CWORK2,W2,CE) PRIVATE(MM) REDUCTION(+:CE) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK2,W2,CE) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CE) IF(OFFLOAD_ON)

       DO M=1,NGVECTOR
          MM =M+ISPINOR*NGVECTOR
! test_
!!          CE=CE+CWORK2(MM __NOOFF_omp_arg(1))*CONJG(W2%CW(MM)) + W1%CW(MM)*CONJG(W2%CW(MM)) * WDES1%DATAKE(M,ISPINOR+1)
          CE=CE+CWORK2(MM __NOOFF_omp_arg(1))*CONJG(W2%CW(MM))
! test_
       ENDDO
    ENDDO

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp_vec_')

  END SUBROUTINE ECCP_VEC_


  SUBROUTINE HAMILTMU_VEC_(HAMILTONIAN, WDES1, W1, CVR)

    USE hamil_struct_def

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)

    COMPLEX(q) :: CVR(:,:)

    ! local variables
    INTEGER :: NP

#ifdef ACC_OFFLOAD
    COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV),     SIZE(W1))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%GRID%MPLWV), SIZE(W1))
#else
    COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))
#endif

    PROFILING_START('hamiltmu_vec_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2) IF(OFFLOAD_ON)

    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL KINHAMIL_VEC_(WDES1, WDES1%GRID, &
             WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
             HAMILTONIAN%AVEC(1,1), W1(NP)%CW(1), &
#ifdef ACC_OFFLOAD
             CWORK1(1,NP), CWORK2(1,NP), &
#else
             CWORK1(1 __NOOFF_omp_arg(1)), CWORK2(1 __NOOFF_omp_arg(1)), &
#endif
             CVR(1,NP))
       ENDIF
    ENDDO

GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_vec_')

    RETURN
  END SUBROUTINE HAMILTMU_VEC_


  SUBROUTINE HAMILT_TAU_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, CH)

    USE hamil_struct_def
    USE lattice

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    TYPE (latt) :: LATT_CUR

    COMPLEX(q) :: CH(:)

    INTEGER :: ISP

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))

    PROFILING_START('hamilt_tau_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    CALL KINHAMIL_TAU_(WDES1, WDES1%GRID, LATT_CUR, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%MU(1,ISP), CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), &
       W1%CW(1), CH(1))

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_tau_')

    RETURN
  END SUBROUTINE HAMILT_TAU_


  SUBROUTINE ECCP_TAU_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1, W2, ISP, CE)

    USE hamil_struct_def
    USE mgrid_struct_def
    USE poscar_struct_def

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (grid_3d) :: GRID
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1, W2
    TYPE (latt) :: LATT_CUR

    INTEGER :: ISP

    COMPLEX(q) :: CE

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
    COMPLEX(q) :: CH(WDES1%NRPLWV)

    INTEGER :: ISPINOR, M, MM, NGVECTOR

    PROFILING_START('eccp_tau_')

    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CH,CWORK1,CWORK2,CWORK3) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CH,CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CH) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CH) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CH=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    CALL KINHAMIL_TAU_(WDES1, WDES1%GRID, LATT_CUR, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%MU(1,ISP), CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), &
       W1%CW(1), CH(1))

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES1,CH,W2,CE) PRIVATE(MM) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:WDES1,CH,W2,CE) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CE) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO COLLAPSE(2) SHARED(WDES1,CH,W2) PRIVATE(MM) REDUCTION(+:CE)
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CE=CE+CH(MM)*CONJG(W2%CW(MM))
       ENDDO
    ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA DELETE(CH,CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CH,CWORK1,CWORK2,CWORK3,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp_tau_')

    RETURN
  END SUBROUTINE ECCP_TAU_


  SUBROUTINE HAMILTMU_TAU_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, CH)

    USE hamil_struct_def
    USE lattice

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (latt) :: LATT_CUR

    COMPLEX(q) :: CH(:,:)

    INTEGER :: ISP

    ! local variables
    INTEGER :: NP

#ifdef ACC_OFFLOAD
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
#endif

    PROFILING_START('hamiltmu_tau_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL KINHAMIL_TAU_(WDES1, WDES1%GRID, LATT_CUR, &
             WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
             HAMILTONIAN%MU(1,ISP), &
#ifdef ACC_OFFLOAD
             CWORK1(1,NP), CWORK2(1,NP), CWORK3(1,NP), &
#else
             CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), &
#endif
             W1(NP)%CW(1), CH(1,NP))
       ENDIF
    ENDDO

GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_tau_')

    RETURN
  END SUBROUTINE HAMILTMU_TAU_


  SUBROUTINE HAMILT_JVEC_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, CH)

    USE hamil_struct_def
    USE lattice

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    TYPE (latt) :: LATT_CUR

    COMPLEX(q) :: CH(:)

    INTEGER :: ISP

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
    COMPLEX(q) :: CWORK4(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))

    PROFILING_START('hamilt_jvec_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3,CWORK4) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3,CWORK4) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    CALL KINHAMIL_JVEC_(WDES1, WDES1%GRID, LATT_CUR, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%MU(1,ISP), HAMILTONIAN%MUJVEC(1,1,1), &
       CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), CWORK4(1 __NOOFF_omp_arg(1)), &
       W1%CW(1), CH(1))

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR%B,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_jvec_')

    RETURN
  END SUBROUTINE HAMILT_JVEC_


  SUBROUTINE ECCP_JVEC_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1, W2, ISP, CE)

    USE hamil_struct_def
    USE mgrid_struct_def
    USE poscar_struct_def

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (grid_3d) :: GRID
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1, W2
    TYPE (latt) :: LATT_CUR

    INTEGER :: ISP

    COMPLEX(q) :: CE

    ! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
    COMPLEX(q) :: CWORK4(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))
    COMPLEX(q) :: CH(WDES1%NRPLWV)

    INTEGER :: ISPINOR, M, MM, NGVECTOR

    PROFILING_START('eccp_jvec_')

    NGVECTOR=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CH,CWORK1,CWORK2,CWORK3,CWORK4) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CH,CWORK1,CWORK2,CWORK3,CWORK4) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CH) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CH) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CH=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    CALL KINHAMIL_JVEC_(WDES1, WDES1%GRID, LATT_CUR, &
       WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
       HAMILTONIAN%MU(1,ISP), HAMILTONIAN%MUJVEC(1,1,1), &
       CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), CWORK4(1 __NOOFF_omp_arg(1)), &
       W1%CW(1), CH(1))

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES1,CH,W2,CE) PRIVATE(MM) REDUCTION(+:CE) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:WDES1,CH,W2,CE) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:CE) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO COLLAPSE(2) SHARED(WDES1,CH,W2) PRIVATE(MM) REDUCTION(+:CE)
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CE=CE+CH(MM)*CONJG(W2%CW(MM))
       ENDDO
    ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA DELETE(CH,CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR%B,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CH,CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('eccp_jvec_')

    RETURN
  END SUBROUTINE ECCP_JVEC_


  SUBROUTINE HAMILTMU_JVEC_(HAMILTONIAN, WDES1, W1, LATT_CUR, ISP, CH)

    USE hamil_struct_def
    USE lattice

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (latt) :: LATT_CUR

    COMPLEX(q) :: CH(:,:)

    INTEGER :: ISP

    ! local variables
    INTEGER :: NP

#ifdef ACC_OFFLOAD
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
    COMPLEX(q) :: CWORK4(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __NOOFF_omp_arg(3))
    COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(6))
    COMPLEX(q) :: CWORK4(__align_DIM(WDES1%GRID%MPLWV) __NOOFF_omp_arg(3))
#endif

    PROFILING_START('hamiltmu_jvec_')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3,CWORK4) IF(OFFLOAD_ON)
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3,CWORK4) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL KINHAMIL_JVEC_(WDES1, WDES1%GRID, LATT_CUR, &
             WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), &
             HAMILTONIAN%MU(1,ISP), HAMILTONIAN%MUJVEC(1,1,1), &
#ifdef ACC_OFFLOAD
             CWORK1(1,NP), CWORK2(1,NP), CWORK3(1,NP), CWORK4(1,NP), &
#else
             CWORK1(1), CWORK2(1 __NOOFF_omp_arg(1)), CWORK3(1 __NOOFF_omp_arg(1)), CWORK4(1 __NOOFF_omp_arg(1)), &
#endif
             W1(NP)%CW(1), CH(1,NP))
       ENDIF
    ENDDO

GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR%B,LATT_CUR) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3,CWORK4,LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

    PROFILING_STOP('hamiltmu_jvec_')

    RETURN
  END SUBROUTINE HAMILTMU_JVEC_


!************************* SUBROUTINE HAMILT_LOCAL *********************
!
!> this subroutine calculates the local and kinetic energy part of H
!> acting onto a wavefuntion
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> LADD=.TRUE.
!> ~~~
!>  CH = CH + SV* W1 + T * W1
!> ~~~
!> LADD=.FALSE.
!>  CH = SV* W1 + T * W1
!> LKIN determines whether the kinetic energy contribution is added or not
!> the default is add the kinetic energy contribution
!
!***********************************************************************

  SUBROUTINE HAMILT_LOCAL_(HAMILTONIAN, LATT_CUR, W1, SV, ISP, CH, LADD, LKIN)

    USE hamil_struct_def
    USE lattice

    IMPLICIT NONE

    TYPE (ham_handle) :: HAMILTONIAN
    TYPE (latt) :: LATT_CUR
    TYPE (wavefun1) :: W1
    RGRID :: SV(:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CH(:)
    LOGICAL :: LADD
    LOGICAL, OPTIONAL :: LKIN

    ! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: CWORK1
#endif

    PROFILING_START('hamilt_local_')

GPACC !$ACC DATA CREATE(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

    ! calculate the local contribution (result in CWORK1)
    CALL VHAMIL(W1%WDES1, W1%WDES1%GRID, SV(1,ISP), W1%CR(1), CWORK1(1))

    IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) &
       CALL VMUJVEC(W1%WDES1, W1%WDES1%GRID, HAMILTONIAN%MUJVEC(1,1,4), W1%CR(1), CWORK1(1))

    ! add the action of the vector potential (if required)
    IF (ASSOCIATED(HAMILTONIAN%AVEC)) &
       CALL HAMILT_VEC_(HAMILTONIAN, W1%WDES1, W1, CWORK1)

    IF (.NOT. PRESENT(LKIN)) THEN
       ! bring CWORK1 to reciprocal space and add T|W1> (result in CH)
       CALL KINHAMIL(W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
           W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))
    ELSE
       IF (LKIN) THEN
          ! bring CWORK1 to reciprocal space and add T|W1> (result in CH)
          CALL KINHAMIL(W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
              W1%WDES1%DATAKE(1,1), 0.0_q, W1%CW(1), CH(1))
       ELSE
          ! bring CWORK1 to reciprocal space (result in CH)
          CALL FFTHAMIL(W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
                                    0.0_q, W1%CW(1), CH(1))
       ENDIF
    ENDIF

    ! add mGGA and spin-current contributions (if required)
    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
       IF (ASSOCIATED(HAMILTONIAN%MUJVEC)) THEN
          ! both mGGA and spin-current
          CALL HAMILT_JVEC_(HAMILTONIAN, W1%WDES1, W1, LATT_CUR, ISP, CH)
       ELSE
          ! only mGGA
          CALL HAMILT_TAU_(HAMILTONIAN, W1%WDES1, W1, LATT_CUR, ISP, CH)
       ENDIF
    ENDIF

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC END DATA
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

    PROFILING_STOP('hamilt_local_')

    RETURN
  END SUBROUTINE HAMILT_LOCAL_
#endif // oldhamil

END MODULE hamil

!************************* SUBROUTINE VHAMIL ***************************
!
!> low level F77 to calculate the product of the local potential
!> times a wavefunction stored in CR and returns the result in CVR
!> ~~~
!>   CVR= SV * CR
!> ~~~
!> this low level routine should be called only from within
!> hamil.F
!
!***********************************************************************

   SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER VHAMIL
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1

     RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
     ! local variables
     REAL(q) RINPLW
     INTEGER ISPINOR, ISPINOR_, M, MM, MM_
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CVR: __align_BTS
#endif

     PROFILING_START('vhamil')

     RINPLW=1._q/GRID%NPLWV

     IF (WDES1%NRSPINORS==1) THEN
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(GRID,CVR,SV,CR) &
GPACC !$ACC __IF_ASYNC__ DEFAULT(PRESENT)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,CVR,SV,CR) &
GPOMP !$OMP IF(OFFLOAD_ON)

__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(GRID, CVR, SV, CR, RINPLW)

         DO M=1,GRID%RL%NP
            CVR(M)= SV(M,1) *CR(M)*RINPLW
         ENDDO

__CPU !$OMP END PARALLEL DO
     ELSE
!!!        CVR(1:GRID%MPLWV*2)=0
!!!!$OMP PARALLEL &
!!!!$OMP SHARED(GRID,CVR,SV,CR,RINPLW) &
!!!!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
!!!        DO ISPINOR =0,1
!!!           DO ISPINOR_=0,1
!!!!DIR$ IVDEP
!!!!OCL NOVREL
!!!!$OMP DO SCHEDULE(STATIC)
!!!              DO M=1,GRID%RL%NP
!!!                 MM =M+ISPINOR *GRID%MPLWV
!!!                 MM_=M+ISPINOR_*GRID%MPLWV
!!!                 CVR(MM)= CVR(MM)+ SV(M,1+ISPINOR_+2*ISPINOR) *CR(MM_)*RINPLW
!!!              ENDDO
!!!!$OMP END DO
!!!           ENDDO
!!!        ENDDO
!!!!$OMP END PARALLEL

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,CR) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRID,CVR,CR) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(GRID,CVR,SV,CR,RINPLW) PRIVATE(ISPINOR,M,MM)
        DO ISPINOR =0,1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)
            DO M=1,GRID%RL%NP
               MM =M+ISPINOR *GRID%MPLWV
               CVR(MM)=SV(M,1+2*ISPINOR) *CR(M)*RINPLW
            ENDDO
__CPU !$OMP END DO
__GPU   ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRID,CVR,SV,CR) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

__GPU   DO ISPINOR =0,1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)

           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              CVR(MM)= CVR(MM)+ SV(M,2+2*ISPINOR) *CR(M+GRID%MPLWV)*RINPLW
           ENDDO

__CPU !$OMP END DO
        ENDDO
__CPU !$OMP END PARALLEL
     ENDIF

     PROFILING_STOP('vhamil')

   END SUBROUTINE VHAMIL

   SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER VHAMIL_TRACE
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1

     COMPLEX(q) :: CR( GRID%MPLWV*WDES1%NRSPINORS), &
                   CVR(GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: SV(GRID%MPLWV)                  ! local potential

     REAL(q) :: WEIGHT

     ! local variables
     INTEGER :: M, MM

     PROFILING_START('vhamil_trace')

__CPU !$OMP PARALLEL SHARED(GRID,CVR,SV,CR,WEIGHT) PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)
GPACC !$ACC PARALLEL LOOP PRESENT(GRID,CVR,SV,CR) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,CVR,SV,CR) IF(OFFLOAD_ON)
     DO M=1,GRID%RL%NP
        CVR(M)= CVR(M)+SV(M) *CR(M)*WEIGHT
     ENDDO
__CPU !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)
GPACC !$ACC PARALLEL LOOP PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,CVR,SV,CR) IF(OFFLOAD_ON)
        DO M=1,GRID%RL%NP
           MM =M+GRID%MPLWV
           CVR(MM)= CVR(MM)+ SV(M) *CR(MM)*WEIGHT
        ENDDO
__CPU !$OMP END DO
     ENDIF
__CPU !$OMP END PARALLEL

     PROFILING_STOP('vhamil_trace')

   END SUBROUTINE VHAMIL_TRACE

   SUBROUTINE VHAMIL_TRACE_MU(WDES1, GRID, SV, LD, CR, WXI, WEIGHT, NSIM)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER VHAMIL_TRACE_MU
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: WXI(NSIM)

     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: SV(LD,NSIM)                     ! local potential

     REAL(q) :: WEIGHT
     INTEGER :: LD,NSIM

     ! local variables
     INTEGER ISPINOR, M, MM, N, IBLOCK
     INTEGER, PARAMETER :: BLOCKSIZE=128

     PROFILING_START('vhamil_trace_mu')

!DIR$ IVDEP
!OCL NOVREL
     ! TODO: Check OpenMP performance
GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(WXI,SV,CR,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WXI,SV,CR,GRID) IF(OFFLOAD_ON)
     DO N=1,NSIM
        DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CR(IBLOCK:MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))
#endif
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              WXI(N)%CR(M)= WXI(N)%CR(M)+SV(M,N) *CR(M)*WEIGHT
           ENDDO
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
      ! TODO: Check OpenMP performance
GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(WXI,SV,CR,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WXI,SV,CR,GRID) IF(OFFLOAD_ON)
        DO N=1,NSIM
           DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CR(GRID%MPLWV+IBLOCK:GRID%MPLWV+MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))
#endif
GPACC !$ACC LOOP VECTOR PRIVATE(MM)
GPOMP !$OMP PARALLEL DO PRIVATE(MM) IF(OFFLOAD_ON)
              DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                 MM =M+GRID%MPLWV
                 WXI(N)%CR(MM)= WXI(N)%CR(MM)+SV(M,N) *CR(MM)*WEIGHT
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     PROFILING_STOP('vhamil_trace_mu')

   END SUBROUTINE VHAMIL_TRACE_MU
#ifdef OFFLOADING
   SUBROUTINE VHAMIL_TRACE_MU_ACC(WDES1, GRID, SV, LD, CR, WXI, WEIGHT, NSIM)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: WXI(NSIM)

     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: SV(LD,NSIM)                     ! local potential

     REAL(q) :: WEIGHT
     INTEGER :: LD,NSIM

     ! local variables
     INTEGER ISPINOR, M, MM, N, IBLOCK
     INTEGER, PARAMETER :: BLOCKSIZE=512
     COMPLEX(q) :: CTMP(BLOCKSIZE)

     PROFILING_START('vhamil_trace_mu_acc')

!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(WXI,SV,CR,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(CTMP) &
! for intel ifx limit of threads is improving performance for reduction
#ifdef IFXWAR6
GPOMP !$OMP THREAD_LIMIT(64) IF(OFFLOAD_ON)
#else
GPOMP !$OMP IF(OFFLOAD_ON)
#endif
     DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
        DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
           CTMP(M-IBLOCK+1) = CR(M)
        ENDDO
GPACC !$ACC LOOP WORKER
        DO N=1,NSIM
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              WXI(N)%CR(M)= WXI(N)%CR(M)+SV(M,N) *CTMP(M-IBLOCK+1)*WEIGHT
           ENDDO
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(WXI,SV,CR,GRID) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WXI,SV,CR,GRID) PRIVATE(CTMP) &
! for intel ifx limit of threads is improving performance for reduction
#ifdef IFXWAR6
GPOMP !$OMP THREAD_LIMIT(64) IF(OFFLOAD_ON)
#else
GPOMP !$OMP IF(OFFLOAD_ON)
#endif
        DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR PRIVATE(MM)
GPOMP !$OMP PARALLEL DO PRIVATE(MM) IF(OFFLOAD_ON)
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              MM =M+GRID%MPLWV
              CTMP(M-IBLOCK+1)= CR(MM)
           ENDDO
GPACC !$ACC LOOP WORKER
           DO N=1,NSIM
GPACC !$ACC LOOP VECTOR PRIVATE(MM)
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
              DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                 MM =M+GRID%MPLWV
                 WXI(N)%CR(MM)= WXI(N)%CR(MM)+SV(M,N) *CTMP(M-IBLOCK+1)*WEIGHT
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     PROFILING_STOP('vhamil_trace_mu_acc')

   END SUBROUTINE VHAMIL_TRACE_MU_ACC
#endif


!************************* SUBROUTINE PW_CHARGE ************************
!
!> low level F77 to calculate charge density from two plane wave states
!> and add the weighted spinor like density to an array
!> this low level routine should be called only from within
!> hamil.F
!> ~~~
!>  CHARGE = CR1*CONJG(CR2) * WEIGHT + CHARGE
!> ~~~
!
!***********************************************************************

!
! this version yields an RGRID array spinor density
!

   SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     RGRID   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
     ! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_,NP,MPLWV

     NP = WDES1%GRID%RL%NP
     MPLWV = WDES1%GRID%MPLWV

     PROFILING_START('pw_charge')

     IF (WDES1%NRSPINORS==1) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(CHARGE,CR1,CR2) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CHARGE,CR1,CR2) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREL
__CPU !$OMP PARALLEL DO PRIVATE(M) SHARED(WDES1,CHARGE,CR1,CR2,WEIGHT)
        DO M=1,NP
           CHARGE(M,1)=CHARGE(M,1)+ &
                CONJG(CR2(M))*CR1(M)*WEIGHT
        ENDDO
__CPU !$OMP END PARALLEL DO
     ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CHARGE,CR1,CR2) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
GPOMP !$OMP MAP(PRESENT,TO:CHARGE,CR1,CR2) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(WDES1,CHARGE,CR1,CR2,WEIGHT) PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP DO
              DO M=1,NP
                 MM =M+ISPINOR *MPLWV
                 MM_=M+ISPINOR_*MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ &
                      CONJG(CR2(MM_))*CR1(MM)*WEIGHT
              ENDDO
__CPU !$OMP END DO NOWAIT
           ENDDO
        ENDDO spinor
__CPU !$OMP END PARALLEL
     ENDIF

     PROFILING_STOP('pw_charge')

   END SUBROUTINE PW_CHARGE

!
!> this version yields a GDEF charge array, spinor density
!> (used by exchange routines)
!
   SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2 )
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     GDEF   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
     ! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_

     IF (WDES1%NRSPINORS==1) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,1)=CONJG(CR2(M))*CR1(M)
        ENDDO
     ELSE
        CHARGE=0
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
              DO M=1,WDES1%GRID%RL%NP
                 MM =M+ISPINOR *WDES1%GRID%MPLWV
                 MM_=M+ISPINOR_*WDES1%GRID%MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ &
                      CONJG(CR2(MM_))*CR1(MM)
              ENDDO
           ENDDO
        ENDDO spinor
     ENDIF

   END SUBROUTINE PW_CHARGE_CMPLX

!
!> this version traces of spinors to yield only the density
!> results array is GDEF (used by exchange routines)
!
   SUBROUTINE PW_CHARGE_TRACE(WDES1, CHARGE, CR1, CR2)

__GPU USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1

     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: CHARGE(WDES1%GRID%MPLWV)

     ! local variables
     REAL(q) RINPLW
     INTEGER M,MM

     PROFILING_START('pw_charge_trace')

NOACC !$OMP PARALLEL SHARED(WDES1,CHARGE,CR1,CR2) PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
NOACC !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) __IF_ASYNC__
     DO M=1,WDES1%GRID%RL%NP
        CHARGE(M)=CONJG(CR2(M))*CR1(M)
     ENDDO
NOACC !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
NOACC !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) PRIVATE(MM) __IF_ASYNC__
        DO M=1,WDES1%GRID%RL%NP
           MM =M+WDES1%GRID%MPLWV
           CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
        ENDDO
NOACC !$OMP END DO
     ENDIF
NOACC !$OMP END PARALLEL

     PROFILING_STOP('pw_charge_trace')

   END SUBROUTINE PW_CHARGE_TRACE

   SUBROUTINE PW_CHARGE_TRACE_MU(WDES1, CHARGE, LD, W1, CR2, NSIM)

__GPU USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: W1(NSIM)

     COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: CHARGE(LD,NSIM)

     INTEGER :: LD,NSIM

     PROFILING_START('pw_charge_trace_mu')

#ifdef OFFLOADING
     IF (OFFLOAD_ON) THEN
        CALL PW_CHARGE_TRACE_MU_ACC(WDES1, CHARGE(1,1), LD, W1, CR2(1), NSIM)
        PROFILING_STOP('pw_charge_trace_mu')
        RETURN
     ENDIF
#endif
     CALL PW_CHARGE_TRACE_MU_(WDES1, CHARGE(1,1), LD, W1, CR2(1), NSIM)

     PROFILING_STOP('pw_charge_trace_mu')

     RETURN
   END SUBROUTINE PW_CHARGE_TRACE_MU

   SUBROUTINE PW_CHARGE_TRACE_MU_(WDES1, CHARGE, LD, W1, CR2, NSIM)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: W1(NSIM)

     COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: CHARGE(LD,NSIM)

     INTEGER :: LD,NSIM

     ! local variables
     INTEGER :: N,M,MM

     PROFILING_START('pw_charge_trace_mu_')

!DIR$ IVDEP
!OCL NOVREL
     DO N=1,NSIM
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,N)=CONJG(CR2(M))*W1(N)%CR(M)
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO N=1,NSIM
           DO M=1,WDES1%GRID%RL%NP
              MM =M+WDES1%GRID%MPLWV
              CHARGE(M,N)=CHARGE(M,N)+CONJG(CR2(MM))*W1(N)%CR(MM)
           ENDDO
        ENDDO
     ENDIF

     PROFILING_STOP('pw_charge_trace_mu_')

   END SUBROUTINE PW_CHARGE_TRACE_MU_

   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU(WDES1, CHARGE, LD, W1, CR2, NSIM)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(NSIM)

    COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    GDEF       :: CHARGE(LD,NSIM)

    INTEGER :: LD,NSIM

    ! local variables
    INTEGER :: N,M,MM,IBLOCK
    INTEGER, PARAMETER :: BLOCKSIZE=128

    PROFILING_START('pw_charge_trace_no_conjg_mu')

!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
    DO N=1,NSIM
       DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
!$ACC CACHE(CR2(IBLOCK:MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))
#endif
!$ACC LOOP VECTOR
          DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
             CHARGE(M,N)=CR2(M)*W1(N)%CR(M)
          ENDDO
       ENDDO
    ENDDO
    IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
       DO N=1,NSIM
          DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
!$ACC CACHE(CR2(WDES1%GRID%MPLWV+IBLOCK:WDES1%GRID%MPLWV+MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))
#endif
!$ACC LOOP VECTOR PRIVATE(MM)
             DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                MM =M+WDES1%GRID%MPLWV
                CHARGE(M,N)=CHARGE(M,N)+CR2(MM)*W1(N)%CR(MM)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    PROFILING_STOP('pw_charge_trace_no_conjg_mu')

  END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU

#ifdef OFFLOADING
   SUBROUTINE PW_CHARGE_TRACE_MU_ACC(WDES1, CHARGE, LD, W1, CR2, NSIM)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER PW_CHARGE_TRACE_MU_ACC
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: W1(NSIM)

     COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     GDEF       :: CHARGE(LD,NSIM)

     INTEGER :: LD,NSIM

     ! local variables
     INTEGER :: N,M,MM,IBLOCK
     INTEGER, PARAMETER :: BLOCKSIZE=512
     COMPLEX(q) :: CTMP(BLOCKSIZE)

     PROFILING_START('pw_charge_trace_mu_acc')

!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:CHARGE,CR2,W1) PRIVATE(CTMP) &
! for intel ifx limit of threads is improving performance for reduction
#ifdef IFXWAR6
GPOMP !$OMP THREAD_LIMIT(64) IF(OFFLOAD_ON)
#else
GPOMP !$OMP IF(OFFLOAD_ON)
#endif
     DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR
GPOMP !$OMP PARALLEL IF(OFFLOAD_ON)
GPOMP !$OMP DO
        DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
           CTMP(M-IBLOCK+1)=CONJG(CR2(M))
        ENDDO

GPACC !$ACC LOOP WORKER
        DO N=1,NSIM
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP DO
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              CHARGE(M,N)=CTMP(M-IBLOCK+1)*W1(N)%CR(M)
           ENDDO
        ENDDO
GPOMP !$OMP END PARALLEL
     ENDDO

     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:CHARGE,CR2,W1) PRIVATE(CTMP) &
! for intel ifx limit of threads is improving performance for reduction
#ifdef IFXWAR6
GPOMP !$OMP THREAD_LIMIT(64) IF(OFFLOAD_ON)
#else
GPOMP !$OMP IF(OFFLOAD_ON)
#endif
        DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR PRIVATE(MM)
GPOMP !$OMP PARALLEL PRIVATE(MM) IF(OFFLOAD_ON)
GPOMP !$OMP DO
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              MM =M+WDES1%GRID%MPLWV
              CTMP(M-IBLOCK+1)=CONJG(CR2(MM))
           ENDDO

GPACC !$ACC LOOP WORKER
           DO N=1,NSIM
GPACC !$ACC LOOP VECTOR PRIVATE(MM)
GPOMP !$OMP DO PRIVATE(MM)
              DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                 MM =M+WDES1%GRID%MPLWV
                 CHARGE(M,N)=CHARGE(M,N)+CTMP(M-IBLOCK+1)*W1(N)%CR(MM)
              ENDDO
           ENDDO
GPOMP !$OMP END PARALLEL
        ENDDO
     ENDIF

     PROFILING_STOP('pw_charge_trace_mu_acc')

   END SUBROUTINE PW_CHARGE_TRACE_MU_ACC

   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU_ACC(WDES1, CHARGE, LD, W1, CR2, NSIM)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER PW_CHARGE_TRACE_NO_CONJG_MU_ACC
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(NSIM)

    COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    GDEF       :: CHARGE(LD,NSIM)

    INTEGER :: LD,NSIM

    ! local variables
    INTEGER :: N,M,MM,IBLOCK
    INTEGER, PARAMETER :: BLOCKSIZE=512
    COMPLEX(q) :: CTMP(BLOCKSIZE)

    PROFILING_START('pw_charge_trace_no_conjg_mu_acc')

!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:CHARGE,CR2,W1) PRIVATE(CTMP) IF(OFFLOAD_ON)
    DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
       DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
          CTMP(M-IBLOCK+1)=CR2(M)
       ENDDO

GPACC !$ACC LOOP WORKER
GPOMP !$OMP PARALLEL DO COLLAPSE(2) IF(OFFLOAD_ON)
       DO N=1,NSIM
GPACC !$ACC LOOP VECTOR
          DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
             CHARGE(M,N)=CTMP(M-IBLOCK+1)*W1(N)%CR(M)
          ENDDO
       ENDDO
    ENDDO
    IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
GPACC !$ACC PARALLEL LOOP GANG PRIVATE(CTMP) PRESENT(CHARGE,CR2,W1) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:CHARGE,CR2,W1) PRIVATE(CTMP) IF(OFFLOAD_ON)
       DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE
#ifndef CRAYHIP
GPACC !$ACC CACHE(CTMP(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR PRIVATE(MM)
GPOMP !$OMP PARALLEL DO PRIVATE(MM) IF(OFFLOAD_ON)
          DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
             MM =M+WDES1%GRID%MPLWV
             CTMP(M-IBLOCK+1)=CR2(MM)
          ENDDO

GPACC !$ACC LOOP WORKER
GPOMP !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(MM) IF(OFFLOAD_ON)
          DO N=1,NSIM
GPACC !$ACC LOOP VECTOR PRIVATE(MM)
             DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                MM =M+WDES1%GRID%MPLWV
                CHARGE(M,N)=CHARGE(M,N)+CTMP(M-IBLOCK+1)*W1(N)%CR(MM)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    PROFILING_STOP('pw_charge_trace_no_conjg_mu_acc')

  END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU_ACC


#endif

   SUBROUTINE PW_CHARGE_TRACE_MU2(WDES1, CHARGE, LDC, CR1, LD1, CR2, LD2, NSIM)

__GPU USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1

     GDEF       :: CHARGE(LDC,NSIM)
     COMPLEX(q) :: CR1(LD1,NSIM),CR2(LD1,NSIM)

     INTEGER :: LDC,LD1,LD2,NSIM

     ! local variables
     INTEGER :: N

     PROFILING_START('pw_charge_trace_mu2')

#ifdef ACC_OFFLOAD
     IF (OFFLOAD_ON) THEN
        CALL PW_CHARGE_TRACE_MU2_ACC(WDES1, CHARGE(1,1), LDC, CR1(1,1), LD1, CR2(1,1), LD2, NSIM)
        PROFILING_STOP('pw_charge_trace_mu2')
        RETURN
     ENDIF
#endif
     DO N=1,NSIM
        CALL PW_CHARGE_TRACE(WDES1, CHARGE(1,N), CR1(1,N), CR2(1,N))
     ENDDO

     PROFILING_STOP('pw_charge_trace_mu2')

     RETURN
   END SUBROUTINE PW_CHARGE_TRACE_MU2

#ifdef ACC_OFFLOAD
   SUBROUTINE PW_CHARGE_TRACE_MU2_ACC(WDES1, CHARGE, LDC, CR1, LD1, CR2, LD2, NSIM)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1

     GDEF       :: CHARGE(LDC,NSIM)
     COMPLEX(q) :: CR1(LD1,NSIM),CR2(LD2,NSIM)

     INTEGER :: LDC,LD1,LD2,NSIM

     ! local variables
     INTEGER :: N,M,MM

     PROFILING_START('pw_charge_trace_mu2_acc')

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR1,CR2,CHARGE) __IF_ASYNC__
     DO N=1,NSIM
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,N)=CR1(M,N)*CONJG(CR2(M,N))
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR1,CR2,CHARGE) PRIVATE(MM) __IF_ASYNC__
        DO N=1,NSIM
           DO M=1,WDES1%GRID%RL%NP
              MM=M+WDES1%GRID%MPLWV
              CHARGE(M,N)=CR1(MM,N)*CONJG(CR2(MM,N))
           ENDDO
        ENDDO
     ENDIF

     PROFILING_STOP('pw_charge_trace_mu2_acc')

   END SUBROUTINE PW_CHARGE_TRACE_MU2_ACC
#endif
!
!> this version traces of spinors to yield only the density
!> result array as well as orbital input is GDEF

   SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     GDEF :: CHARGE(WDES1%GRID%MPLWV)
     GDEF :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     ! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=GCONJG(CR2(M))*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+GCONJG(CR2(MM))*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_GDEF

!
!> this version traces of spinors to yield only the density
!> result array with a REAL(q) layout (used by chi_super.F)
!
   SUBROUTINE PW_CHARGE_TRACE_REAL(WDES1,  CHARGE, CR1, CR2)

__GPU USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q)    :: CHARGE(WDES1%GRID%MPLWV)

     ! local variables
     REAL(q) RINPLW
     INTEGER M,MM

     PROFILING_START('pw_charge_trace_real')

NOACC !$OMP PARALLEL SHARED(WDES1,CHARGE,CR1,CR2) !$OMP PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
NOACC !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) __IF_ASYNC__
     DO M=1,WDES1%GRID%RL%NP
        CHARGE(M)=CONJG(CR2(M))*CR1(M)
     ENDDO
NOACC !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
NOACC !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) PRIVATE(MM) __IF_ASYNC__
        DO M=1,WDES1%GRID%RL%NP
           MM =M+WDES1%GRID%MPLWV
           CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
        ENDDO
NOACC !$OMP END DO
     ENDIF
NOACC !$OMP END PARALLEL

     PROFILING_STOP('pw_charge_trace_real')

   END SUBROUTINE PW_CHARGE_TRACE_REAL


   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     GDEF   CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     ! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CR2(M)*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CR2(MM)*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG


!************************* SUBROUTINE KINHAMIL *************************
!
!> low level F77 to calculate
!> the kinetic energy part of a Hamiltonian
!> possibly substract the wavefunction times an eigenvalue
!> and adds a component after FFT to real space
!> this low level routine should be called only from within
!> hamil.F
!> ~~~
!> LADD=.TRUE.
!> ~~~
!> ~~~
!>  CH= CH + FFT(CVR) + DATAKE*CW - EVALUE*CW
!> ~~~
!> ~~~
!> LADD=.FALSE.
!> ~~~
!> ~~~
!>  CH= FFT(CVR) + DATAKE*CW - EVALUE*CW
!> ~~~
!>
!> FFTHAMIL exclude the kinetic energy contribution
!
!***********************************************************************

   SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CW, CH)

#ifdef CRAYWAR4
! in CCE20 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL
#endif
__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
   ! local
     INTEGER ISPINOR, M, MM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CVR: __align_BTS
#endif

     PROFILING_START('kinhamil')

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)

GPACC !$ACC PARALLEL LOOP PRESENT(CH,CW,DATAKE,WDES1) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CH,CW,DATAKE,WDES1) PRIVATE(MM) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREL
__CPU !$OMP PARALLEL DO SHARED(WDES1,CH,ISPINOR,DATAKE,EVALUE) PRIVATE (M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)+CW(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
__CPU !$OMP END PARALLEL DO
     ENDDO

     PROFILING_STOP('kinhamil')

   END SUBROUTINE KINHAMIL


   SUBROUTINE FFTHAMIL( WDES1, GRID, CVR, LADD, EVALUE, CW, CH)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
   ! local
     INTEGER ISPINOR, M, MM

     PROFILING_START('ffthamil')

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
GPACC !$ACC PARALLEL LOOP PRESENT(CH,CW) PRIVATE(MM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CH,CW) PRIVATE(MM) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREL
__CPU !$OMP PARALLEL DO SHARED(WDES1,CH,ISPINOR,EVALUE) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CW(MM)*EVALUE
        ENDDO
__CPU !$OMP END PARALLEL DO
     ENDDO

     PROFILING_STOP('ffthamil')

   END SUBROUTINE FFTHAMIL


!************************* SUBROUTINE KINHAMIL_VEC *********************
!
!> low level F77 to calculate
!> similar to #KINHAMIL
!> but adds the action of the vector potential on the wavefunction
!>
!> ~~~
!> -(e A p) e(iGr)
!>  |e|A p  e(iGr) = (|e| A hbar) -i nabla e(iGr)
!> (|e| A hbar) G e(iGr)
!> ~~~
!>
!> the coefficients |e| A hbar are included in the stored AVEC
!>
!> @details @ref openmp :
!> the workarrays CWORK1 and CWORK2 acquire a second dimension (size=3).
!> The additional workspace is used to reorder the routine. First fill
!> all three components of CWORK1, then do all FFTs, then add the
!> relevant contributions of all components of CWORK2 to CVR inside
!> a single loop.
!> This improves vectorization of the loops over the elements of these
!> arrays. Aforementioned loops are distributed over all available
!> OpenMP threads.
!
!***********************************************************************

   SUBROUTINE KINHAMIL_VEC( WDES1, GRID, CVR,  LADD, &
                  DATAKE, IGX, IGY, IGZ, VKPT, AVEC, CWORK1, CWORK2, EVALUE, CW, CH)
#ifdef CRAYWAR4
! in CCE20 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_VEC
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: AVEC(DIMREAL(GRID%MPLWV),3)     !< magnetic vector potential
     COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK2(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(3))
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
   ! local
     COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     REAL(q)    :: RINPLW
     INTEGER    :: ISPINOR,M,MM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED    CVR: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
#endif

     PROFILING_START('kinhamil_vec')

     RINPLW=2*PI/GRID%NPLWV

     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
      ! component along first axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGX) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGX) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         MM=M+ISPINOR*WDES1%NGVECTOR
         CWORK1(M)=CW(MM)*(IGX(M)+VKPT(1))
      ENDDO

      ! FFT to real space
      CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,GRID%RL%NP
         MM =M+ISPINOR *WDES1%GRID%MPLWV
         CVR(MM)= CVR(MM)+AVEC(M,1) *CWORK2(M) *RINPLW
      ENDDO

      ! component along second axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGY) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGY) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         MM=M+ISPINOR*WDES1%NGVECTOR
         CWORK1(M)=CW(MM)*(IGY(M)+VKPT(2))
      ENDDO

      ! FFT to real space
      CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,GRID%RL%NP
         MM =M+ISPINOR *WDES1%GRID%MPLWV
         CVR(MM)= CVR(MM)+AVEC(M,2) *CWORK2(M)*RINPLW
      ENDDO

      ! component along third axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGZ) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGZ) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         MM=M+ISPINOR*WDES1%NGVECTOR
         CWORK1(M)=CW(MM)*(IGZ(M)+VKPT(3))
      ENDDO

      ! FFT to real space
      CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,GRID%RL%NP
         MM =M+ISPINOR *WDES1%GRID%MPLWV
         CVR(MM)= CVR(MM)+AVEC(M,3) *CWORK2(M)*RINPLW
      ENDDO

      CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW,DATAKE) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CW,DATAKE) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         MM=M+ISPINOR*WDES1%NGVECTOR
         CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
      ENDDO
#else
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CWORK1,CW,IGX,IGY,IGZ,VKPT) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M,1)=CW(MM)*(IGX(M)+VKPT(1))
           CWORK1(M,2)=CW(MM)*(IGY(M)+VKPT(2))
           CWORK1(M,3)=CW(MM)*(IGZ(M)+VKPT(3))
        ENDDO
!$OMP END PARALLEL DO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,1),CWORK1(1,1),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,2),CWORK1(1,2),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,3),CWORK1(1,3),GRID)
!$OMP PARALLEL DO SHARED(GRID,WDES1,ISPINOR,AVEC,CWORK2,CVR,RINPLW) &
!$OMP PRIVATE(M,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CTMP1=AVEC(M,1) *CWORK2(M,1)
           CTMP2=AVEC(M,2) *CWORK2(M,2)
           CTMP3=AVEC(M,3) *CWORK2(M,3)
           CTMP4=CTMP1+CTMP2+CTMP3
           CVR(MM)= CVR(MM)+CTMP4*RINPLW
        ENDDO
!$OMP END PARALLEL DO
        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CH,EVALUE,CW,DATAKE) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)+CW(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_vec')

   END SUBROUTINE KINHAMIL_VEC


!************************* SUBROUTINE KINHAMIL_TAU *********************
!
! low level F77 to calculate
!> similar to #KINHAMIL
!> but adds the action of   - nabla mu  nabla
!>
!> - nabla mu nabla e(iGr)
!>
!> @details @ref openmp :
!> the workarrays CWORK1 and CWORK2 acquire a second dimension
!> (size=3 and 6, respectively).
!> The additional workspace is used to reorder the routine. First fill
!> all three components of CWORK1, FFT(CWORK1), build the product
!> with the field MU in CWORK2, FFT(CWORK2), and finally add the
!> relevant contributions to CH.
!> This improves vectorization of the loops over the elements of these
!> arrays. Aforementioned loops are distributed over all available
!> OpenMP threads.
!
!***********************************************************************

   SUBROUTINE KINHAMIL_TAU( WDES1, GRID, CVR, LADD, LKINETIC, DATAKE, &
                  IGX, IGY, IGZ, VKPT, LATT_CUR, MU, CWORK1, CWORK2, EVALUE, CW, CH)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_TAU
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: MU(DIMREAL(GRID%MPLWV)*WDES1%NRSPINORS*WDES1%NRSPINORS) ! \mu = d f_xc / d \tau
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     LOGICAL    :: LKINETIC
   ! local
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
     COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK2(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(6))
!$   COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
#ifdef memalign
!!DIR$ ASSUME_ALIGNED    CVR: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
#endif

     PROFILING_START('kinhamil_tau')

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)

      CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)

      IF (LKINETIC) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW,DATAKE) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CW,DATAKE) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
         ENDDO

      ELSE

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CW) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CH(MM)=CH(MM)-CW(MM)*EVALUE
         ENDDO

      ENDIF

      ! component along x-axis
GPACC !$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVR,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GX,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%NGVECTOR
            GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
            MM=M+ISPINOR_*WDES1%NGVECTOR
            CWORK1(M)=CW(MM)*GX*CITPI
         ENDDO

         ! FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,MU,CWORK2) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
            CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M) *RINPLW
         ENDDO

      ENDDO

      ! FFT to reciprocal space
      CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(GX,MM)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
         MM=M+ISPINOR*WDES1%NGVECTOR
         CH(MM)=CH(MM)-CWORK1(M)*GX*CITPI
      ENDDO

      ! component along y-axis
GPACC !$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVR,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GY,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%NGVECTOR
            GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
            MM=M+ISPINOR_*WDES1%NGVECTOR
            CWORK1(M)=CW(MM)*GY*CITPI
         ENDDO

         ! FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,MU,CWORK2) PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
            CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
         ENDDO

      ENDDO

      ! FFT to reciprocal space
      CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(GY,MM)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
         MM=M+ISPINOR*WDES1%NGVECTOR
         CH(MM)=CH(MM)-CWORK1(M)*GY*CITPI
      ENDDO

      ! component along z-axis
GPACC !$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVR,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GZ,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%NGVECTOR
            GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
            MM=M+ISPINOR_*WDES1%NGVECTOR
            CWORK1(M)=CW(MM)*GZ*CITPI
         ENDDO

         ! FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,MU,CWORK2) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
            CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
         ENDDO

      ENDDO

      ! FFT to reciprocal space
      CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(GZ,MM)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%NGVECTOR
         GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
         MM=M+ISPINOR*WDES1%NGVECTOR
         CH(MM)=CH(MM)-CWORK1(M)*GZ*CITPI
      ENDDO

#else
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
        IF (LKINETIC) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CH,CW,EVALUE,DATAKE) PRIVATE(M,MM)
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
           ENDDO
!$OMP END PARALLEL DO
        ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CH,CW,EVALUE) PRIVATE(M,MM)
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE
           ENDDO
!$OMP END PARALLEL DO
        ENDIF

        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,ISPINOR_,CWORK1,CW) PRIVATE(M,GX,GY,GZ,MM)
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M,1)=CW(MM)*GX*CITPI
              CWORK1(M,2)=CW(MM)*GY*CITPI
              CWORK1(M,3)=CW(MM)*GZ*CITPI
           ENDDO
!$OMP END PARALLEL DO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,1),CWORK1(1,1),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,2),CWORK1(1,2),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,3),CWORK1(1,3),GRID)

           IF (ISPINOR_==0) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK2) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK2(M,4)=MU(MM) *CWORK2(M,1)
                 CWORK2(M,5)=MU(MM) *CWORK2(M,2)
                 CWORK2(M,6)=MU(MM) *CWORK2(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK2) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK2(M,4)=CWORK2(M,4)+MU(MM) *CWORK2(M,1)
                 CWORK2(M,5)=CWORK2(M,5)+MU(MM) *CWORK2(M,2)
                 CWORK2(M,6)=CWORK2(M,6)+MU(MM) *CWORK2(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ENDIF
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,4),CWORK1(1,1),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,5),CWORK1(1,2),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,6),CWORK1(1,3),GRID,.FALSE.)

!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO &
!$OMP SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,VKPT,ISPINOR,CWORK1,CH,RINPLW,CW) &
!$OMP PRIVATE(M,GX,GY,GZ,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           CTMP1=CWORK1(M,1)*GX
           CTMP2=CWORK1(M,2)*GY
           CTMP3=CWORK1(M,3)*GZ
           CTMP4=CTMP1+CTMP2+CTMP3
           CTMP4=CTMP4*CITPI*RINPLW
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CTMP4
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_tau')

   END SUBROUTINE KINHAMIL_TAU


!************************* SUBROUTINE PSCURRENT ************************
!
!> gradient psi contribution to current: Im[psi^*(r) nabla psi(r)]
!> density contribution to current: psi^*(r) * psi(r)
!> for a supplied weighted k-point and psi
!> with psi(r) = e(iGr)
!> we get contributions proportional to (iG) e(iGr)
!
!***********************************************************************

   SUBROUTINE PSCURRENT( WDES1, GRID, CR, &
                  IGX, IGY, IGZ, VKPT, WEIGHT, AVEC, CWORK1, CWORK2, CW, J, C)
     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)  !< wavefunction in real space
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     REAL(q)    :: WEIGHT
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: AVEC(DIMREAL(GRID%MPLWV),3)     !< magnetic vector potential
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
     COMPLEX(q) :: CWORK2(GRID%MPLWV)
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     RGRID      :: J(GRID%RL%NP,3)                 !< three component vector storing current density in real space
     RGRID      :: C(GRID%RL%NP)                   !< charge density
   ! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
        ! component along first axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGX(M)+VKPT(1))*CITPI
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,1)=J(M,1)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

        ! component along second axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGY(M)+VKPT(2))*CITPI
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,2)=J(M,2)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

        ! component along third axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGZ(M)+VKPT(3))*CITPI
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,3)=J(M,3)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           C(M)=C(M)+CR(MM)*CONJG(CR(MM))*WEIGHT
        ENDDO
     ENDDO
   END SUBROUTINE PSCURRENT


!************************* SUBROUTINE KINHAMIL_C ***********************
!
!> low level routine for complex eigenvalues
!
!***********************************************************************

   SUBROUTINE KINHAMIL_C( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CW, CH)

#ifdef CRAYWAR4
! in CCE20 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_C
#endif
__GPU USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     COMPLEX(q) :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
   ! local
     INTEGER ISPINOR, M, MM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CVR: __align_BTS
#endif

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(CH,CW,DATAKE,WDES1) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
NOACC !$OMP PARALLEL DO DEFAULT(NONE) SHARED(WDES1,ISPINOR,CH,CW,EVALUE,DATAKE) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
!          CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
           CH(MM)=CH(MM)+CW(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
NOACC !$OMP END PARALLEL DO
     ENDDO
   END SUBROUTINE KINHAMIL_C


!************************* SUBROUTINE KINHAMIL_VEC_C *******************
!
!> low level routine for complex eigenvalues
!
!***********************************************************************

   SUBROUTINE KINHAMIL_VEC_C( WDES1, GRID, CVR,  LADD, &
                  DATAKE, IGX, IGY, IGZ, VKPT, AVEC, CWORK1, CWORK2, EVALUE, CW, CH)

#ifdef CRAYWAR4
! in CCE20 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_VEC_C
#endif
GPACC USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: AVEC(DIMREAL(GRID%MPLWV),3)     !< magnetic vector potential
     COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK2(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(3))
     COMPLEX(q) :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
   ! local
     COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     REAL(q)    :: RINPLW
     INTEGER    :: ISPINOR,M,MM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED    CVR: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
#endif

     PROFILING_START('kinhamil_vec_c')

     RINPLW=2*PI/GRID%NPLWV

     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
        ! component along first axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGX) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGX(M)+VKPT(1))
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,1) *CWORK2(M) *RINPLW
        ENDDO

        ! component along second axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGY) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGY(M)+VKPT(2))
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,2) *CWORK2(M)*RINPLW
        ENDDO

        ! component along third axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGZ) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGZ(M)+VKPT(3))
        ENDDO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,3) *CWORK2(M)*RINPLW
        ENDDO

        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW,DATAKE) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
        ENDDO
#else
!$OMP PARALLEL DO &
!$OMP SHARED(WDES1,ISPINOR,CWORK1,CW,IGX,IGY,IGZ,VKPT) &
!$OMP PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M,1)=CW(MM)*(IGX(M)+VKPT(1))
           CWORK1(M,2)=CW(MM)*(IGY(M)+VKPT(2))
           CWORK1(M,3)=CW(MM)*(IGZ(M)+VKPT(3))
        ENDDO
!$OMP END PARALLEL DO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,1),CWORK1(1,1),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,2),CWORK1(1,2),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,3),CWORK1(1,3),GRID)
!$OMP PARALLEL DO &
!$OMP SHARED(GRID,WDES1,ISPINOR,AVEC,CWORK2,CVR,RINPLW) &
!$OMP PRIVATE(M,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CTMP1=AVEC(M,1) *CWORK2(M,1)
           CTMP2=AVEC(M,2) *CWORK2(M,2)
           CTMP3=AVEC(M,3) *CWORK2(M,3)
           CTMP4=CTMP1+CTMP2+CTMP3
           CVR(MM)= CVR(MM)+CTMP4*RINPLW
        ENDDO
!$OMP END PARALLEL DO
        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO &
!$OMP SHARED(WDES1,ISPINOR,CH,EVALUE,CW,DATAKE) &
!$OMP PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)+CW(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_vec_c')

   END SUBROUTINE KINHAMIL_VEC_C


!************************* SUBROUTINE KINHAMIL_TAU_C *******************
!
!> low level routine for complex eigenvalues
!
!***********************************************************************

   SUBROUTINE KINHAMIL_TAU_C( WDES1, GRID, CVR, LADD, LKINETIC, DATAKE, &
                  IGX, IGY, IGZ, VKPT, LATT_CUR, MU, CWORK1, CWORK2, EVALUE, CW, CH)
#ifdef CRAYWAR4
! in CCE20 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_TAU_C
#endif

GPACC USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: MU(DIMREAL(GRID%MPLWV)*WDES1%NRSPINORS*WDES1%NRSPINORS) ! \mu = d f_xc / d \tau
     COMPLEX(q) :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     LOGICAL    :: LKINETIC
   ! local
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
     COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK2(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(6))
!$   COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
#ifdef memalign
!!DIR$ ASSUME_ALIGNED    CVR: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
#endif

     PROFILING_START('kinhamil_tau_c')

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
        IF (LKINETIC) THEN
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW,DATAKE) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
           ENDDO
        ELSE
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CW) PRIVATE(MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE
           ENDDO
        ENDIF

        ! component along x-axis
!$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GX,MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CW(MM)*GX*CITPI
           ENDDO
           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M) *RINPLW
           ENDDO
        ENDDO
        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
!$ACC PRIVATE(GX,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GX*CITPI
        ENDDO

        ! component along y-axis
!$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GY,MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CW(MM)*GY*CITPI
           ENDDO
           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO
        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
!$ACC PRIVATE(GY,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GY*CITPI
        ENDDO

        ! component along z-axis
!$ACC KERNELS PRESENT(CVR) __IF_ASYNC__
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CW) PRIVATE(GZ,MM) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CW(MM)*GZ*CITPI
           ENDDO
           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO
        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CW) __IF_ASYNC__ &
!$ACC PRIVATE(GZ,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GZ*CITPI
        ENDDO
#else
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
        IF (LKINETIC) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CH,CW,EVALUE,DATAKE) PRIVATE(M,MM)
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE+CW(MM)*DATAKE(M,ISPINOR+1)
           ENDDO
!$OMP END PARALLEL DO
        ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CH,CW,EVALUE) PRIVATE(M,MM)
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CW(MM)*EVALUE
           ENDDO
!$OMP END PARALLEL DO
        ENDIF

        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,ISPINOR_,CWORK1,CW) PRIVATE(M,GX,GY,GZ,MM)
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M,1)=CW(MM)*GX*CITPI
              CWORK1(M,2)=CW(MM)*GY*CITPI
              CWORK1(M,3)=CW(MM)*GZ*CITPI
           ENDDO
!$OMP END PARALLEL DO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,1),CWORK1(1,1),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,2),CWORK1(1,2),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,3),CWORK1(1,3),GRID)

           IF (ISPINOR_==0) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK2) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK2(M,4)=MU(MM) *CWORK2(M,1)
                 CWORK2(M,5)=MU(MM) *CWORK2(M,2)
                 CWORK2(M,6)=MU(MM) *CWORK2(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK2) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK2(M,4)=CWORK2(M,4)+MU(MM) *CWORK2(M,1)
                 CWORK2(M,5)=CWORK2(M,5)+MU(MM) *CWORK2(M,2)
                 CWORK2(M,6)=CWORK2(M,6)+MU(MM) *CWORK2(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ENDIF
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,4),CWORK1(1,1),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,5),CWORK1(1,2),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,6),CWORK1(1,3),GRID,.FALSE.)

!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO &
!$OMP SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,VKPT,ISPINOR,CWORK1,CH,RINPLW,CW) &
!$OMP PRIVATE(M,GX,GY,GZ,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           CTMP1=CWORK1(M,1)*GX
           CTMP2=CWORK1(M,2)*GY
           CTMP3=CWORK1(M,3)*GZ
           CTMP4=CTMP1+CTMP2+CTMP3
           CTMP4=CTMP4*CITPI*RINPLW
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CTMP4
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_tau_c')

   END SUBROUTINE KINHAMIL_TAU_C


!***********************************************************************
!
!> low level F77 version
!> determine norm of a wavefunction and norm with a metric
!> (if supplied)
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER PW_NORM_WITH_METRIC
#endif

__GPU USE moffload

    USE wave
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavedes1)    WDES1
    COMPLEX(q) :: C(WDES1%NPL)
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(WDES1%NPL)
  ! local
    INTEGER ISPINOR, M, MM

    IF (PRESENT (METRIC).AND. PRESENT(FMETRIC)) THEN

GPACC !$ACC DATA COPYOUT(FMETRIC,FNORM) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:FMETRIC,FNORM) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNORM, FMETRIC) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
       FNORM=0
       FMETRIC=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) REDUCTION(+:FNORM,FMETRIC) PRIVATE(MM) &
GPACC !$ACC PRESENT(C,METRIC,FMETRIC,FNORM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C,METRIC,FMETRIC,FNORM) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:FNORM,FMETRIC) IF(OFFLOAD_ON)

       DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(WDES1,ISPINOR,C,METRIC) REDUCTION(+:FNORM,FMETRIC)
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
             FMETRIC =FMETRIC+C(MM)*CONJG(C(MM))*METRIC(MM)
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

GPACC !$ACC END DATA
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:FMETRIC,FNORM) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.(.NOT.(IS_ON_DEVICE(FNORM,1).AND.IS_ON_DEVICE(FMETRIC,1))))
       CALLMPI( M_sum_2(WDES1%COMM_INB, FNORM, FMETRIC))

    ELSE

GPACC !$ACC DATA COPYOUT(FNORM) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:FNORM) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNORM) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
       FNORM=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) REDUCTION(+:FNORM) PRIVATE(MM) &
GPACC !$ACC PRESENT(C,FNORM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C,FNORM) &
GPOMP !$OMP PRIVATE(MM) REDUCTION(+:FNORM) IF(OFFLOAD_ON)

       DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(WDES1,ISPINOR,C) REDUCTION(+:FNORM)
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
          ENDDO
__CPU !$OMP END PARALLEL DO
       ENDDO

GPACC !$ACC END DATA
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:FNORM) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.(.NOT.IS_ON_DEVICE(FNORM,1)))
       CALLMPI( M_sum_d(WDES1%COMM_INB, FNORM, 1))
    ENDIF

  END SUBROUTINE PW_NORM_WITH_METRIC

!***********************************************************************
!
!> low level F77 routine
!> truncate high frequency components if LDELAY is set
!> or if use_eini_without_ldelay is defined
!> first version is scheduled for removal
!
!***********************************************************************


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY( WDES1, C, LDELAY, ENINI)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavedes1)    WDES1
    COMPLEX(q) C(WDES1%NPL)
    LOGICAL LDELAY     !< usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      !< cutoff at which the wavefunction is truncated
    ! local
    INTEGER ISPINOR

    PROFILING_START('truncate_high_frequency')

    IF (LDELAY.OR.WDES1%LSPIRAL) THEN
       DO ISPINOR=1,WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(WDES1%NGVECTOR, &
               C(1+(ISPINOR-1)*WDES1%NGVECTOR), WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF

    PROFILING_STOP('truncate_high_frequency')

  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE(NP, CW, DATAKE, ENINI)
#ifdef CRAYWAR1
! inlining this function with cray cce19 gives memory errors
!DIR$ INLINENEVER TRUNCATE_HIGH_FREQUENCY_ONE
#endif
__GPU USE moffload

    USE prec
    IMPLICIT NONE
    INTEGER NP
    COMPLEX(q) :: CW(NP)
    REAL(q)    :: DATAKE(NP)
    REAL(q)    :: ENINI
    ! local
    INTEGER N

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DATAKE,CW) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DATAKE,CW) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO PRIVATE(N) SHARED(NP,DATAKE,ENINI,CW)
    DO N=1,NP
       IF (DATAKE(N)>ENINI) CW(N)=0
    ENDDO
__CPU !$OMP END PARALLEL DO
  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE


!************************* SUBROUTINE ECCP_NL_ *************************
!
!> this subroutine calculates the expectation value of
!>    <c| D- evalue Q |cp>
!> where c and cp are two wavefunctions; non local part only
!> for one ion only
!
!***********************************************************************

  SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    COMPLEX(q)  CNL
    INTEGER LMMAXC, LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM),CQIJ(LMDIM,LMDIM)
    REAL(q) EVALUE
    GDEF CPROJ1(LMMAXC),CPROJ2(LMMAXC)
    ! local variables
    INTEGER L, LP

!$ACC LOOP COLLAPSE(2) REDUCTION(+:CNL) VECTOR
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+(CDIJ(LP,L)-EVALUE*CQIJ(LP,L))*CPROJ1(LP)*GCONJG(CPROJ2(L))
       ENDDO
    ENDDO
  END SUBROUTINE ECCP_NL_

#ifndef oldhamil
   SUBROUTINE VMUJVEC(WDES1,GRID,SV,CR,CVR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER VHAMIL
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1

     RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
     ! local variables
     REAL(q) RINPLW
     INTEGER ISPINOR, ISPINOR_, M, MM, MM_

     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CVR: __align_BTS
#endif

     PROFILING_START('vmujvec')

     RINPLW=1._q/GRID%NPLWV

     IF (WDES1%NRSPINORS==1) THEN
        PROFILING_STOP('vmujvec')
        CALL vtutor%bug('VMUJVEC: only implemented for spinor wavefunctions', __FILE__, __LINE__)
        RETURN
     ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRID,CVR,SV,CR) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL SHARED(GRID,CVR,SV,CR,RINPLW) PRIVATE(ISPINOR,M,MM)
        DO ISPINOR =0,1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)
            DO M=1,GRID%RL%NP
               MM =M+ISPINOR *GRID%MPLWV
               CVR(MM)=CVR(MM)-CI*SV(M,1+2*ISPINOR)*CR(M)*RINPLW
            ENDDO
__CPU !$OMP END DO
__GPU   ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRID,CVR,SV,CR) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

__GPU   DO ISPINOR =0,1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
__CPU !$OMP DO SCHEDULE(STATIC)
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              CVR(MM)= CVR(MM)-CI*SV(M,2+2*ISPINOR)*CR(M+GRID%MPLWV)*RINPLW
           ENDDO
__CPU !$OMP END DO
        ENDDO
__CPU !$OMP END PARALLEL
     ENDIF

     PROFILING_STOP('vmujvec')

   END SUBROUTINE VMUJVEC


   SUBROUTINE KINHAMIL_VEC_( WDES1, GRID, IGX, IGY, IGZ, VKPT, AVEC, CW, CWORK1, CWORK2, CVR)

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: AVEC(DIMREAL(GRID%MPLWV),3)     !< magnetic vector potential
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CWORK1(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK2(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(3))
   ! local
     COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     REAL(q)    :: RINPLW
     INTEGER    :: ISPINOR,M,MM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED    CVR: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
#endif

     PROFILING_START('kinhamil_vec_')

     RINPLW=2*PI/GRID%NPLWV

     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
        ! component along first axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGX) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGX) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGX(M)+VKPT(1))
        ENDDO

        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,1) *CWORK2(M) *RINPLW
        ENDDO

        ! component along second axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGY) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGY) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGY(M)+VKPT(2))
        ENDDO

        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,2) *CWORK2(M)*RINPLW
        ENDDO

        ! component along third axis
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CW,IGZ) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CWORK1,CW,IGZ) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CW(MM)*(IGZ(M)+VKPT(3))
        ENDDO

        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CVR,AVEC,CWORK2) &
GPOMP !$OMP PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,3) *CWORK2(M)*RINPLW
        ENDDO
#else
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CWORK1,CW,IGX,IGY,IGZ,VKPT) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M,1)=CW(MM)*(IGX(M)+VKPT(1))
           CWORK1(M,2)=CW(MM)*(IGY(M)+VKPT(2))
           CWORK1(M,3)=CW(MM)*(IGZ(M)+VKPT(3))
        ENDDO
!$OMP END PARALLEL DO
        ! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,1),CWORK1(1,1),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,2),CWORK1(1,2),GRID)
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1,3),CWORK1(1,3),GRID)
!$OMP PARALLEL DO SHARED(GRID,WDES1,ISPINOR,AVEC,CWORK2,CVR,RINPLW) &
!$OMP PRIVATE(M,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CTMP1=AVEC(M,1) *CWORK2(M,1)
           CTMP2=AVEC(M,2) *CWORK2(M,2)
           CTMP3=AVEC(M,3) *CWORK2(M,3)
           CTMP4=CTMP1+CTMP2+CTMP3
           CVR(MM)= CVR(MM)+CTMP4*RINPLW
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_vec_')

   END SUBROUTINE KINHAMIL_VEC_


   SUBROUTINE KINHAMIL_TAU_( WDES1, GRID, LATT_CUR, IGX, IGY, IGZ, VKPT, MU, CWORK1, CWORK2, CWORK3, CW, CH)

#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_TAU_
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: MU(DIMREAL(GRID%MPLWV)*WDES1%NRSPINORS*WDES1%NRSPINORS) ! \mu = d f_xc / d \tau
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
   ! work arrays
     COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK3(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(6))
   ! local
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
!$   COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK3: __align_BTS
#endif

     PROFILING_START('kinhamil_tau_')

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
        ! component along x-axis
GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GX,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GX*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MU,CWORK3) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MU,CWORK3) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MU(MM_) *CWORK3(M) *RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GX,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GX*CITPI
        ENDDO

        ! component along y-axis
GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GY,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GY*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MU,CWORK3) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MU,CWORK3) PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MU(MM_) *CWORK3(M)*RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GY,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GY*CITPI
        ENDDO

        ! component along z-axis
GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GZ,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GZ*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MU,CWORK3) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MU,CWORK3) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MU(MM_) *CWORK3(M)*RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GZ,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GZ*CITPI
        ENDDO
#else
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,ISPINOR_,CWORK2,CW) PRIVATE(M,GX,GY,GZ,MM)
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M,1)=CW(MM)*GX*CITPI
              CWORK2(M,2)=CW(MM)*GY*CITPI
              CWORK2(M,3)=CW(MM)*GZ*CITPI
           ENDDO
!$OMP END PARALLEL DO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,1),CWORK2(1,1),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,2),CWORK2(1,2),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,3),CWORK2(1,3),GRID)

           IF (ISPINOR_==0) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK3) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK3(M,4)=MU(MM) *CWORK3(M,1)
                 CWORK3(M,5)=MU(MM) *CWORK3(M,2)
                 CWORK3(M,6)=MU(MM) *CWORK3(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MU,CWORK3) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK3(M,4)=CWORK3(M,4)+MU(MM) *CWORK3(M,1)
                 CWORK3(M,5)=CWORK3(M,5)+MU(MM) *CWORK3(M,2)
                 CWORK3(M,6)=CWORK3(M,6)+MU(MM) *CWORK3(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ENDIF
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,4),CWORK2(1,1),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,5),CWORK2(1,2),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,6),CWORK2(1,3),GRID,.FALSE.)

!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,VKPT,ISPINOR,CWORK2,CH,RINPLW,CW) &
!$OMP PRIVATE(M,GX,GY,GZ,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           CTMP1=CWORK2(M,1)*GX
           CTMP2=CWORK2(M,2)*GY
           CTMP3=CWORK2(M,3)*GZ
           CTMP4=CTMP1+CTMP2+CTMP3
           CTMP4=CTMP4*CITPI*RINPLW
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CTMP4
        ENDDO
!$OMP END PARALLEL DO
#endif
     ENDDO

     PROFILING_STOP('kinhamil_tau_')

   END SUBROUTINE KINHAMIL_TAU_


   SUBROUTINE KINHAMIL_JVEC_( WDES1, GRID, LATT_CUR, IGX, IGY, IGZ, VKPT, MUTAU, MUJVEC, CWORK1, CWORK2, CWORK3, CWORK4, CW, CH)

#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER KINHAMIL_JVEC_
#endif

__GPU USE moffload

     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     RGRID      :: MUTAU (DIMREAL(GRID%MPLWV)*WDES1%NRSPINORS*WDES1%NRSPINORS)   !< \mu_tau  = d f_xc / d \tau
     RGRID      :: MUJVEC(DIMREAL(GRID%MPLWV)*WDES1%NRSPINORS*WDES1%NRSPINORS,3) !< \mu_jvec = d f_xc / d \jvec
     COMPLEX(q) :: CW(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
   ! work arrays
     COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV) __NOOFF_omp_arg(3))
     COMPLEX(q) :: CWORK3(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(6))
#if defined(OFFLOADING) || ! defined(_OPENMP)
     COMPLEX(q) :: CWORK4(GRID%MPLWV*WDES1%NRSPINORS)
#else
     COMPLEX(q) :: CWORK4(__align_DIM(GRID%MPLWV)   __NOOFF_omp_arg(3))
#endif
   ! local
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
!$   COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CWORK1: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK2: __align_BTS
!!DIR$ ASSUME_ALIGNED CWORK3: __align_BTS
#endif

     PROFILING_START('kinhamil_jvec_')

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1
#if defined(OFFLOADING) || ! defined(_OPENMP)
        ! component along x-axis
GPACC !$ACC KERNELS PRESENT(CWORK1,CWORK4) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,CWORK4,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        CWORK4(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GX,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GX*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MUTAU (MM_  ) *CWORK3(M) *RINPLW
              CWORK4(MM)=CWORK4(MM)+MUJVEC(MM_,1) *CWORK3(M) *RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GX,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GX,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GX*CITPI
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CWORK2) PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*CI*2._q
        ENDDO


        ! component along y-axis
GPACC !$ACC KERNELS PRESENT(CWORK1,CWORK4) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,CWORK4,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        CWORK4(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GY,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GY*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MUTAU (MM_  ) *CWORK3(M)*RINPLW
              CWORK4(MM)=CWORK4(MM)+MUJVEC(MM_,2) *CWORK3(M)*RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GY,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GY,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GY*CITPI
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CWORK2) PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*CI*2._q
        ENDDO


        ! component along z-axis
GPACC !$ACC KERNELS PRESENT(CWORK1,CWORK4) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1,CWORK4,WDES1%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CWORK1(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
        CWORK4(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        DO ISPINOR_=0,WDES1%NRSPINORS-1

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) PRIVATE(GZ,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK2,CW) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M)=CW(MM)*GZ*CITPI
           ENDDO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1),CWORK2(1),GRID)

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) PRIVATE(MM,MM_) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,WDES1,WDES1%GRID,CWORK1,MUTAU,CWORK3,MUJVEC,CWORK4) &
GPOMP !$OMP PRIVATE(MM,MM_) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
              CWORK1(MM)=CWORK1(MM)+MUTAU (MM_  ) *CWORK3(M)*RINPLW
              CWORK4(MM)=CWORK4(MM)+MUJVEC(MM_,3) *CWORK3(M)*RINPLW
           ENDDO
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) PRIVATE(GZ,MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK2) &
GPOMP !$OMP PRIVATE(GZ,MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*GZ*CITPI
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1+ISPINOR*WDES1%GRID%MPLWV),CWORK2(1),GRID,.FALSE.)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CH,CWORK2) PRIVATE(MM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,CH,CWORK2) PRIVATE(MM) IF(OFFLOAD_ON)

!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK2(M)*CI*2._q
        ENDDO
#else // defined(OFFLOADING) || ! defined(_OPENMP)
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,ISPINOR_,CWORK2,CW) PRIVATE(M,GX,GY,GZ,MM)
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK2(M,1)=CW(MM)*GX*CITPI
              CWORK2(M,2)=CW(MM)*GY*CITPI
              CWORK2(M,3)=CW(MM)*GZ*CITPI
           ENDDO
!$OMP END PARALLEL DO

           ! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,1),CWORK2(1,1),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,2),CWORK2(1,2),GRID)
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,3),CWORK2(1,3),GRID)

           IF (ISPINOR_==0) THEN
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MUTAU,CWORK3,MUJVEC,CWORK4) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK3(M,4)=MUTAU(MM) *CWORK3(M,1)
                 CWORK3(M,5)=MUTAU(MM) *CWORK3(M,2)
                 CWORK3(M,6)=MUTAU(MM) *CWORK3(M,3)

                 CWORK4(M,1)=MUJVEC(MM,1) *CWORK3(M,1)
                 CWORK4(M,2)=MUJVEC(MM,2) *CWORK3(M,2)
                 CWORK4(M,3)=MUJVEC(MM,3) *CWORK3(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ELSE
!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(GRID,ISPINOR,ISPINOR_,MUTAU,CWORK3) PRIVATE(M,MM)
              DO M=1,GRID%RL%NP
                 MM=M+(ISPINOR_+2*ISPINOR)*DIMREAL(GRID%MPLWV)
                 CWORK3(M,4)=CWORK3(M,4)+MUTAU(MM) *CWORK3(M,1)
                 CWORK3(M,5)=CWORK3(M,5)+MUTAU(MM) *CWORK3(M,2)
                 CWORK3(M,6)=CWORK3(M,6)+MUTAU(MM) *CWORK3(M,3)

                 CWORK4(M,1)=CWORK4(M,1)+MUJVEC(MM,1) *CWORK3(M,1)
                 CWORK4(M,2)=CWORK4(M,2)+MUJVEC(MM,2) *CWORK3(M,2)
                 CWORK4(M,3)=CWORK4(M,3)+MUJVEC(MM,3) *CWORK3(M,3)
              ENDDO
!$OMP END PARALLEL DO
           ENDIF
        ENDDO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,4),CWORK2(1,1),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,5),CWORK2(1,2),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK3(1,6),CWORK2(1,3),GRID,.FALSE.)

!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,IGX,IGY,IGZ,LATT_CUR,VKPT,ISPINOR,CWORK2,CH,RINPLW) &
!$OMP PRIVATE(M,GX,GY,GZ,MM,CTMP1,CTMP2,CTMP3,CTMP4)
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           CTMP1=CWORK2(M,1)*GX
           CTMP2=CWORK2(M,2)*GY
           CTMP3=CWORK2(M,3)*GZ
           CTMP4=CTMP1+CTMP2+CTMP3
           CTMP4=CTMP4*CITPI*RINPLW
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CTMP4
        ENDDO
!$OMP END PARALLEL DO

        ! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1,1),CWORK2(1,1),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1,2),CWORK2(1,2),GRID,.FALSE.)
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK4(1,3),CWORK2(1,3),GRID,.FALSE.)

!DIR$ IVDEP
!OCL NOVREL
!$OMP PARALLEL DO SHARED(WDES1,ISPINOR,CWORK2,CH,RINPLW) PRIVATE(MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-(CWORK2(M,1)+CWORK2(M,2)+CWORK2(M,3))*RINPLW*CI*2._q
        ENDDO
!$OMP END PARALLEL DO
#endif // defined(OFFLOADING) || ! defined(_OPENMP)
     ENDDO

     PROFILING_STOP('kinhamil_jvec_')

   END SUBROUTINE KINHAMIL_JVEC_
#endif // oldhamil
