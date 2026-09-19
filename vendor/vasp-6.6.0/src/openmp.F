#include "symbol.inc"
    MODULE mopenmp

      USE mopenmp_struct_def

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_OMP

#ifdef OMP_OFFLOAD
      PUBLIC :: INIT_OFFLOAD, SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX
#ifndef IFX_BEFORE_20241008
      PUBLIC :: OMPX_DUMP_MAPPING_TABLES
#endif

      INTERFACE
         SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(SPLIT_CMPLX,TO_ADD)
            USE prec
            COMPLEX(q) :: SPLIT_CMPLX, TO_ADD
         END SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX

#ifndef IFX_BEFORE_20241008
         SUBROUTINE OMPX_DUMP_MAPPING_TABLES() BIND(C)
         END SUBROUTINE OMPX_DUMP_MAPPING_TABLES
#endif

      END INTERFACE
#endif //OMP_OFFLOAD

    CONTAINS

!****************** SUBROUTINE INIT_OMP ********************************
!
!> @details @ref openmp :
!> this subroutine reads
!>   mopenmp::omp_nthreads
!>   mopenmp::omp_nthreads_alltoall
!>   mopenmp::omp_nthreads_nonlr_rspace
!>   mopenmp::omp_nonlr_planewise
!>
!> When none of these tags are specified in the INCAR file,
!>   mopenmp::omp_nthreads = omp_get_max_threads() (OMP_NUM_THREADS),
!>   mopenmp::omp_nthreads_nonlr_rspace = mopenmp::omp_nthreads
!
!***********************************************************************

      SUBROUTINE INIT_OMP(IO)
      USE prec
      USE base
      USE reader_tags
      TYPE (in_struct) IO
      ! local
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LDUM
      CHARACTER(1) CHARAC
#ifdef _OPENMP
      INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

      ! prevent nesting in openMP
      CALL OMP_SET_MAX_ACTIVE_LEVELS(1)

      omp_nthreads=OMP_GET_MAX_THREADS()

      IF (omp_nthreads<2) THEN
         omp_nthreads=1 ; RETURN
      ENDIF

      omp_nthreads_alltoall=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NTHREADS_A2A',omp_nthreads_alltoall,IERR)

      IF (omp_nthreads_alltoall<1)            omp_nthreads_alltoall=1
      IF (omp_nthreads_alltoall>omp_nthreads) omp_nthreads_alltoall=omp_nthreads

      omp_nthreads_nonlr_rspace=omp_nthreads
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NTHREADS_RPROJ',omp_nthreads_nonlr_rspace,IERR)

      IF (omp_nthreads_nonlr_rspace<1.OR.omp_nthreads_nonlr_rspace>omp_nthreads)  &
     &   omp_nthreads_nonlr_rspace=omp_nthreads
#ifdef OFFLOADING
      omp_nthreads_nonlr_rspace=1
#endif
#endif
      RETURN
      END SUBROUTINE INIT_OMP

#ifdef OMP_OFFLOAD
!****************** SUBROUTINE INIT_OFFLOAD ****************************
!
!***********************************************************************

      SUBROUTINE INIT_OFFLOAD(IO)
      USE base
      USE main_mpi

      USE moffload_struct_def, ONLY : OFFLOAD_ALLOW, OFFLOAD_ON, &
                                      OFFLOAD_NUM_DEVICES, OFFLOAD_MAX_RANKS_PER_DEVICE

      USE mdevice_info, ONLY : GET_NUM_DEVICES_GLOBAL, GET_NUM_RANKS_ON_DEVICE

      USE tutor, ONLY: vtutor

      TYPE(in_struct) IO
      ! local variables
      INTEGER :: NUM_DEVICES, DEVICE_NUM
      LOGICAL :: LALL_RANKS_HAVE_DEVICE

      CHARACTER(LEN=256) :: ENVVAR_VALUE
      INTEGER :: ENVVAR_STAT

      INTEGER, EXTERNAL :: OMP_GET_NUM_DEVICES

      PROFILING_START('init_offload')

      NUM_DEVICES = OMP_GET_NUM_DEVICES()

      IF (.NOT.ALL_RANKS_HAVE_DEVICE(NUM_DEVICES)) THEN
         OFFLOAD_ALLOW = .FALSE.
         PROFILING_STOP('init_offload')
         RETURN
      ENDIF

      CALL GET_ENVIRONMENT_VARIABLE("OMP_TARGET_OFFLOAD", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
      IF (ENVVAR_STAT==0 .AND. (ENVVAR_VALUE=='DISABLED' .OR. ENVVAR_VALUE=='disabled')) THEN
         OFFLOAD_ALLOW = .FALSE.
#ifdef USENCCL
         LUSENCCL = .FALSE.
#endif
         PROFILING_STOP('init_offload')
         RETURN
      ENDIF

      OFFLOAD_ALLOW = .TRUE.


#ifdef MPI
      DEVICE_NUM=MOD(COMM_WORLD%NODE_ME-1, NUM_DEVICES)
#endif
!!      GPU_NUM = OMP_GET_DEVICE_NUM()
      CALL OMP_SET_DEFAULT_DEVICE(DEVICE_NUM)

      OFFLOAD_NUM_DEVICES = GET_NUM_DEVICES_GLOBAL(DEVICE_NUM)
      OFFLOAD_MAX_RANKS_PER_DEVICE = MAX_RANKS_PER_DEVICE(DEVICE_NUM)

      ! Store DEVICE_NUM (might be needed elsewhere lateron)
      OMP_MY_DEVICE_NUM = DEVICE_NUM

      IF (IO%IU0>=0) WRITE(IO%IU0,'(" Offloading initialized ... ",I4," GPUs detected")') OFFLOAD_NUM_DEVICES

#if defined(USENCCL) && defined(MPI)
    IF (LUSENCCL) THEN
       CALL nccl_init_comm(COMM_INTER)
       CALL nccl_init_comm(COMM_INB)
       CALL nccl_init_comm(COMM_KIN)
       CALL nccl_init_comm(COMM_KINTER)

       ! were all NCCL communicators created succesfully?
       LUSENCCL = COMM_INTER%LUSENCCL .AND. COMM_INB%LUSENCCL .AND. COMM_KIN%LUSENCCL .AND. COMM_KINTER%LUSENCCL
       CALL M_and(COMM_WORLD, LUSENCCL, 1)

       IF (IO%IU0>=0 .AND. LUSENCCL) WRITE(IO%IU0,'(" RCCL MPI communication initialized ...")')

       ! if not write a warning that NCCL has been switched off (completely for now) ...
       IF (IO%IU0>=0 .AND. .NOT. LUSENCCL) &
          CALL vtutor%warning('INIT_OFFLOAD: several MPI-ranks need to share a GPU, which is not supported ' // &
                              'by NCCL. The use of NCCL will be switched off. To avoid this, reduce '// &
                              'the number of MPI-ranks: #-of-ranks <= #-of-GPUs (on every node!).')

    ENDIF
#endif

      PROFILING_STOP('init_offload')

      CONTAINS

         LOGICAL FUNCTION ALL_RANKS_HAVE_DEVICE(NUM_DEVICES)
         INTEGER :: NUM_DEVICES
         ALL_RANKS_HAVE_DEVICE = NUM_DEVICES > 0
         CALLMPI( M_and(COMM_WORLD, ALL_RANKS_HAVE_DEVICE, 1) )
         END FUNCTION ALL_RANKS_HAVE_DEVICE


         INTEGER FUNCTION MAX_RANKS_PER_DEVICE(DEVICE_NUM)
         INTEGER :: DEVICE_NUM
         MAX_RANKS_PER_DEVICE = GET_NUM_RANKS_ON_DEVICE(DEVICE_NUM)
         CALLMPI( M_max_i(COMM_inter_node_world, MAX_RANKS_PER_DEVICE, 1) )
         END FUNCTION MAX_RANKS_PER_DEVICE

      END SUBROUTINE INIT_OFFLOAD

#endif //OMP_OFFLOAD
    END MODULE mopenmp
#ifdef OMP_OFFLOAD

!****************** SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX *******
!
! This routine takes one complex argument and accepts this as two real
! arguments for the real and imaginary parts, respectively. It then
! performs two separate OpenMP atomic adds, one for each part.
!
!***********************************************************************
      SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(SPLIT_CMPLX,TO_ADD)
      USE prec
!$OMP DECLARE TARGET (SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX)
      REAL(q),DIMENSION(2) :: SPLIT_CMPLX
      COMPLEX(q) :: TO_ADD
!$OMP ATOMIC UPDATE
      SPLIT_CMPLX(1)=SPLIT_CMPLX(1)+REAL(TO_ADD,q) ! real part
!$OMP END ATOMIC
!$OMP ATOMIC UPDATE
      SPLIT_CMPLX(2)=SPLIT_CMPLX(2)+AIMAG(TO_ADD) ! imaginary part
!$OMP END ATOMIC
      END SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX


!***********************************************************************
!
!***********************************************************************

    MODULE mopenmp_deepcopy

      USE moffload_struct_def

      USE mgrid_struct_def
      USE wave_struct_def
      USE nonlr_struct_def
      USE nonl_struct_def
      USE pseudo_struct_def
      USE mkpoints_struct_def
      USE setexm_struct_def
      USE chi_glb
      USE pot_struct_def

      USE tutor, ONLY : vtutor

      USE iso_c_binding, ONLY : c_loc

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: COPYIN_TYPED_VAR,CREATE_TYPED_VAR,DELETE_TYPED_VAR,COPYOUT_TYPED_VAR

      INTERFACE COPYIN_TYPED_VAR
         MODULE PROCEDURE COPYIN_NONL_STRUCT
         MODULE PROCEDURE COPYIN_NONLR_STRUCT
         MODULE PROCEDURE COPYIN_NONLR_PROJ
         MODULE PROCEDURE COPYIN_LAYOUT
         MODULE PROCEDURE COPYIN_GRID_MAP
         MODULE PROCEDURE COPYIN_GRID_3D
         MODULE PROCEDURE COPYIN_TRANSIT
         MODULE PROCEDURE COPYIN_KERNEL_TRUNCATION
         MODULE PROCEDURE COPYIN_COULOMB_POTENTIAL
         MODULE PROCEDURE COPYIN_WAVEDES
         MODULE PROCEDURE COPYIN_WAVEDES1
         MODULE PROCEDURE COPYIN_WAVESPIN
         MODULE PROCEDURE COPYIN_WAVEFUNA
         MODULE PROCEDURE COPYIN_WAVEFUN1
         MODULE PROCEDURE COPYIN_WAVEFUN1_ARRAY
         MODULE PROCEDURE COPYIN_POTCARS
         MODULE PROCEDURE COPYIN_KPOINTS_FULL
         MODULE PROCEDURE COPYIN_GREENSFDES
         MODULE PROCEDURE COPYIN_XC
      END INTERFACE COPYIN_TYPED_VAR

      INTERFACE CREATE_TYPED_VAR
         MODULE PROCEDURE CREATE_WAVESPIN
      END INTERFACE CREATE_TYPED_VAR

      INTERFACE DELETE_TYPED_VAR
         MODULE PROCEDURE DELETE_NONL_STRUCT
         MODULE PROCEDURE DELETE_NONLR_STRUCT
         MODULE PROCEDURE DELETE_NONLR_PROJ
         MODULE PROCEDURE DELETE_LAYOUT
         MODULE PROCEDURE DELETE_GRID_MAP
         MODULE PROCEDURE DELETE_GRID_3D
         MODULE PROCEDURE DELETE_TRANSIT
         MODULE PROCEDURE DELETE_KERNEL_TRUNCATION
         MODULE PROCEDURE DELETE_COULOMB_POTENTIAL
         MODULE PROCEDURE DELETE_WAVEDES
         MODULE PROCEDURE DELETE_WAVEDES1
         MODULE PROCEDURE DELETE_WAVESPIN
         MODULE PROCEDURE DELETE_WAVEFUNA
         MODULE PROCEDURE DELETE_WAVEFUN1
         MODULE PROCEDURE DELETE_WAVEFUN1_ARRAY
         MODULE PROCEDURE DELETE_POTCARS
         MODULE PROCEDURE DELETE_KPOINTS_FULL
         MODULE PROCEDURE DELETE_GREENSFDES
         MODULE PROCEDURE DELETE_XC
      END INTERFACE DELETE_TYPED_VAR

      INTERFACE COPYOUT_TYPED_VAR
         MODULE PROCEDURE COPYOUT_WAVESPIN
      END INTERFACE COPYOUT_TYPED_VAR

    CONTAINS

!****************** SUBROUTINE COPYIN_NONL_STRUCT **********************
!
!***********************************************************************

      SUBROUTINE COPYIN_NONL_STRUCT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_NONL_STRUCT
#endif
      TYPE(nonl_struct) :: VAR
      ! early exit if possible
      IF (.NOT.OFFLOAD_ON) RETURN

      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%NITYP, &
!$OMP                           VAR%ITYP, &
!$OMP                           VAR%LMMAX, &
!$OMP                           VAR%LMBASE, &
!$OMP                           VAR%POSION, &
!$OMP                           VAR%VKPT_SHIFT, &
!$OMP                           VAR%CQFAK, &
!$OMP                           VAR%QPROJ, &
!$OMP                           VAR%CREXP, &
!$OMP                           VAR%CPROJK)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NITYP, &
!$OMP                                      VAR%ITYP, &
!$OMP                                      VAR%LMMAX, &
!$OMP                                      VAR%LMBASE, &
!$OMP                                      VAR%POSION, &
!$OMP                                      VAR%VKPT_SHIFT, &
!$OMP                                      VAR%CQFAK, &
!$OMP                                      VAR%QPROJ, &
!$OMP                                      VAR%CREXP, &
!$OMP                                      VAR%CPROJK)
#endif
      END SUBROUTINE COPYIN_NONL_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_STRUCT *********************
!
!***********************************************************************

      SUBROUTINE COPYIN_NONLR_STRUCT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_NONLR_STRUCT
#endif
      TYPE(nonlr_struct) :: VAR
      ! local variables
      INTEGER :: I
      ! early exit if possible
      IF (.NOT.OFFLOAD_ON) RETURN

      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%NITYP, &
!$OMP                           VAR%ITYP, &
!$OMP                           VAR%LMAX, &
!$OMP                           VAR%LMMAX, &
!$OMP                           VAR%CHANNELS, &
!$OMP                           VAR%PSRMAX, &
!$OMP                           VAR%RSMOOTH, &
!$OMP                           VAR%POSION, &
!$OMP                           VAR%VKPT_SHIFT, &
!$OMP                           VAR%NLIMAX, &
!$OMP                           VAR%NLI, &
!$OMP                           VAR%LMBASE, &
!$OMP                           VAR%NLIBASE, &
!$OMP                           VAR%RPROJ, &
!$OMP                           VAR%CRREXP)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NITYP, &
!$OMP                                      VAR%ITYP, &
!$OMP                                      VAR%LMAX, &
!$OMP                                      VAR%LMMAX, &
!$OMP                                      VAR%CHANNELS, &
!$OMP                                      VAR%PSRMAX, &
!$OMP                                      VAR%RSMOOTH, &
!$OMP                                      VAR%POSION, &
!$OMP                                      VAR%VKPT_SHIFT, &
!$OMP                                      VAR%NLIMAX, &
!$OMP                                      VAR%NLI, &
!$OMP                                      VAR%LMBASE, &
!$OMP                                      VAR%NLIBASE, &
!$OMP                                      VAR%RPROJ, &
!$OMP                                      VAR%CRREXP)
#endif
      ! derived type pointers
      DO I=1,SIZE(VAR%BETA)
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%BETA(I))
         CALL COPYIN_TYPED_VAR(VAR%BETA(I))
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%BETA(I))
      ENDDO
      END SUBROUTINE COPYIN_NONLR_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_PROJ ***********************
!
!***********************************************************************

      SUBROUTINE COPYIN_NONLR_PROJ(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_NONLR_PROJ
#endif
      TYPE(nonlr_proj) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%PSPRNL, &
!$OMP                           VAR%LPS)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%PSPRNL, &
!$OMP                                      VAR%LPS)
#endif
      END SUBROUTINE COPYIN_NONLR_PROJ


!****************** SUBROUTINE COPYIN_LAYOUT ***************************
!
!***********************************************************************

      SUBROUTINE COPYIN_LAYOUT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_LAYOUT
#endif
      TYPE(layout) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%I2, &
!$OMP                           VAR%I3, &
!$OMP                           VAR%INDEX)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%I2, &
!$OMP                                      VAR%I3, &
!$OMP                                      VAR%INDEX)
#endif
      END SUBROUTINE COPYIN_LAYOUT


!****************** SUBROUTINE COPYIN_GRID_MAP *************************
!
!***********************************************************************

      SUBROUTINE COPYIN_GRID_MAP(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_GRID_MAP
#endif
      TYPE(grid_map) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%N, &
!$OMP                           VAR%PTR, &
!$OMP                           VAR%TBL, &
!$OMP                           VAR%NI, &
!$OMP                           VAR%PTRI, &
!$OMP                           VAR%TBLI)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%N, &
!$OMP                                      VAR%PTR, &
!$OMP                                      VAR%TBL, &
!$OMP                                      VAR%NI, &
!$OMP                                      VAR%PTRI, &
!$OMP                                      VAR%TBLI)
#endif
      END SUBROUTINE COPYIN_GRID_MAP


!****************** SUBROUTINE COPYIN_GRID_3D **************************
!
!***********************************************************************

      SUBROUTINE COPYIN_GRID_3D(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_GRID_3D
#endif
      TYPE(grid_3d) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%LPCTX, &
!$OMP                           VAR%LPCTY, &
!$OMP                           VAR%LPCTZ, &
!$OMP                           VAR%LPCTX_, &
!$OMP                           VAR%LPCTY_, &
!$OMP                           VAR%LPCTZ_, &
!$OMP                           VAR%FFTSCA, &
!$OMP                           VAR%FFTWEIGHT, &
!$OMP                           VAR%NINDPWCONJG, &
!$OMP                           VAR%IND_IN_SPHERE, &
!$OMP                           VAR%COMM, &
!$OMP                           VAR%COMM_KIN, &
!$OMP                           VAR%COMM_KINTER)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%LPCTX, &
!$OMP                                      VAR%LPCTY, &
!$OMP                                      VAR%LPCTZ, &
!$OMP                                      VAR%LPCTX_, &
!$OMP                                      VAR%LPCTY_, &
!$OMP                                      VAR%LPCTZ_, &
!$OMP                                      VAR%FFTSCA, &
!$OMP                                      VAR%FFTWEIGHT, &
!$OMP                                      VAR%NINDPWCONJG, &
!$OMP                                      VAR%IND_IN_SPHERE, &
!$OMP                                      VAR%COMM, &
!$OMP                                      VAR%COMM_KIN, &
!$OMP                                      VAR%COMM_KINTER)
#endif
      ! derived type pointers
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%RC, VAR%IN, VAR%RL, VAR%RC_IN, VAR%IN_RL)
      CALL COPYIN_TYPED_VAR(VAR%RC)
      CALL COPYIN_TYPED_VAR(VAR%IN)
      CALL COPYIN_TYPED_VAR(VAR%RL)
      CALL COPYIN_TYPED_VAR(VAR%RC_IN)
      CALL COPYIN_TYPED_VAR(VAR%IN_RL)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%RC, VAR%IN, VAR%RL, VAR%RC_IN, VAR%IN_RL)
      ! if associated
      IF (ASSOCIATED(VAR%RL_FFT)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%RL_FFT)
         CALL COPYIN_TYPED_VAR(VAR%RL_FFT)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%RL_FFT)
      ENDIF
      END SUBROUTINE COPYIN_GRID_3D


!****************** SUBROUTINE COPYIN_TRANSIT **************************
!
!***********************************************************************

      SUBROUTINE COPYIN_TRANSIT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_TRANSIT
#endif
      TYPE(transit) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR)
!$OMP TARGET ENTER DATA MAP(TO: VAR%IND1, &
!$OMP                           VAR%INDCOL)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%IND1, &
!$OMP                                      VAR%INDCOL)
#endif
      END SUBROUTINE COPYIN_TRANSIT


!****************** SUBROUTINE COPYIN_KERNEL_TRUNCATION ****************
!
!***********************************************************************

      SUBROUTINE COPYIN_KERNEL_TRUNCATION(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_KERNEL_TRUNCATION
#endif
      TYPE(kernel_truncation) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%FINE_GRID)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%FINE_GRID)
         CALL COPYIN_TYPED_VAR(VAR%FINE_GRID)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%FINE_GRID)
      ENDIF
      IF (ASSOCIATED(VAR%COARSE_GRID)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%COARSE_GRID)
         CALL COPYIN_TYPED_VAR(VAR%COARSE_GRID)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COARSE_GRID)
      ENDIF
      IF (ASSOCIATED(VAR%TRANSLATION_TABLE)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%TRANSLATION_TABLE)
         CALL COPYIN_TYPED_VAR(VAR%TRANSLATION_TABLE)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%TRANSLATION_TABLE)
      ENDIF
      END SUBROUTINE COPYIN_KERNEL_TRUNCATION


!****************** SUBROUTINE COPYIN_COULOMB_POTENTIAL ****************
!
!***********************************************************************

      SUBROUTINE COPYIN_COULOMB_POTENTIAL(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_COULOMB_POTENTIAL
#endif
      TYPE(coulomb_potential) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! derived type pointers
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%KERNEL_TRUNCATE)
      CALL COPYIN_TYPED_VAR(VAR%KERNEL_TRUNCATE)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%KERNEL_TRUNCATE)
      END SUBROUTINE COPYIN_COULOMB_POTENTIAL


!****************** SUBROUTINE COPYIN_WAVEDES **************************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVEDES(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVEDES
#endif
      TYPE(wavedes) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%NPLWKP, &
!$OMP                           VAR%NGVECTOR, &
!$OMP                           VAR%NGVECTOR_POS, &
!$OMP                           VAR%NPLWKP_TOT, &
!$OMP                           VAR%NB_TOTK, &
!$OMP                           VAR%PL_INDEX, &
!$OMP                           VAR%PL_COL, &
!$OMP                           VAR%NPRO_POS, &
!$OMP                           VAR%LMMAX, &
!$OMP                           VAR%LMBASE, &
!$OMP                           VAR%NITYP, &
!$OMP                           VAR%ITYP, &
!$OMP                           VAR%NT_GLOBAL, &
!$OMP                           VAR%VKPT, &
!$OMP                           VAR%WTKPT, &
!$OMP                           VAR%NINDPW, &
!$OMP                           VAR%LUSEINV, &
!$OMP                           VAR%NINDPW_INV, &
!$OMP                           VAR%FFTSCA, &
!$OMP                           VAR%MAP_TO_FULL, &
!$OMP                           VAR%IGX, &
!$OMP                           VAR%IGY, &
!$OMP                           VAR%IGZ, &
!$OMP                           VAR%DATAKE, &
!$OMP                           VAR%AT_GAMMA)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NPLWKP, &
!$OMP                                      VAR%NGVECTOR, &
!$OMP                                      VAR%NGVECTOR_POS, &
!$OMP                                      VAR%NPLWKP_TOT, &
!$OMP                                      VAR%NB_TOTK, &
!$OMP                                      VAR%PL_INDEX, &
!$OMP                                      VAR%PL_COL, &
!$OMP                                      VAR%NPRO_POS, &
!$OMP                                      VAR%LMMAX, &
!$OMP                                      VAR%LMBASE, &
!$OMP                                      VAR%NITYP, &
!$OMP                                      VAR%ITYP, &
!$OMP                                      VAR%NT_GLOBAL, &
!$OMP                                      VAR%VKPT, &
!$OMP                                      VAR%WTKPT, &
!$OMP                                      VAR%NINDPW, &
!$OMP                                      VAR%LUSEINV, &
!$OMP                                      VAR%NINDPW_INV, &
!$OMP                                      VAR%FFTSCA, &
!$OMP                                      VAR%MAP_TO_FULL, &
!$OMP                                      VAR%IGX, &
!$OMP                                      VAR%IGY, &
!$OMP                                      VAR%IGZ, &
!$OMP                                      VAR%DATAKE, &
!$OMP                                      VAR%AT_GAMMA)
#endif
      ! native type pointers if associated
      IF (ASSOCIATED(VAR%COMM)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_INTER)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_INTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_INB)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_INB)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_KINTER)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_KINTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_KIN)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_KIN)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_SHMEM)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_SHMEM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_intra_node)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_intra_node)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_inter_node)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_inter_node)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%GRID)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%GRID)
         CALL COPYIN_TYPED_VAR(VAR%GRID)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%GRID)
      ENDIF
      END SUBROUTINE COPYIN_WAVEDES


!****************** SUBROUTINE COPYIN_WAVEDES1 *************************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVEDES1(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVEDES1
#endif
      TYPE(wavedes1) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%PL_INDEX, &
!$OMP                           VAR%PL_COL, &
!$OMP                           VAR%NPRO_POS, &
!$OMP                           VAR%LMMAX, &
!$OMP                           VAR%LMBASE, &
!$OMP                           VAR%NITYP, &
!$OMP                           VAR%ITYP, &
!$OMP                           VAR%NT_GLOBAL, &
!$OMP                           VAR%NINDPW, &
!$OMP                           VAR%NINDPW_INV, &
!$OMP                           VAR%FFTSCA, &
!$OMP                           VAR%MAP_TO_FULL, &
!$OMP                           VAR%IGX, &
!$OMP                           VAR%IGY, &
!$OMP                           VAR%IGZ, &
!$OMP                           VAR%VKPT, &
!$OMP                           VAR%DATAKE)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%PL_INDEX, &
!$OMP                                      VAR%PL_COL, &
!$OMP                                      VAR%NPRO_POS, &
!$OMP                                      VAR%LMMAX, &
!$OMP                                      VAR%LMBASE, &
!$OMP                                      VAR%NITYP, &
!$OMP                                      VAR%ITYP, &
!$OMP                                      VAR%NT_GLOBAL, &
!$OMP                                      VAR%NINDPW, &
!$OMP                                      VAR%NINDPW_INV, &
!$OMP                                      VAR%FFTSCA, &
!$OMP                                      VAR%MAP_TO_FULL, &
!$OMP                                      VAR%IGX, &
!$OMP                                      VAR%IGY, &
!$OMP                                      VAR%IGZ, &
!$OMP                                      VAR%VKPT, &
!$OMP                                      VAR%DATAKE)
#endif
      ! native type pointers if associated
      IF (ASSOCIATED(VAR%COMM)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_INTER)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_INTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_INB)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_INB)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_KINTER)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_KINTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_KIN)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_KIN)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_SHMEM)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_SHMEM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_intra_node)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_intra_node)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%COMM_inter_node)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%COMM_inter_node)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%GRID)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%GRID)
         CALL COPYIN_TYPED_VAR(VAR%GRID)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%GRID)
      ENDIF
      END SUBROUTINE COPYIN_WAVEDES1


!****************** SUBROUTINE COPYIN_WAVESPIN *************************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVESPIN(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVESPIN
#endif
      TYPE(wavespin) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELTOT needs to be created first, the same holds for FERWE and FERTOT
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%FERTOT, &
!$OMP                           VAR%AUXTOT, &
!$OMP                           VAR%CELTOT)
!$OMP TARGET ENTER DATA MAP(TO: VAR%FERWE, &
!$OMP                           VAR%CELEN)
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%CPROJ)
      ENDIF
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
      ! for crayhip we have to this explicitly after copyin was called
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%FERTOT, &
!$OMP                                      VAR%FERWE, &
!$OMP                                      VAR%AUXTOT, &
!$OMP                                      VAR%CELTOT, &
!$OMP                                      VAR%CELEN)
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: VAR%CPROJ)
      ENDIF
#endif
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%WDES)
         CALL COPYIN_TYPED_VAR(VAR%WDES)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%WDES)
      ENDIF
      END SUBROUTINE COPYIN_WAVESPIN


!****************** SUBROUTINE COPYIN_WAVEFUNA *************************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVEFUNA(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVEFUNA
#endif
      TYPE(wavefuna) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      ! native type pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR%CPTWFP, &
!$OMP                           VAR%CPROJ, &
!$OMP                           VAR%CW_RED, &
!$OMP                           VAR%CPROJ_RED)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CPTWFP, &
!$OMP                                      VAR%CPROJ, &
!$OMP                                      VAR%CW_RED, &
!$OMP                                      VAR%CPROJ_RED)
#endif
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES1)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%WDES1)
         CALL COPYIN_TYPED_VAR(VAR%WDES1)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%WDES1)
      ENDIF
      END SUBROUTINE COPYIN_WAVEFUNA


!****************** SUBROUTINE COPYIN_WAVEFUN1 *************************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVEFUN1(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVEFUN1
#endif
      TYPE(wavefun1) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      !native type pointers if associated
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%CPTWFP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%CPROJ)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CPROJ)
      ENDIF
      IF (ASSOCIATED(VAR%CR)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%CR)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CR)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES1)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%WDES1)
         CALL COPYIN_TYPED_VAR(VAR%WDES1)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%WDES1)
      ENDIF
      END SUBROUTINE COPYIN_WAVEFUN1


!****************** SUBROUTINE COPYIN_WAVEFUN1_ARRAY *******************
!
!***********************************************************************

      SUBROUTINE COPYIN_WAVEFUN1_ARRAY(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_WAVEFUN1_ARRAY
#endif
      INTEGER  :: NB
      TYPE(wavefun1), TARGET :: VAR(:)
FTNW3 TYPE(wavefun1), POINTER :: CURRENT_ENTRY => NULL()
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      DO NB=1,SIZE(VAR)
        !derived type pointers if associated
        IF (ASSOCIATED(VAR(NB)%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR(NB)%CPTWFP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR(NB)%CPTWFP)
        ENDIF
        IF (ASSOCIATED(VAR(NB)%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR(NB)%CPROJ)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR(NB)%CPROJ)
        ENDIF
        IF (ASSOCIATED(VAR(NB)%CR)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR(NB)%CR)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR(NB)%CR)
        ENDIF
        ! derived type pointers if associated
        IF (ASSOCIATED(VAR(NB)%WDES1)) THEN
FTNW3      CURRENT_ENTRY => VAR(NB)
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: CURRENT_ENTRY%WDES1)
           CALL COPYIN_TYPED_VAR(VAR(NB)%WDES1)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR(NB)%WDES1)
        ENDIF
      ENDDO
      END SUBROUTINE COPYIN_WAVEFUN1_ARRAY


!****************** SUBROUTINE COPYIN_POTCARS **************************
!
!***********************************************************************

      SUBROUTINE COPYIN_POTCARS(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_POTCARS
#endif
      TYPE(potcar) :: VAR(:)
      INTEGER :: NT
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      !native type pointers if associated
      ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
      DO NT = 1, SIZE(VAR)
!$OMP TARGET ENTER DATA MAP(TO: VAR(NT)%PSP, &
!$OMP                           VAR(NT)%PSPCOR, &
!$OMP                           VAR(NT)%PSPRHO, &
!$OMP                           VAR(NT)%LPS)
#ifdef IFXWAR2
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR(NT)%PSP, &
!$OMP                                      VAR(NT)%PSPCOR, &
!$OMP                                      VAR(NT)%PSPRHO, &
!$OMP                                      VAR(NT)%LPS)
#endif
      ENDDO
      END SUBROUTINE COPYIN_POTCARS


!****************** SUBROUTINE COPYIN_KPOINTS_FULL *********************
!
!***********************************************************************

      SUBROUTINE COPYIN_KPOINTS_FULL(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_KPOINTS_FULL
#endif
      TYPE(skpoints_full) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! native type variables
!$OMP TARGET ENTER DATA MAP(TO: VAR)
      !native type pointers if associated
      IF (ASSOCIATED(VAR%NG_INDEX)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%NG_INDEX)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NG_INDEX)
      ENDIF
      IF (ASSOCIATED(VAR%PHASE)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%PHASE)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%PHASE)
      ENDIF
      IF (ASSOCIATED(VAR%LINV)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%LINV)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%LINV)
      ENDIF
      IF (ASSOCIATED(VAR%NEQUIV)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%NEQUIV)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NEQUIV)
      ENDIF
      IF (ASSOCIATED(VAR%SPINFLIP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%SPINFLIP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%SPINFLIP)
      ENDIF
      IF (ASSOCIATED(VAR%ROTMAP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%ROTMAP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%ROTMAP)
      ENDIF
      IF (ASSOCIATED(VAR%ISYMOP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%ISYMOP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%ISYMOP)
      ENDIF
      IF (ASSOCIATED(VAR%MAP_INTO_BZ)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%MAP_INTO_BZ)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%MAP_INTO_BZ)
      ENDIF
      IF (ASSOCIATED(VAR%VKPT)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%VKPT)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%VKPT)
      ENDIF
#ifndef oldsym
      IF (ASSOCIATED(VAR%RSSYMOP)) THEN
!$OMP TARGET ENTER DATA MAP(TO: VAR%RSSYMOP)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%RSSYMOP)
      ENDIF
#endif \\ oldsym
      END SUBROUTINE COPYIN_KPOINTS_FULL


!****************** SUBROUTINE COPYIN_GREENSFDES ***********************
!
!***********************************************************************

      SUBROUTINE COPYIN_GREENSFDES(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_GREENSFDES
#endif
      TYPE(greensfdes) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR)
!$OMP TARGET ENTER DATA MAP(TO: VAR%NIONS_CORE, &
!$OMP                           VAR%NITYP, &
!$OMP                           VAR%ITYP, &
!$OMP                           VAR%NT_GLOBAL, &
!$OMP                           VAR%NI_GLOBAL, &
!$OMP                           VAR%POSION, &
!$OMP                           VAR%NPRO_LMMAX, &
!$OMP                           VAR%NPRO_LMBASE, &
!$OMP                           VAR%NPRO_POS, &
!$OMP                           VAR%NPRO_ENTRIES, &
!$OMP                           VAR%NLM_LMMAX, &
!$OMP                           VAR%NLM_POS, &
!$OMP                           VAR%NLM_ENTRIES, &
!$OMP                           VAR%COMM, &
!$OMP                           VAR%LUSEINV, &
!$OMP                           VAR%NGVECTOR_INV, &
!$OMP                           VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$OMP                           VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$OMP                           VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$OMP                           VAR%MAP_TO_FULL)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%NIONS_CORE, &
!$OMP                                      VAR%NITYP, &
!$OMP                                      VAR%ITYP, &
!$OMP                                      VAR%NT_GLOBAL, &
!$OMP                                      VAR%NI_GLOBAL, &
!$OMP                                      VAR%POSION, &
!$OMP                                      VAR%NPRO_LMMAX, &
!$OMP                                      VAR%NPRO_LMBASE, &
!$OMP                                      VAR%NPRO_POS, &
!$OMP                                      VAR%NPRO_ENTRIES, &
!$OMP                                      VAR%NLM_LMMAX, &
!$OMP                                      VAR%NLM_POS, &
!$OMP                                      VAR%NLM_ENTRIES, &
!$OMP                                      VAR%COMM, &
!$OMP                                      VAR%LUSEINV, &
!$OMP                                      VAR%NGVECTOR_INV, &
!$OMP                                      VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$OMP                                      VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$OMP                                      VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$OMP                                      VAR%MAP_TO_FULL)
#endif
      END SUBROUTINE COPYIN_GREENSFDES


!****************** SUBROUTINE COPYIN_XC *******************************
!
!***********************************************************************

      SUBROUTINE COPYIN_XC(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYIN_XC
#endif
      TYPE(xc_info) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ALLOW .OR. IS_ON_DEVICE(VAR)) RETURN
      !native type variables and pointers
!$OMP TARGET ENTER DATA MAP(TO: VAR)
!$OMP TARGET ENTER DATA MAP(TO: VAR%ID, &
!$OMP                           VAR%FAMILY, &
!$OMP                           VAR%KIND, &
!$OMP                           VAR%ID_LDA_TABLE, &
!$OMP                           VAR%NPARAM, &
!$OMP                           VAR%COEFF, &
!$OMP                           VAR%PARAM, &
!$OMP                           VAR%NAME)
#ifdef IFXWAR2
      ! Map derived type pointers with always, alloc to ensure update of pointee
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%ID, &
!$OMP                                      VAR%FAMILY, &
!$OMP                                      VAR%KIND, &
!$OMP                                      VAR%ID_LDA_TABLE, &
!$OMP                                      VAR%NPARAM, &
!$OMP                                      VAR%COEFF, &
!$OMP                                      VAR%PARAM, &
!$OMP                                      VAR%NAME)
#endif
      END SUBROUTINE COPYIN_XC


!****************** SUBROUTINE CREATE_WAVESPIN *************************
!
!***********************************************************************

      SUBROUTINE CREATE_WAVESPIN(VAR)
! in CCE19 we have to use INLINENEVER to avoid inlining
#ifdef CRAYWAR1
!DIR$ INLINENEVER CREATE_WAVESPIN
#endif
      TYPE(wavespin) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
      !note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
      !CELTOT needs to be created first, the same holds for FERWE and FERTOT
!$OMP TARGET ENTER DATA MAP(ALLOC: VAR)
!$OMP TARGET ENTER DATA MAP(ALLOC: VAR%FERTOT, &
!$OMP                              VAR%AUXTOT, &
!$OMP                              VAR%CELTOT)
!$OMP TARGET ENTER DATA MAP(ALLOC: VAR%FERWE, &
!$OMP                              VAR%CELEN)
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(ALLOC: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(ALLOC: VAR%CPROJ)
      ENDIF
#ifdef IFXWAR2
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%FERTOT, &
!$OMP                                      VAR%FERWE, &
!$OMP                                      VAR%AUXTOT, &
!$OMP                                      VAR%CELTOT, &
!$OMP                                      VAR%CELEN)
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
!$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%CPROJ)
      ENDIF
#endif
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES)) THEN
FTNW3 !$OMP TARGET ENTER DATA MAP(TO: VAR%WDES)
         CALL COPYIN_TYPED_VAR(VAR%WDES)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: VAR%WDES)
      ENDIF
      END SUBROUTINE CREATE_WAVESPIN


!****************** SUBROUTINE DELETE_NONL_STRUCT **********************
!
!***********************************************************************

      SUBROUTINE DELETE_NONL_STRUCT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_NONL_STRUCT
#endif
      TYPE(nonl_struct) :: VAR
      ! early exit if possible
      IF (.NOT.OFFLOAD_ON) RETURN
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%VKPT_SHIFT, &
!$OMP                               VAR%CQFAK, &
!$OMP                               VAR%QPROJ, &
!$OMP                               VAR%CREXP, &
!$OMP                               VAR%CPROJK)
#endif
      ! native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%VKPT_SHIFT, &
!$OMP                               VAR%CQFAK, &
!$OMP                               VAR%QPROJ, &
!$OMP                               VAR%CREXP, &
!$OMP                               VAR%CPROJK)

!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_NONL_STRUCT


!****************** SUBROUTINE DELETE_NONLR_STRUCT *********************
!
!***********************************************************************

      SUBROUTINE DELETE_NONLR_STRUCT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_NONLR_STRUCT
#endif
      TYPE(nonlr_struct) :: VAR
      ! logical variables
      INTEGER :: I
      ! early exit if possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers
      DO I=1,SIZE(VAR%BETA)
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%BETA(I))
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%BETA(I))
         CALL DELETE_TYPED_VAR(VAR%BETA(I))
      ENDDO
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%LMAX, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%CHANNELS, &
!$OMP                               VAR%PSRMAX, &
!$OMP                               VAR%RSMOOTH, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%VKPT_SHIFT, &
!$OMP                               VAR%NLIMAX, &
!$OMP                               VAR%NLI, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NLIBASE, &
!$OMP                               VAR%RPROJ, &
!$OMP                               VAR%CRREXP)
#endif

!native type pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%LMAX, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%CHANNELS, &
!$OMP                               VAR%PSRMAX, &
!$OMP                               VAR%RSMOOTH, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%VKPT_SHIFT, &
!$OMP                               VAR%NLIMAX, &
!$OMP                               VAR%NLI, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NLIBASE, &
!$OMP                               VAR%RPROJ, &
!$OMP                               VAR%CRREXP)

!native type variables
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)

      END SUBROUTINE DELETE_NONLR_STRUCT


!****************** SUBROUTINE DELETE_NONLR_PROJ ***********************
!
!***********************************************************************

      SUBROUTINE DELETE_NONLR_PROJ(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_NONLR_PROJ
#endif
      TYPE(nonlr_proj) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PSPRNL, &
!$OMP                               VAR%LPS)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PSPRNL, &
!$OMP                               VAR%LPS)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_NONLR_PROJ


!****************** SUBROUTINE DELETE_LAYOUT ***************************
!
!***********************************************************************

      SUBROUTINE DELETE_LAYOUT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_LAYOUT
#endif
      TYPE(layout) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%I2, &
!$OMP                               VAR%I3, &
!$OMP                               VAR%INDEX)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%I2, &
!$OMP                               VAR%I3, &
!$OMP                               VAR%INDEX)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_LAYOUT


!****************** SUBROUTINE DELETE_GRID_MAP *************************
!
!***********************************************************************

      SUBROUTINE DELETE_GRID_MAP(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_GRID_MAP
#endif
      TYPE(grid_map) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%N, &
!$OMP                               VAR%PTR, &
!$OMP                               VAR%TBL, &
!$OMP                               VAR%NI, &
!$OMP                               VAR%PTRI, &
!$OMP                               VAR%TBLI)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%N, &
!$OMP                               VAR%PTR, &
!$OMP                               VAR%TBL, &
!$OMP                               VAR%NI, &
!$OMP                               VAR%PTRI, &
!$OMP                               VAR%TBLI)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_GRID_MAP


!****************** SUBROUTINE DELETE_GRID_3D **************************
!
!***********************************************************************

      SUBROUTINE DELETE_GRID_3D(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_GRID_3D
#endif
      TYPE(grid_3d) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RC, VAR%IN, VAR%RL, VAR%RC_IN, VAR%IN_RL)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RC, VAR%IN, VAR%RL, VAR%RC_IN, VAR%IN_RL)
      CALL DELETE_TYPED_VAR(VAR%RC)
      CALL DELETE_TYPED_VAR(VAR%IN)
      CALL DELETE_TYPED_VAR(VAR%RL)
      CALL DELETE_TYPED_VAR(VAR%RC_IN)
      CALL DELETE_TYPED_VAR(VAR%IN_RL)
      ! native type pointer if associated
      IF (ASSOCIATED(VAR%RL_FFT)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RL_FFT)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RL_FFT)
         CALL DELETE_TYPED_VAR(VAR%RL_FFT)
      ENDIF
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%LPCTX, &
!$OMP                               VAR%LPCTY, &
!$OMP                               VAR%LPCTZ, &
!$OMP                               VAR%LPCTX_, &
!$OMP                               VAR%LPCTY_, &
!$OMP                               VAR%LPCTZ_, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%FFTWEIGHT, &
!$OMP                               VAR%NINDPWCONJG, &
!$OMP                               VAR%IND_IN_SPHERE, &
!$OMP                               VAR%COMM, &
!$OMP                               VAR%COMM_KIN, &
!$OMP                               VAR%COMM_KINTER)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%LPCTX, &
!$OMP                               VAR%LPCTY, &
!$OMP                               VAR%LPCTZ, &
!$OMP                               VAR%LPCTX_, &
!$OMP                               VAR%LPCTY_, &
!$OMP                               VAR%LPCTZ_, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%FFTWEIGHT, &
!$OMP                               VAR%NINDPWCONJG, &
!$OMP                               VAR%IND_IN_SPHERE, &
!$OMP                               VAR%COMM, &
!$OMP                               VAR%COMM_KIN, &
!$OMP                               VAR%COMM_KINTER)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_GRID_3D


!****************** SUBROUTINE DELETE_TRANSIT **************************
!
!***********************************************************************

      SUBROUTINE DELETE_TRANSIT(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_TRANSIT
#endif
      TYPE(transit) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%IND1, &
!$OMP                               VAR%INDCOL)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%IND1, &
!$OMP                               VAR%INDCOL)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_TRANSIT


!****************** SUBROUTINE DELETE_KERNEL_TRUNCATION ****************
!
!***********************************************************************

      SUBROUTINE DELETE_KERNEL_TRUNCATION(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_KERNEL_TRUNCATION
#endif
      TYPE(kernel_truncation) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%FINE_GRID)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FINE_GRID)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FINE_GRID)
         CALL DELETE_TYPED_VAR(VAR%FINE_GRID)
      ENDIF
      IF (ASSOCIATED(VAR%COARSE_GRID)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COARSE_GRID)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COARSE_GRID)
         CALL DELETE_TYPED_VAR(VAR%COARSE_GRID)
      ENDIF
      IF (ASSOCIATED(VAR%TRANSLATION_TABLE)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%TRANSLATION_TABLE)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%TRANSLATION_TABLE)
         CALL DELETE_TYPED_VAR(VAR%TRANSLATION_TABLE)
      ENDIF
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_KERNEL_TRUNCATION


!****************** SUBROUTINE DELETE_COULOMB_POTENTIAL ****************
!
!***********************************************************************

      SUBROUTINE DELETE_COULOMB_POTENTIAL(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_COULOMB_POTENTIAL
#endif
      TYPE(coulomb_potential) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%KERNEL_TRUNCATE)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%KERNEL_TRUNCATE)
      CALL DELETE_TYPED_VAR(VAR%KERNEL_TRUNCATE)
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_COULOMB_POTENTIAL


!****************** SUBROUTINE DELETE_WAVEDES **************************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVEDES(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVEDES
#endif
      TYPE(wavedes) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%COMM)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INTER)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INTER)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INB)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INB)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INB)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KINTER)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KINTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KIN)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KIN)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KIN)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_SHMEM)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_SHMEM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_intra_node)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_intra_node)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_inter_node)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_inter_node)
      ENDIF
      IF (ASSOCIATED(VAR%GRID)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%GRID)
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%GRID)
         CALL DELETE_TYPED_VAR(VAR%GRID)
      ENDIF
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NPLWKP, &
!$OMP                               VAR%NGVECTOR, &
!$OMP                               VAR%NGVECTOR_POS, &
!$OMP                               VAR%NPLWKP_TOT, &
!$OMP                               VAR%NB_TOTK, &
!$OMP                               VAR%PL_INDEX, &
!$OMP                               VAR%PL_COL, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%VKPT, &
!$OMP                               VAR%WTKPT, &
!$OMP                               VAR%NINDPW, &
!$OMP                               VAR%LUSEINV, &
!$OMP                               VAR%NINDPW_INV, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%MAP_TO_FULL, &
!$OMP                               VAR%IGX, &
!$OMP                               VAR%IGY, &
!$OMP                               VAR%IGZ, &
!$OMP                               VAR%DATAKE, &
!$OMP                               VAR%AT_GAMMA)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NPLWKP, &
!$OMP                               VAR%NGVECTOR, &
!$OMP                               VAR%NGVECTOR_POS, &
!$OMP                               VAR%NPLWKP_TOT, &
!$OMP                               VAR%NB_TOTK, &
!$OMP                               VAR%PL_INDEX, &
!$OMP                               VAR%PL_COL, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%VKPT, &
!$OMP                               VAR%WTKPT, &
!$OMP                               VAR%NINDPW, &
!$OMP                               VAR%LUSEINV, &
!$OMP                               VAR%NINDPW_INV, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%MAP_TO_FULL, &
!$OMP                               VAR%IGX, &
!$OMP                               VAR%IGY, &
!$OMP                               VAR%IGZ, &
!$OMP                               VAR%DATAKE, &
!$OMP                               VAR%AT_GAMMA)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVEDES


!****************** SUBROUTINE DELETE_WAVEDES1 *************************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVEDES1(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVEDES1
#endif
      TYPE(wavedes1) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%COMM)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INTER)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INTER)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_INB)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INB)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_INB)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KINTER)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KINTER)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_KIN)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KIN)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_KIN)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_SHMEM)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_SHMEM)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_intra_node)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_intra_node)
      ENDIF
      IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_inter_node)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%COMM_inter_node)
      ENDIF
      IF (ASSOCIATED(VAR%GRID)) THEN
IFXW2 !$OMP    TARGET EXIT DATA MAP(RELEASE: VAR%GRID)
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%GRID)
         CALL DELETE_TYPED_VAR(VAR%GRID)
      ENDIF
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PL_INDEX, &
!$OMP                               VAR%PL_COL, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%NINDPW, &
!$OMP                               VAR%NINDPW_INV, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%MAP_TO_FULL, &
!$OMP                               VAR%IGX, &
!$OMP                               VAR%IGY, &
!$OMP                               VAR%IGZ, &
!$OMP                               VAR%VKPT, &
!$OMP                               VAR%DATAKE)
#endif
      !native type variables and pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PL_INDEX, &
!$OMP                               VAR%PL_COL, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%LMMAX, &
!$OMP                               VAR%LMBASE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%NINDPW, &
!$OMP                               VAR%NINDPW_INV, &
!$OMP                               VAR%FFTSCA, &
!$OMP                               VAR%MAP_TO_FULL, &
!$OMP                               VAR%IGX, &
!$OMP                               VAR%IGY, &
!$OMP                               VAR%IGZ, &
!$OMP                               VAR%VKPT, &
!$OMP                               VAR%DATAKE)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVEDES1


!****************** SUBROUTINE DELETE_WAVESPIN *************************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVESPIN(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVESPIN
#endif
      TYPE(wavespin) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES)
         CALL DELETE_TYPED_VAR(VAR%WDES)
      ENDIF
      !native type variables and pointers
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPROJ) IF(IS_ON_DEVICE(VAR%CPROJ))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP) IF(IS_ON_DEVICE(VAR%CPTWFP))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FERWE, &
!$OMP                               VAR%FERTOT, &
!$OMP                               VAR%AUXTOT, &
!$OMP                               VAR%CELEN, &
!$OMP                               VAR%CELTOT)
#endif
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPROJ) IF(IS_ON_DEVICE(VAR%CPROJ))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP) IF(IS_ON_DEVICE(VAR%CPTWFP))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FERWE, &
!$OMP                               VAR%CELEN)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%AUXTOT, &
!$OMP                               VAR%FERTOT, &
!$OMP                               VAR%CELTOT)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVESPIN


!****************** SUBROUTINE DELETE_WAVEFUNA *************************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVEFUNA(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVEFUNA
#endif
      TYPE(wavefuna) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES1)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES1)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES1)
         CALL DELETE_TYPED_VAR(VAR%WDES1)
      ENDIF
      !native type variables and pointers
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP, &
!$OMP                               VAR%CPROJ, &
!$OMP                               VAR%CW_RED, &
!$OMP                               VAR%CPROJ_RED)
#endif
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP, &
!$OMP                               VAR%CPROJ, &
!$OMP                               VAR%CW_RED, &
!$OMP                               VAR%CPROJ_RED)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVEFUNA


!****************** SUBROUTINE DELETE_WAVEFUN1 *************************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVEFUN1(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVEFUN1
#endif
      TYPE(wavefun1) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type pointers if associated
      IF (ASSOCIATED(VAR%CPTWFP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR%CPROJ)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPROJ)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPROJ)
      ENDIF
      IF (ASSOCIATED(VAR%CR)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CR)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CR)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES1)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES1)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES1)
         CALL DELETE_TYPED_VAR(VAR%WDES1)
      ENDIF
      ! native type variables
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVEFUN1


!****************** SUBROUTINE DELETE_WAVEFUN1_ARRAY *******************
!
!***********************************************************************

      SUBROUTINE DELETE_WAVEFUN1_ARRAY(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_WAVEFUN1_ARRAY
#endif
      INTEGER  :: NB
      TYPE(wavefun1) :: VAR(:)
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type pointers if associated
      DO NB=1,SIZE(VAR)
        IF (ASSOCIATED(VAR(NB)%CPTWFP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CPTWFP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CPTWFP)
        ENDIF
        IF (ASSOCIATED(VAR(NB)%CPROJ)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CPROJ)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CPROJ)
        ENDIF
        IF (ASSOCIATED(VAR(NB)%CR)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CR)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%CR)
        ENDIF
        ! derived type pointers if associated
        IF (ASSOCIATED(VAR(NB)%WDES1)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%WDES1)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NB)%WDES1)
           CALL DELETE_TYPED_VAR(VAR(NB)%WDES1)
        ENDIF
      ENDDO
      ! native type variables
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_WAVEFUN1_ARRAY


!****************** SUBROUTINE DELETE_POTCARS **************************
!
!***********************************************************************

      SUBROUTINE DELETE_POTCARS(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_POTCARS
#endif
      TYPE(potcar) :: VAR(:)
      INTEGER :: NT
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
      ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
      DO NT = 1, SIZE(VAR)
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NT)%PSP, &
!$OMP                               VAR(NT)%PSPCOR, &
!$OMP                               VAR(NT)%PSPRHO, &
!$OMP                               VAR(NT)%LPS)
#endif
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR(NT)%PSP, &
!$OMP                               VAR(NT)%PSPCOR, &
!$OMP                               VAR(NT)%PSPRHO, &
!$OMP                               VAR(NT)%LPS)
      ENDDO
      ! native type variables
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_POTCARS


!****************** SUBROUTINE DELETE_KPOINTS_FULL *********************
!
!***********************************************************************

      SUBROUTINE DELETE_KPOINTS_FULL(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_KPOINTS_FULL
#endif
      TYPE(skpoints_full) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers if associated
      IF (ASSOCIATED(VAR%NG_INDEX).AND.IS_ON_DEVICE(VAR%NG_INDEX)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NG_INDEX)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NG_INDEX)
      ENDIF
      IF (ASSOCIATED(VAR%PHASE).AND.IS_ON_DEVICE(VAR%PHASE)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PHASE)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%PHASE)
      ENDIF
      IF (ASSOCIATED(VAR%LINV)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%LINV)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%LINV)
      ENDIF
      IF (ASSOCIATED(VAR%NEQUIV)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NEQUIV)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NEQUIV)
      ENDIF
      IF (ASSOCIATED(VAR%SPINFLIP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%SPINFLIP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%SPINFLIP)
      ENDIF
      IF (ASSOCIATED(VAR%ROTMAP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ROTMAP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ROTMAP)
      ENDIF
      IF (ASSOCIATED(VAR%ISYMOP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ISYMOP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ISYMOP)
      ENDIF
      IF (ASSOCIATED(VAR%MAP_INTO_BZ)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%MAP_INTO_BZ)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%MAP_INTO_BZ)
      ENDIF
      IF (ASSOCIATED(VAR%VKPT)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%VKPT)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%VKPT)
      ENDIF
#ifndef oldsym
      IF (ASSOCIATED(VAR%RSSYMOP)) THEN
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RSSYMOP)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%RSSYMOP)
      ENDIF
#endif \\ oldsym
      ! native type variables
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_KPOINTS_FULL


!****************** SUBROUTINE DELETE_GREENSFDES ***********************
!
!***********************************************************************

      SUBROUTINE DELETE_GREENSFDES(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_GREENSFDES
#endif
      TYPE(greensfdes) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      !native type variables and pointers
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NIONS_CORE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%NI_GLOBAL, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%NPRO_LMMAX, &
!$OMP                               VAR%NPRO_LMBASE, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%NPRO_ENTRIES, &
!$OMP                               VAR%NLM_LMMAX, &
!$OMP                               VAR%NLM_POS, &
!$OMP                               VAR%NLM_ENTRIES, &
!$OMP                               VAR%COMM, &
!$OMP                               VAR%LUSEINV, &
!$OMP                               VAR%NGVECTOR_INV, &
!$OMP                               VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$OMP                               VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$OMP                               VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$OMP                               VAR%MAP_TO_FULL)
#endif
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%NIONS_CORE, &
!$OMP                               VAR%NITYP, &
!$OMP                               VAR%ITYP, &
!$OMP                               VAR%NT_GLOBAL, &
!$OMP                               VAR%NI_GLOBAL, &
!$OMP                               VAR%POSION, &
!$OMP                               VAR%NPRO_LMMAX, &
!$OMP                               VAR%NPRO_LMBASE, &
!$OMP                               VAR%NPRO_POS, &
!$OMP                               VAR%NPRO_ENTRIES, &
!$OMP                               VAR%NLM_LMMAX, &
!$OMP                               VAR%NLM_POS, &
!$OMP                               VAR%NLM_ENTRIES, &
!$OMP                               VAR%COMM, &
!$OMP                               VAR%LUSEINV, &
!$OMP                               VAR%NGVECTOR_INV, &
!$OMP                               VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$OMP                               VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$OMP                               VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$OMP                               VAR%MAP_TO_FULL)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      END SUBROUTINE DELETE_GREENSFDES


!****************** SUBROUTINE DELETE_XC *******************************
!
!***********************************************************************

      SUBROUTINE DELETE_XC(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELETE_XC
#endif
      TYPE(xc_info), TARGET :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ALLOW .OR. .NOT.IS_ON_DEVICE(VAR)) RETURN
      !native type variables and pointers
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ID, &
!$OMP                               VAR%FAMILY, &
!$OMP                               VAR%KIND, &
!$OMP                               VAR%ID_LDA_TABLE, &
!$OMP                               VAR%NPARAM, &
!$OMP                               VAR%COEFF, &
!$OMP                               VAR%PARAM, &
!$OMP                               VAR%NAME)
#endif
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%ID, &
!$OMP                               VAR%FAMILY, &
!$OMP                               VAR%KIND, &
!$OMP                               VAR%ID_LDA_TABLE, &
!$OMP                               VAR%NPARAM, &
!$OMP                               VAR%COEFF, &
!$OMP                               VAR%PARAM, &
!$OMP                               VAR%NAME)
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR)
      ! safety check
      IF (IS_ON_DEVICE(VAR)) &
         CALL vtutor%bug("DELETE_XC: VAR still present on device.", __FILE__, __LINE__)
      END SUBROUTINE DELETE_XC


!****************** SUBROUTINE COPYOUT_WAVESPIN ************************
!
!***********************************************************************

      SUBROUTINE COPYOUT_WAVESPIN(VAR)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER COPYOUT_WAVESPIN
#endif
      TYPE(wavespin) :: VAR
      ! early exit of possible
      IF (.NOT.OFFLOAD_ON) RETURN
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR%WDES)) THEN
FTNW3 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE: VAR%WDES)
         CALL DELETE_TYPED_VAR(VAR%WDES)
      ENDIF
#ifdef IFXWAR2
      ! Extra release of pointers
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPTWFP) IF(IS_ON_DEVICE(VAR%CPTWFP))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%CPROJ) IF(IS_ON_DEVICE(VAR%CPROJ))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FERTOT, &
!$OMP                               VAR%FERWE, &
!$OMP                               VAR%AUXTOT, &
!$OMP                               VAR%CELTOT, &
!$OMP                               VAR%CELEN)
#endif
!native type variables and pointers
!$OMP TARGET EXIT DATA MAP(FROM: VAR%CPROJ) IF(IS_ON_DEVICE(VAR%CPROJ))
!$OMP TARGET EXIT DATA MAP(FROM: VAR%CPTWFP) IF(IS_ON_DEVICE(VAR%CPTWFP))
!$OMP TARGET EXIT DATA MAP(RELEASE: VAR%FERWE, &
!$OMP                               VAR%CELEN)
!$OMP TARGET EXIT DATA MAP(FROM: VAR%AUXTOT, &
!$OMP                            VAR%FERTOT, &
!$OMP                            VAR%CELTOT)      
!$OMP TARGET EXIT DATA MAP(FROM: VAR)
      END SUBROUTINE COPYOUT_WAVESPIN

    END MODULE mopenmp_deepcopy

#endif // OMP_OFFLOAD

!> @page openmp OpenMP
!! @tableofcontents
!! @section general All changes
!! All datatypes, subroutines, and functions that have been changed
!! to enable parallelization under OpenMP.
!!
!! \li ::apply_gfac
!! \li ::apply_gfac_2_encuts
!! \li ::apply_gfac_der
!! \li ::apply_gfac_encuts
!! \li ::apply_gfac_exchange
!! \li ::apply_gfac_exchange_2
!! \li augfast::calc_dllmm_trans
!! \li augfast::depsum_two_bands_rholm_trace
!! \li augfast::depsum_vector
!! \li augfast::setup_aug_des
!! \li broyden::brmix
!! \li chi_base::add_responsefunction_cache
!! \li chi_base::add_responsefunction_int
!! \li chi_base::add_xi
!! \li chi_base::add_xi_spectral
!! \li chi_base::allocate_responsefun_cache
!! \li chi_base::clean_responsefunction_cache
!! \li chi_base::clean_responsefunction_int
!! \li chi_base::deallocate_responsefun
!! \li chi_base::responsefunction
!! \li chi_base::screened_two_electron_cached
!! \li choleski::orthch
!! \li ::crrexp_mul_work_gadd
!! \li ::crexp_mul_wave
!! \li ::complex_add
!! \li ::cmplx_cmplx_cmplx_mul
!! \li ::cmplx_cmplx_real_mul
!! \li ::cmplx_real_real_mul
!! \li ::dllmm_kernel
!! \li ::eccp_nl
!! \li ::eccp_nl_fock
!! \li ::exchange_gfac
!! \li ::fexcg_
!! \li ::fftbas
!! \li ::fftbas_plan_mpi
!! \li ::fftbrc
!! \li ::fftext_mpi
!! \li ::ffthamil
!! \li ::fftmakeplan
!! \li ::fftmakeplan_mpi
!! \li ::fftwav
!! \li ::fftwav_mpi
!! \li fock::allocate_fock_handle
!! \li fock::deallocate_fock_handle
!! \li fock::fock_acc
!! \li fock::fock_charge_mu
!! \li fock::fock_force
!! \li fock::fock_handle
!! \li fock::xc_fock_reader
!! \li fock_dbl::fock_all_dblbuf
!! \li force::forhar
!! \li force::forloc
!! \li greens_real_space::transpose_g_r_response
!! \li greens_real_space::transpose_g_proj_response
!! \li greens_real_space::transpose_r_g
!! \li greens_real_space::transpose_r_proj
!! \li greens_real_space_k::calculate_response_super
!! \li greens_real_space_k::fft_g_super
!! \li hamil::eccp
!! \li hamil::eccp_tau
!! \li hamil::eccp_vec
!! \li hamil::hamiltmu
!! \li hamil::hamiltmu_c
!! \li hamil::hamiltmu_tau
!! \li hamil::hamiltmu_vec
!! \li hamil::setup_precond
!! \li hamil::simple_precond
!! \li hamil_lrf::hamiltmu_commutator
!! \li hamil_lrf::hamiltmu_lrf
!! \li hamil_lrf::lrf_commutator
!! \li hamil_lrf::lrf_hamil
!! \li ::kinhamil
!! \li ::kinhamil_c
!! \li ::kinhamil_tau
!! \li ::kinhamil_vec
!! \li main_mpi::init_mpi
!! \li main_mpi::wrt_distr
!! \li ::map_backward
!! \li ::map_forward
!! \li ::map_gather
!! \li ::map_scatter
!! \li mpimy::m_init
!! \li nonl::phase
!! \li nonl::proj1
!! \li nonl::projxyz
!! \li nonl::projxyz_wa
!! \li nonl::spher
!! \li nonl::strenl
!! \li nonlr::fornlr
!! \li nonlr::nonlr_alloc
!! \li nonlr::nonlr_alloc_crrexp
!! \li nonlr::nonlr_dealloc
!! \li nonlr::nonlr_set_single_ion
!! \li nonlr::phaser
!! \li nonlr::phaserr
!! \li nonlr::phaser_hf
!! \li nonlr::real_optlay
!! \li nonlr::real_optlay_grid
!! \li nonlr::rnlpr
!! \li nonlr::rpro1
!! \li nonlr::rpro1_hf
!! \li nonlr::rpromu
!! \li nonlr::rpromu_hf
!! \li nonlr::rspher_all
!! \li nonlr_struct_def::nonlr_struct
!! \li nonl_high::w1_projall
!! \li mopenmp::init_omp
!! \li ::overl
!! \li ::overl1
!! \li ::overl1_c
!! \li ::overl_fock
!! \li pawfock::coloumb_4term
!! \li pawfock::coloumb_4term_ps
!! \li pawfock_inter::ntyp_slater
!! \li pawfock_inter::s
!! \li pawm::set_dd_paw
!! \li ::pw_charge
!! \li ::pw_charge_trace
!! \li ::pw_norm_with_metric
!! \li ::racc0
!! \li ::racc0mu
!! \li ::racc0mu_hf
!! \li ::racc0_hf
!! \li radial::rad_lda_xc
!! \li ::real_add
!! \li ::real_cmplx_cmplx_mul
!! \li ::real_real_cmplx_mul
!! \li ::real_real_real_mul
!! \li ::rholm_kernel
!! \li ::rholm_kernel_aux
!! \li ::rholm_kernel_dgemm
!! \li ::rholm_one_center_kernel
!! \li rmm_diis::eddrmm
!! \li rmm_diis_lr::linear_response_diis
!! \li rot::edwav
!! \li rs_greensfunc_kernel::set_rsgf_all
!! \li rs_greensfunc_kernel::rs_coulomb_green_func
!! \li rs_greensfunc_kernel::rs_greensfunc_type
!! \li screened_2e::determine_slot
!! \li screened_2e::determine_slot_inter
!! \li screened_2e::determine_slot_inter_weight
!! \li screened_2e::integrate_w_2e_simple
!! \li screened_2e::integrate_w_2e_spectral
!! \li screened_2e::integrate_w_2e_spectral_imag
!! \li screened_2e::qp_shift
!! \li screened_2e::qp_shift_pade
!! \li ::setdij_
!! \li ::start_profiling
!! \li ::stop_profiling
!! \li subrot::eddiag
!! \li ::truncate_high_frequency_one
!! \li twoelectron4o::apply_phase
!! \li us::depsum
!! \li us::fordep
!! \li ::vhamil
!! \li ::vhamil_trace
!! \li wave::delwav_omp
!! \li wave::newwav_omp
!! \li wave_cacher::allocate_cacher
!! \li wave_cacher::deallocate_cacher
!! \li wave_cacher::find_new_cacher
!! \li wave_cacher::remove_cacher
!! \li wave_cacher::store_cacher
!! \li wave_cacher::store_gw_acc
!! \li wave_high::cnorma
!! \li wave_high::w1_copy
!! \li wave_high::w1_daxpy
!! \li wave_high::w1_dot
!! \li wave_high::w1_dscal
!! \li wave_high::w1_gaxpy
!! \li ::work_mul_crexp
!! \li xi::calculate_xi
!!
!! @section simple Simple loop parallelism
!! In the following cases only basic loop parallelism under OpenMP was added.
!! Sometimes in conjunction with a reduction operation, and sometimes
!! conditional upon a certain loop size.
!! These loops often run over all points in a real or reciprocal space grid,
!! or over all basis vectors in an orbital.
!!
!! \li ::apply_gfac
!! \li ::apply_gfac_2_encuts
!! \li ::apply_gfac_der
!! \li ::apply_gfac_encuts
!! \li ::apply_gfac_exchange
!! \li ::apply_gfac_exchange_2
!! \li augfast::depsum_vector
!! \li broyden::brmix
!! \li ::crrexp_mul_work_gadd
!! \li ::crexp_mul_wave
!! \li ::complex_add
!! \li ::dllmm_kernel
!! \li ::cmplx_cmplx_cmplx_mul
!! \li ::cmplx_cmplx_real_mul
!! \li ::cmplx_real_real_mul
!! \li ::eccp_nl
!! \li ::eccp_nl_fock
!! \li ::exchange_gfac
!! \li ::fexcg_
!! \li ::ffthamil
!! \li force::forhar
!! \li force::forloc
!! \li ::kinhamil
!! \li nonl::phase
!! \li nonl::strenl
!! \li pawfock::coloumb_4term
!! \li pawfock::coloumb_4term_ps
!! \li ::pw_charge
!! \li ::pw_charge_trace
!! \li ::pw_norm_with_metric
!! \li radial::rad_lda_xc
!! \li ::real_add
!! \li ::real_cmplx_cmplx_mul
!! \li ::real_real_cmplx_mul
!! \li ::real_real_real_mul
!! \li ::rholm_kernel
!! \li ::rholm_kernel_aux
!! \li ::rholm_kernel_dgemm
!! \li ::rholm_one_center_kernel
!! \li rs_greensfunc_kernel::rs_coulomb_green_func
!! \li ::truncate_high_frequency_one
!! \li twoelectron4o::apply_phase
!! \li ::vhamil
!! \li ::vhamil_trace
!! \li ::work_mul_crexp
