#include "symbol.inc"
#ifdef ACC_OFFLOAD
!***********************************************************************
!
!***********************************************************************

  MODULE mopenacc
  USE prec

  USE mopenacc_struct_def
  USE moffload_struct_def, ONLY : OFFLOAD_ALLOW, OFFLOAD_ON, &
                                  OFFLOAD_NUM_DEVICES, OFFLOAD_MAX_RANKS_PER_DEVICE

  USE mdevice_info

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INIT_OFFLOAD, &
            ACC_SET_ASYNC_Q,ACC_SYNC_ASYNC_Q, &
            GET_FREE_MEM_DEVICE,ACC_GET_BATCH_SIZE, &
            SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX

  INTEGER :: ACC_MAX_STREAMS=1024

#ifdef ACC_DEBUG
  LOGICAL, PUBLIC :: ACC_DISABLE_SET_ASYNC_Q=.FALSE.
#endif

  INTERFACE ACC_SYNC_ASYNC_Q
    MODULE PROCEDURE SYNC_ASYNC_Q_ALL
    MODULE PROCEDURE SYNC_ASYNC_Q_RANGE
    MODULE PROCEDURE SYNC_ASYNC_Q_MASKED
  END INTERFACE ACC_SYNC_ASYNC_Q

  INTERFACE
    SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(SPLIT_CMPLX, TO_ADD)
!$ACC ROUTINE SEQ
    USE prec
    COMPLEX(q) :: SPLIT_CMPLX, TO_ADD
    END SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX
  END INTERFACE

  CONTAINS

!****************** SUBROUTINE INIT_OFFLOAD ****************************
!
!***********************************************************************

  SUBROUTINE INIT_OFFLOAD(IO)
    USE base
    USE main_mpi
!$  USE mopenmp_struct_def

    USE openacc, ONLY : ACC_DEVICE_KIND,ACC_GET_DEVICE_TYPE,ACC_DEVICE_NVIDIA,ACC_DEVICE_HOST, &
                        ACC_GET_NUM_DEVICES,ACC_SET_DEVICE_NUM,ACC_DEVICE_KIND,ACC_INIT

    TYPE(in_struct) IO
    ! local variables
    INTEGER(KIND=ACC_DEVICE_KIND) :: DEVICE_TYPE
    INTEGER :: DEVICE_NUM,NUM_DEVICES
    INTEGER :: NBLOCK_RANKS_PER_DEVICE=1

    PROFILING_START('init_offload')

    DEVICE_NUM=0

    DEVICE_TYPE=ACC_GET_DEVICE_TYPE()

    IF (DEVICE_TYPE==ACC_DEVICE_HOST) THEN
       OFFLOAD_ALLOW=.FALSE.
#ifdef USENCCL
       LUSENCCL=.FALSE.
#endif
       PROFILING_STOP('init_offload')
       RETURN
    ENDIF

#ifdef MPI
    CALL CHECK_MPI()
#endif

    OFFLOAD_ALLOW=.TRUE.

    IF (DEVICE_TYPE==ACC_DEVICE_NVIDIA) THEN
       NUM_DEVICES=ACC_GET_NUM_DEVICES(ACC_DEVICE_NVIDIA)
#ifdef MPI
       DEVICE_NUM=MOD((COMM_WORLD%NODE_ME-1)/NBLOCK_RANKS_PER_DEVICE,NUM_DEVICES)
#endif
!$OMP PARALLEL
       CALL ACC_SET_DEVICE_NUM(DEVICE_NUM,ACC_DEVICE_NVIDIA)
!$OMP END PARALLEL
    ENDIF

!$OMP PARALLEL
    CALL ACC_INIT(DEVICE_TYPE)
!$OMP END PARALLEL

#ifdef MPI
    OFFLOAD_NUM_DEVICES = GET_NUM_DEVICES_GLOBAL(DEVICE_NUM)

    OFFLOAD_MAX_RANKS_PER_DEVICE = GET_NUM_RANKS_ON_DEVICE(DEVICE_NUM)
    CALL M_max_i(COMM_inter_node_world, OFFLOAD_MAX_RANKS_PER_DEVICE, 1)
#else
    OFFLOAD_NUM_DEVICES = NUM_DEVICES
    OFFLOAD_MAX_RANKS_PER_DEVICE = 1
#endif

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

       ! if not write a warning that NCCL has been switched off (completely for now) ...
       IF (IO%IU0>=0 .AND. .NOT. LUSENCCL) THEN
          CALL vtutor%warning('INIT_OFFLOAD: several MPI-ranks need to share a GPU, which is not supported ' // &
                              'by NCCL. The use of NCCL will be switched off. To avoid this, reduce '// &
                              'the number of MPI-ranks: #-of-ranks <= #-of-GPUs (on every node!).')
#if ! PGI_BEFORE_25_9 && defined(CUSOLVERMP)
          CALL vtutor%warning('INIT_OFFLOAD: cuSOLVERmp requires NCCL. The use of cuSOLVERmp will be switched off.')
#endif
       ENDIF
    ENDIF
#endif

#ifdef debug
    CALL ACC_SET_WAIT_ON_EXCEPTION(.TRUE.)
#endif

    PROFILING_STOP('init_offload')

  END SUBROUTINE INIT_OFFLOAD

#ifdef debug
  SUBROUTINE SLEEP_ON_EXCEPTION()
     INTERFACE
     FUNCTION GETPID() BIND(C, NAME='getpid')
        USE ISO_C_BINDING
        INTEGER(C_INT) :: GETPID
     END FUNCTION
     END INTERFACE

        PRINT *, "EXCEPTION ENCOUNTERED! waiting for the debugger to be attached to PID ", getpid()
        DO WHILE (.TRUE.)
           CALL SLEEP(1)
        END DO
  END SUBROUTINE SLEEP_ON_EXCEPTION

  SUBROUTINE ACC_SET_WAIT_ON_EXCEPTION(ENABLE)
    USE ISO_C_BINDING
    IMPLICIT NONE
    LOGICAL :: ENABLE
    INTEGER :: IERR
    INTERFACE
    FUNCTION ATEXIT(FPTR) BIND(C, NAME='atexit')
       USE ISO_C_BINDING
       TYPE(C_FUNPTR), VALUE, INTENT(IN) :: FPTR
       INTEGER(c_int) :: ATEXIT
    END FUNCTION
    END INTERFACE

    TYPE(C_FUNPTR) :: fptr
    fptr = c_funloc(SLEEP_ON_EXCEPTION)
!$OMP PARALLEL
    IF ( ENABLE ) IERR = ATEXIT(fptr)
!$OMP END PARALLEL

  END SUBROUTINE ACC_SET_WAIT_ON_EXCEPTION
#endif // debug

!****************** SUBROUTINE ACC_SET_ASYNC_Q *************************
!
!***********************************************************************

  SUBROUTINE ACC_SET_ASYNC_Q(VAR)
    INTEGER :: VAR
#ifdef ACC_DEBUG
    IF (ACC_DISABLE_SET_ASYNC_Q) RETURN
#endif
    IF (OFFLOAD_ON) ACC_ASYNC_Q=MOD(VAR-1,ACC_MAX_STREAMS)+1
  END SUBROUTINE ACC_SET_ASYNC_Q


!****************** SUBROUTINE ACC_SYNC_ASYNC_Q ************************
!
!***********************************************************************

  SUBROUTINE SYNC_ASYNC_Q_ALL
!!#define use_cudaDeviceSynchronize
#ifdef use_cudaDeviceSynchronize
    USE cudafor, ONLY : cudaDeviceSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
    IF (OFFLOAD_ON) THEN
       ierr = cudaDeviceSynchronize()
       IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_ALL: internal error:',cudaGetErrorString(cudaGetLastError())
    ENDIF
    RETURN
#endif
!$ACC WAIT IF(OFFLOAD_ON)
  END SUBROUTINE SYNC_ASYNC_Q_ALL

  SUBROUTINE SYNC_ASYNC_Q_RANGE(NQ)
!!#define use_cudaStreamSynchronize
#ifdef use_cudaStreamSynchronize
    USE openacc, ONLY : acc_get_cuda_stream
    USE cudafor, ONLY : cudaStreamSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
#endif
    INTEGER :: NQ
    ! local variables
    INTEGER :: I
    IF (.NOT.OFFLOAD_ON) RETURN
    DO I=1,NQ
#ifdef use_cudaStreamSynchronize
       ierr = cudaStreamSynchronize(acc_get_cuda_stream(I))
       IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_RANGE: internal error:',cudaGetErrorString(cudaGetLastError())
#else
!$ACC WAIT(I)
#endif
    ENDDO
  END SUBROUTINE SYNC_ASYNC_Q_RANGE

  SUBROUTINE SYNC_ASYNC_Q_MASKED(NQ,MASK)
!!#define use_cudaStreamSynchronize
#ifdef use_cudaStreamSynchronize
    USE openacc, ONLY : acc_get_cuda_stream
    USE cudafor, ONLY : cudaStreamSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
#endif
    INTEGER :: NQ
    LOGICAL :: MASK(:)
    ! local variables
    INTEGER :: I
    IF (.NOT.OFFLOAD_ON) RETURN
    DO I=1,NQ
       IF (MASK(I)) THEN
#ifdef use_cudaStreamSynchronize
          ierr = cudaStreamSynchronize(acc_get_cuda_stream(I))
          IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_MASKED: internal error:',cudaGetErrorString(cudaGetLastError())
#else
!$ACC WAIT(I)
#endif
       ENDIF
    ENDDO
  END SUBROUTINE SYNC_ASYNC_Q_MASKED


!****************** FUNCTION ACC_GET_BATCH_SIZE ************************
!
!***********************************************************************

  FUNCTION ACC_GET_BATCH_SIZE(RSIZE,RBASE)
    REAL(q) :: RSIZE
    REAL(q), OPTIONAL :: RBASE

    INTEGER :: ACC_GET_BATCH_SIZE

    ! local variables
    REAL(q) :: MFREE
    INTEGER :: IERR

    REAL(q), PARAMETER :: MAXIMUM_FRACTION = 0.8_q

    ACC_GET_BATCH_SIZE=1

    ! when we have (at most) a single MPI-rank per device
    IF (OFFLOAD_MAX_RANKS_PER_DEVICE == 1) THEN
       ! get available memory (in bytes)
       IERR = GET_FREE_MEM_DEVICE(MFREE)
       IF (IERR == 0) THEN
          ! optionally remove a chunk of RBASE bytes
          IF (PRESENT(RBASE)) MFREE = MFREE - RBASE
          ! only take a fraction of the remaining available memory
          MFREE = MFREE * MAXIMUM_FRACTION
          ! and determine the corresponding batch size
          ACC_GET_BATCH_SIZE = MAX(FLOOR(MFREE/RSIZE), 1)
       ENDIF
    ENDIF
  END FUNCTION ACC_GET_BATCH_SIZE

  END MODULE mopenacc


!****************** SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX *******
!
! This routine takes one complex argument and accepts this as two real
! arguments for the real and imaginary parts, respectively. It then
! performs two separate OpenACC atomic adds, one for each part.
!
!***********************************************************************

  SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(SPLIT_CMPLX,TO_ADD)
!$ACC ROUTINE SEQ
    USE prec
    REAL(q),DIMENSION(2) :: SPLIT_CMPLX
    COMPLEX(q) :: TO_ADD
!$ACC ATOMIC UPDATE
    SPLIT_CMPLX(1)=SPLIT_CMPLX(1)+REAL(TO_ADD,q) ! real part
!$ACC END ATOMIC
!$ACC ATOMIC UPDATE
    SPLIT_CMPLX(2)=SPLIT_CMPLX(2)+AIMAG(TO_ADD) ! imaginary part
!$ACC END ATOMIC
  END SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX


!***********************************************************************
!
!***********************************************************************

  MODULE mopenacc_deepcopy

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
    TYPE(nonl_struct) :: VAR
    ! early exit if possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%POSION, &
!$ACC&                  VAR%VKPT_SHIFT, &
!$ACC&                  VAR%CQFAK, &
!$ACC&                  VAR%QPROJ, &
!$ACC&                  VAR%CREXP) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_NONL_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_STRUCT *********************
!
!***********************************************************************

  SUBROUTINE COPYIN_NONLR_STRUCT(VAR)
    TYPE(nonlr_struct) :: VAR
    ! local variables
    INTEGER :: I
    ! early exit if possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%LMAX, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%CHANNELS, &
!$ACC&                  VAR%PSRMAX, &
!$ACC&                  VAR%RSMOOTH, &
!$ACC&                  VAR%POSION, &
!$ACC&                  VAR%VKPT_SHIFT, &
!$ACC&                  VAR%NLIMAX, &
!$ACC&                  VAR%NLI, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NLIBASE, &
!$ACC&                  VAR%RPROJ, &
!$ACC&                  VAR%CRREXP) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers
    DO I=1,SIZE(VAR%BETA)
       CALL COPYIN_TYPED_VAR(VAR%BETA(I))
#if PGI_BEFORE_24_5
       CALL ACC_ATTACH_ASYNC(VAR%BETA(I),ACC_ASYNC_Q)
#endif
    ENDDO
  END SUBROUTINE COPYIN_NONLR_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_PROJ ***********************
!
!***********************************************************************

  SUBROUTINE COPYIN_NONLR_PROJ(VAR)
    TYPE(nonlr_proj) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%PSPRNL, &
!$ACC&                  VAR%LPS) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_NONLR_PROJ


!****************** SUBROUTINE COPYIN_LAYOUT ***************************
!
!***********************************************************************

  SUBROUTINE COPYIN_LAYOUT(VAR)
    TYPE(layout) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%I2, &
!$ACC&                  VAR%I3, &
!$ACC&                  VAR%INDEX) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_LAYOUT


!****************** SUBROUTINE COPYIN_GRID_MAP *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_GRID_MAP(VAR)
    TYPE(grid_map) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%N, &
!$ACC&                  VAR%PTR, &
!$ACC&                  VAR%TBL, &
!$ACC&                  VAR%NI, &
!$ACC&                  VAR%PTRI, &
!$ACC&                  VAR%TBLI) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_GRID_MAP


!****************** SUBROUTINE COPYIN_GRID_3D **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_GRID_3D(VAR)
    TYPE(grid_3d) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%LPCTX, &
!$ACC&                  VAR%LPCTY, &
!$ACC&                  VAR%LPCTZ, &
!$ACC&                  VAR%LPCTX_, &
!$ACC&                  VAR%LPCTY_, &
!$ACC&                  VAR%LPCTZ_, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%FFTWEIGHT, &
!$ACC&                  VAR%NINDPWCONJG, &
!$ACC&                  VAR%IND_IN_SPHERE, &
!$ACC&                  VAR%COMM, &
!$ACC&                  VAR%COMM_KIN, &
!$ACC&                  VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)

    ! native type pointer if associated
    IF (ASSOCIATED(VAR%RL_FFT)) THEN
!!$ACC ENTER DATA COPYIN(VAR%RL_FFT) ASYNC(ACC_ASYNC_Q)
       CALL COPYIN_TYPED_VAR(VAR%RL_FFT)
       CALL ACC_ATTACH_ASYNC(VAR%RL_FFT,ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers

    CALL COPYIN_TYPED_VAR(VAR%RC)
    CALL ACC_ATTACH_ASYNC(VAR%RC,ACC_ASYNC_Q)

    CALL COPYIN_TYPED_VAR(VAR%IN)
    CALL ACC_ATTACH_ASYNC(VAR%IN,ACC_ASYNC_Q)

    CALL COPYIN_TYPED_VAR(VAR%RL)
    CALL ACC_ATTACH_ASYNC(VAR%RL,ACC_ASYNC_Q)

    CALL COPYIN_TYPED_VAR(VAR%RC_IN)
    CALL ACC_ATTACH_ASYNC(VAR%RC_IN,ACC_ASYNC_Q)

    CALL COPYIN_TYPED_VAR(VAR%IN_RL)
    CALL ACC_ATTACH_ASYNC(VAR%IN_RL,ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_GRID_3D


!****************** SUBROUTINE COPYIN_TRANSIT **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_TRANSIT(VAR)
    TYPE(transit) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%IND1, &
!$ACC&                  VAR%INDCOL) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_TRANSIT

  SUBROUTINE COPYIN_KERNEL_TRUNCATION(VAR)
    TYPE(kernel_truncation) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%FINE_GRID)) THEN
       CALL COPYIN_TYPED_VAR(VAR%FINE_GRID)
       CALL ACC_ATTACH_ASYNC(VAR%FINE_GRID,ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COARSE_GRID)) THEN
       CALL COPYIN_TYPED_VAR(VAR%COARSE_GRID)
       CALL ACC_ATTACH_ASYNC(VAR%COARSE_GRID,ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%TRANSLATION_TABLE)) THEN
       CALL COPYIN_TYPED_VAR(VAR%TRANSLATION_TABLE)
       CALL ACC_ATTACH_ASYNC(VAR%TRANSLATION_TABLE,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE

  SUBROUTINE COPYIN_COULOMB_POTENTIAL(VAR)
    TYPE(coulomb_potential) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    !
    CALL COPYIN_TYPED_VAR(VAR%KERNEL_TRUNCATE)
  END SUBROUTINE


!****************** SUBROUTINE COPYIN_WAVEDES **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEDES(VAR)
    TYPE(wavedes) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NPLWKP, &
!$ACC&                  VAR%NGVECTOR, &
!$ACC&                  VAR%NGVECTOR_POS, &
!$ACC&                  VAR%NPLWKP_TOT, &
!$ACC&                  VAR%NB_TOTK, &
!$ACC&                  VAR%PL_INDEX, &
!$ACC&                  VAR%PL_COL, &
!$ACC&                  VAR%NPRO_POS, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%NT_GLOBAL, &
!$ACC&                  VAR%VKPT, &
!$ACC&                  VAR%WTKPT, &
!$ACC&                  VAR%NINDPW, &
!$ACC&                  VAR%LUSEINV, &
!$ACC&                  VAR%NINDPW_INV, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%MAP_TO_FULL, &
!$ACC&                  VAR%IGX, &
!$ACC&                  VAR%IGY, &
!$ACC&                  VAR%IGZ, &
!$ACC&                  VAR%DATAKE, &
!$ACC&                  VAR%AT_GAMMA) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL COPYIN_TYPED_VAR(VAR%GRID)
       CALL ACC_ATTACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEDES


!****************** SUBROUTINE COPYIN_WAVEDES1 *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEDES1(VAR)
    TYPE(wavedes1) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%PL_INDEX, &
!$ACC&                  VAR%PL_COL, &
!$ACC&                  VAR%NPRO_POS, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%NT_GLOBAL, &
!$ACC&                  VAR%NINDPW, &
!$ACC&                  VAR%NINDPW_INV, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%MAP_TO_FULL, &
!$ACC&                  VAR%IGX, &
!$ACC&                  VAR%IGY, &
!$ACC&                  VAR%IGZ, &
!$ACC&                  VAR%VKPT, &
!$ACC&                  VAR%DATAKE) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL COPYIN_TYPED_VAR(VAR%GRID)
       CALL ACC_ATTACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEDES1


!****************** SUBROUTINE COPYIN_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELTOT needs to be created first, the same holds for FERWE and FERTOT
!$ACC ENTER DATA COPYIN(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%FERWE, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL COPYIN_TYPED_VAR(VAR%WDES)
       CALL ACC_ATTACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVESPIN


!****************** SUBROUTINE COPYIN_WAVEFUNA *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEFUNA(VAR)
    TYPE(wavefuna) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%CW_RED, &
!$ACC&                  VAR%CPROJ_RED) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL COPYIN_TYPED_VAR(VAR%WDES1)
       CALL ACC_ATTACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEFUNA


!****************** SUBROUTINE COPYIN_WAVEFUN1 *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEFUN1(VAR)
    TYPE(wavefun1) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$ACC ENTER DATA COPYIN(VAR%CPTWFP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CPROJ)) THEN
!$ACC ENTER DATA COPYIN(VAR%CPROJ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CR)) THEN
!$ACC ENTER DATA COPYIN(VAR%CR) ASYNC(ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL COPYIN_TYPED_VAR(VAR%WDES1)
       CALL ACC_ATTACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEFUN1


!****************** SUBROUTINE COPYIN_WAVEFUN1_ARRAY *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEFUN1_ARRAY(VAR)
    INTEGER  :: NB
    TYPE(wavefun1) :: VAR(:)
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR)
    DO NB=1,SIZE(VAR)
      ! native type pointers if associated
      IF (ASSOCIATED(VAR(NB)%CPTWFP)) THEN
!$ACC ENTER DATA COPYIN(VAR(NB)%CPTWFP)
      ENDIF
      IF (ASSOCIATED(VAR(NB)%CPROJ)) THEN
!$ACC ENTER DATA COPYIN(VAR(NB)%CPROJ)
      ENDIF
      IF (ASSOCIATED(VAR(NB)%CR)) THEN
!$ACC ENTER DATA COPYIN(VAR(NB)%CR)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR(NB)%WDES1)) THEN
         CALL COPYIN_TYPED_VAR(VAR(NB)%WDES1)
         CALL ACC_ATTACH_ASYNC(VAR(NB)%WDES1,ACC_ASYNC_Q)
      ENDIF
    ENDDO
  END SUBROUTINE COPYIN_WAVEFUN1_ARRAY


!****************** SUBROUTINE COPYIN_POTCARS **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_POTCARS(VAR)
    TYPE(potcar) :: VAR(:)
    INTEGER :: NT
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
    DO NT = 1, SIZE(VAR)
!$ACC ENTER DATA COPYIN(VAR(NT)%PSP, &
!$ACC&                  VAR(NT)%PSPCOR, &
!$ACC&                  VAR(NT)%PSPRHO, &
!$ACC&                  VAR(NT)%LPS) ASYNC(ACC_ASYNC_Q)
    ENDDO
  END SUBROUTINE COPYIN_POTCARS


!****************** SUBROUTINE COPYIN_KPOINTS_FULL *********************
!
!***********************************************************************

  SUBROUTINE COPYIN_KPOINTS_FULL(VAR)
    TYPE(skpoints_full) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%NG_INDEX)) THEN
!$ACC ENTER DATA COPYIN(VAR%NG_INDEX) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%PHASE)) THEN
!$ACC ENTER DATA COPYIN(VAR%PHASE) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%LINV)) THEN
!$ACC ENTER DATA COPYIN(VAR%LINV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%NEQUIV)) THEN
!$ACC ENTER DATA COPYIN(VAR%NEQUIV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%SPINFLIP)) THEN
!$ACC ENTER DATA COPYIN(VAR%SPINFLIP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ROTMAP)) THEN
!$ACC ENTER DATA COPYIN(VAR%ROTMAP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ISYMOP)) THEN
!$ACC ENTER DATA COPYIN(VAR%ISYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%MAP_INTO_BZ)) THEN
!$ACC ENTER DATA COPYIN(VAR%MAP_INTO_BZ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%VKPT)) THEN
!$ACC ENTER DATA COPYIN(VAR%VKPT) ASYNC(ACC_ASYNC_Q)
    ENDIF
#ifndef oldsym
    IF (ASSOCIATED(VAR%RSSYMOP)) THEN
!$ACC ENTER DATA COPYIN(VAR%RSSYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#endif
  END SUBROUTINE COPYIN_KPOINTS_FULL


!****************** SUBROUTINE COPYIN_GREENSFDES ***********************
!
!***********************************************************************

  SUBROUTINE COPYIN_GREENSFDES(VAR)
    TYPE(greensfdes) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NIONS_CORE, &
!$ACC&                  VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%NT_GLOBAL, &
!$ACC&                  VAR%NI_GLOBAL, &
!$ACC&                  VAR%POSION, &
!$ACC&                  VAR%NPRO_LMMAX, &
!$ACC&                  VAR%NPRO_LMBASE, &
!$ACC&                  VAR%NPRO_POS, &
!$ACC&                  VAR%NPRO_ENTRIES, &
!$ACC&                  VAR%NLM_LMMAX, &
!$ACC&                  VAR%NLM_POS, &
!$ACC&                  VAR%NLM_ENTRIES, &
!$ACC&                  VAR%COMM, &
!$ACC&                  VAR%LUSEINV, &
!$ACC&                  VAR%NGVECTOR_INV, &
!$ACC&                  VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$ACC&                  VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$ACC&                  VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$ACC&                  VAR%MAP_TO_FULL) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_GREENSFDES


!****************** SUBROUTINE COPYIN_XC *******************************
!
!***********************************************************************

  SUBROUTINE COPYIN_XC(VAR)
    TYPE(xc_info) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ALLOW .OR. IS_ON_DEVICE(VAR,1)) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%ID, &
!$ACC&                  VAR%FAMILY, &
!$ACC&                  VAR%KIND, &
!$ACC&                  VAR%ID_LDA_TABLE, &
!$ACC&                  VAR%NPARAM, &
!$ACC&                  VAR%COEFF, &
!$ACC&                  VAR%PARAM, &
!$ACC&                  VAR%NAME) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_XC


!****************** SUBROUTINE CREATE_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE CREATE_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC ENTER DATA CREATE(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELTOT needs to be created first, the same holds for FERWE and FERTOT
!$ACC ENTER DATA CREATE(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%FERWE, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL COPYIN_TYPED_VAR(VAR%WDES)
       CALL ACC_ATTACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE CREATE_WAVESPIN


!****************** SUBROUTINE DELETE_NONL_STRUCT **********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONL_STRUCT(VAR)
    TYPE(nonl_struct) :: VAR
    ! early exit if possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%POSION, &
!$ACC&                 VAR%VKPT_SHIFT, &
!$ACC&                 VAR%CQFAK, &
!$ACC&                 VAR%QPROJ, &
!$ACC&                 VAR%CREXP) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONL_STRUCT


!****************** SUBROUTINE DELETE_NONLR_STRUCT *********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONLR_STRUCT(VAR)
    TYPE(nonlr_struct) :: VAR
    ! logical variables
    INTEGER :: I
    ! early exit if possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers
    DO I=1,SIZE(VAR%BETA)
       CALL DELETE_TYPED_VAR(VAR%BETA(I))
#if PGI_BEFORE_24_5
       CALL ACC_DETACH_ASYNC(VAR%BETA(I),ACC_ASYNC_Q)
#endif
    ENDDO
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%LMAX, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%CHANNELS, &
!$ACC&                 VAR%PSRMAX, &
!$ACC&                 VAR%RSMOOTH, &
!$ACC&                 VAR%POSION, &
!$ACC&                 VAR%VKPT_SHIFT, &
!$ACC&                 VAR%NLIMAX, &
!$ACC&                 VAR%NLI, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NLIBASE, &
!$ACC&                 VAR%RPROJ, &
!$ACC&                 VAR%CRREXP) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONLR_STRUCT


!****************** SUBROUTINE DELETE_NONLR_PROJ ***********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONLR_PROJ(VAR)
    TYPE(nonlr_proj) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%PSPRNL, &
!$ACC&                 VAR%LPS) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONLR_PROJ


!****************** SUBROUTINE DELETE_LAYOUT ***************************
!
!***********************************************************************

  SUBROUTINE DELETE_LAYOUT(VAR)
    TYPE(layout) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%I2, &
!$ACC&                 VAR%I3, &
!$ACC&                 VAR%INDEX) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_LAYOUT


!****************** SUBROUTINE DELETE_GRID_MAP *************************
!
!***********************************************************************

  SUBROUTINE DELETE_GRID_MAP(VAR)
    TYPE(grid_map) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%N, &
!$ACC&                 VAR%PTR, &
!$ACC&                 VAR%TBL, &
!$ACC&                 VAR%NI, &
!$ACC&                 VAR%PTRI, &
!$ACC&                 VAR%TBLI) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_GRID_MAP


!****************** SUBROUTINE DELETE_GRID_3D **************************
!
!***********************************************************************

  SUBROUTINE DELETE_GRID_3D(VAR)
    TYPE(grid_3d) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers
    CALL DELETE_TYPED_VAR(VAR%RC)
    IF (.NOT.IS_ON_DEVICE(VAR%RC,1)) CALL ACC_DETACH_ASYNC(VAR%RC,ACC_ASYNC_Q)

    CALL DELETE_TYPED_VAR(VAR%IN)
    IF (.NOT.IS_ON_DEVICE(VAR%IN,1)) CALL ACC_DETACH_ASYNC(VAR%IN,ACC_ASYNC_Q)

    CALL DELETE_TYPED_VAR(VAR%RL)
    IF (.NOT.IS_ON_DEVICE(VAR%RL,1)) CALL ACC_DETACH_ASYNC(VAR%RL,ACC_ASYNC_Q)

    CALL DELETE_TYPED_VAR(VAR%RC_IN)
    IF (.NOT.IS_ON_DEVICE(VAR%RC_IN,1)) CALL ACC_DETACH_ASYNC(VAR%RC_IN,ACC_ASYNC_Q)

    CALL DELETE_TYPED_VAR(VAR%IN_RL)
    IF (.NOT.IS_ON_DEVICE(VAR%IN_RL,1)) CALL ACC_DETACH_ASYNC(VAR%IN_RL,ACC_ASYNC_Q)

    ! native type pointer if associated
    IF (ASSOCIATED(VAR%RL_FFT)) THEN
!!$ACC EXIT DATA DELETE(VAR%RL_FFT) ASYNC(ACC_ASYNC_Q)
       CALL DELETE_TYPED_VAR(VAR%RL_FFT)
       IF (.NOT.IS_ON_DEVICE(VAR%RL_FFT,1)) CALL ACC_DETACH_ASYNC(VAR%RL_FFT,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%LPCTX, &
!$ACC&                 VAR%LPCTY, &
!$ACC&                 VAR%LPCTZ, &
!$ACC&                 VAR%LPCTX_, &
!$ACC&                 VAR%LPCTY_, &
!$ACC&                 VAR%LPCTZ_, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%FFTWEIGHT, &
!$ACC&                 VAR%NINDPWCONJG, &
!$ACC&                 VAR%IND_IN_SPHERE, &
!$ACC&                 VAR%COMM, &
!$ACC&                 VAR%COMM_KIN, &
!$ACC&                 VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_GRID_3D


!****************** SUBROUTINE DELETE_TRANSIT **************************
!
!***********************************************************************

  SUBROUTINE DELETE_TRANSIT(VAR)
    TYPE(transit) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%IND1, &
!$ACC&                  VAR%INDCOL) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_TRANSIT


!****************** SUBROUTINE DELETE_WAVEDES **************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEDES(VAR)
    TYPE(wavedes) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL DELETE_TYPED_VAR(VAR%GRID)
       CALL ACC_DETACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NPLWKP, &
!$ACC&                 VAR%NGVECTOR, &
!$ACC&                 VAR%NGVECTOR_POS, &
!$ACC&                 VAR%NPLWKP_TOT, &
!$ACC&                 VAR%NB_TOTK, &
!$ACC&                 VAR%PL_INDEX, &
!$ACC&                 VAR%PL_COL, &
!$ACC&                 VAR%NPRO_POS, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%NT_GLOBAL, &
!$ACC&                 VAR%VKPT, &
!$ACC&                 VAR%WTKPT, &
!$ACC&                 VAR%NINDPW, &
!$ACC&                 VAR%LUSEINV, &
!$ACC&                 VAR%NINDPW_INV, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%MAP_TO_FULL, &
!$ACC&                 VAR%IGX, &
!$ACC&                 VAR%IGY, &
!$ACC&                 VAR%IGZ, &
!$ACC&                 VAR%DATAKE, &
!$ACC&                 VAR%AT_GAMMA) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEDES

  SUBROUTINE DELETE_KERNEL_TRUNCATION(VAR)
    TYPE(kernel_truncation) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%FINE_GRID)) THEN
       CALL DELETE_TYPED_VAR(VAR%FINE_GRID)
       CALL ACC_DETACH_ASYNC(VAR%FINE_GRID,ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COARSE_GRID)) THEN
       CALL DELETE_TYPED_VAR(VAR%COARSE_GRID)
       CALL ACC_DETACH_ASYNC(VAR%COARSE_GRID,ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%TRANSLATION_TABLE)) THEN
       CALL DELETE_TYPED_VAR(VAR%TRANSLATION_TABLE)
       CALL ACC_DETACH_ASYNC(VAR%TRANSLATION_TABLE,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE

  SUBROUTINE DELETE_COULOMB_POTENTIAL(VAR)
    TYPE(coulomb_potential) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
    !
    CALL DELETE_TYPED_VAR(VAR%KERNEL_TRUNCATE)
  END SUBROUTINE


!****************** SUBROUTINE DELETE_WAVEDES1 *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEDES1(VAR)
    TYPE(wavedes1) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL DELETE_TYPED_VAR(VAR%GRID)
       CALL ACC_DETACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%PL_INDEX, &
!$ACC&                 VAR%PL_COL, &
!$ACC&                 VAR%NPRO_POS, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%NT_GLOBAL, &
!$ACC&                 VAR%NINDPW, &
!$ACC&                 VAR%NINDPW_INV, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%MAP_TO_FULL, &
!$ACC&                 VAR%IGX, &
!$ACC&                 VAR%IGY, &
!$ACC&                 VAR%IGZ, &
!$ACC&                 VAR%VKPT, &
!$ACC&                 VAR%DATAKE) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEDES1


!****************** SUBROUTINE DELETE_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL DELETE_TYPED_VAR(VAR%WDES)
       CALL ACC_DETACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%CPTWFP) IF(IS_ON_DEVICE(VAR%CPTWFP,1)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(VAR%CPROJ) IF(IS_ON_DEVICE(VAR%CPROJ,1)) ASYNC(ACC_ASYNC_Q)
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELEN needs to be freed first, the same holds for FERWE and FERTOT.
!$ACC EXIT DATA DELETE(VAR%FERWE, &
!$ACC&                 VAR%FERTOT, &
!$ACC&                 VAR%AUXTOT, &
!$ACC&                 VAR%CELEN, &
!$ACC&                 VAR%CELTOT) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVESPIN


!****************** SUBROUTINE DELETE_WAVEFUNA *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEFUNA(VAR)
    TYPE(wavefuna) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL DELETE_TYPED_VAR(VAR%WDES1)
       CALL ACC_DETACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%CW_RED, &
!$ACC&                  VAR%CPROJ_RED) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEFUNA


!****************** SUBROUTINE DELETE_WAVEFUN1 *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEFUN1(VAR)
    TYPE(wavefun1) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$ACC EXIT DATA DELETE(VAR%CPTWFP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CPROJ)) THEN
!$ACC EXIT DATA DELETE(VAR%CPROJ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CR)) THEN
!$ACC EXIT DATA DELETE(VAR%CR) ASYNC(ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL DELETE_TYPED_VAR(VAR%WDES1)
       CALL ACC_DETACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEFUN1


!****************** SUBROUTINE DELETE_WAVEFUN1_ARRAY *******************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEFUN1_ARRAY(VAR)
    INTEGER  :: NB
    TYPE(wavefun1) :: VAR(:)
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers if associated
    DO NB=1,SIZE(VAR)
      IF (ASSOCIATED(VAR(NB)%CPTWFP)) THEN
!$ACC EXIT DATA DELETE(VAR(NB)%CPTWFP) ASYNC(ACC_ASYNC_Q)
      ENDIF
      IF (ASSOCIATED(VAR(NB)%CPROJ)) THEN
!$ACC EXIT DATA DELETE(VAR(NB)%CPROJ) ASYNC(ACC_ASYNC_Q)
      ENDIF
      IF (ASSOCIATED(VAR(NB)%CR)) THEN
!$ACC EXIT DATA DELETE(VAR(NB)%CR) ASYNC(ACC_ASYNC_Q)
      ENDIF
      ! derived type pointers if associated
      IF (ASSOCIATED(VAR(NB)%WDES1)) THEN
         CALL DELETE_TYPED_VAR(VAR(NB)%WDES1)
         CALL ACC_DETACH_ASYNC(VAR(NB)%WDES1,ACC_ASYNC_Q)
      ENDIF
    ENDDO
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEFUN1_ARRAY


!****************** SUBROUTINE DELETE_POTCARS **************************
!
!***********************************************************************

  SUBROUTINE DELETE_POTCARS(VAR)
    TYPE(potcar) :: VAR(:)
    INTEGER :: NT
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
    ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
    DO NT = 1, SIZE(VAR)
!$ACC EXIT DATA DELETE(VAR(NT)%PSP, &
!$ACC&                 VAR(NT)%PSPCOR, &
!$ACC&                 VAR(NT)%PSPRHO, &
!$ACC&                 VAR(NT)%LPS) ASYNC(ACC_ASYNC_Q)
    ENDDO
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_POTCARS


!****************** SUBROUTINE DELETE_KPOINTS_FULL *********************
!
!***********************************************************************

  SUBROUTINE DELETE_KPOINTS_FULL(VAR)
    TYPE(skpoints_full) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%NG_INDEX).AND.IS_ON_DEVICE(VAR%NG_INDEX,1)) THEN
!$ACC EXIT DATA DELETE(VAR%NG_INDEX) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%PHASE).AND.IS_ON_DEVICE(VAR%PHASE,1)) THEN
!$ACC EXIT DATA DELETE(VAR%PHASE) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%LINV)) THEN
!$ACC EXIT DATA DELETE(VAR%LINV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%NEQUIV)) THEN
!$ACC EXIT DATA DELETE(VAR%NEQUIV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%SPINFLIP)) THEN
!$ACC EXIT DATA DELETE(VAR%SPINFLIP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ROTMAP)) THEN
!$ACC EXIT DATA DELETE(VAR%ROTMAP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ISYMOP)) THEN
!$ACC EXIT DATA DELETE(VAR%ISYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%MAP_INTO_BZ)) THEN
!$ACC EXIT DATA DELETE(VAR%MAP_INTO_BZ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%VKPT)) THEN
!$ACC EXIT DATA DELETE(VAR%VKPT) ASYNC(ACC_ASYNC_Q)
    ENDIF
#ifndef oldsym
    IF (ASSOCIATED(VAR%RSSYMOP)) THEN
!$ACC EXIT DATA DELETE(VAR%RSSYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#endif
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_KPOINTS_FULL


!****************** SUBROUTINE DELETE_GREENSFDES ***********************
!
!***********************************************************************

  SUBROUTINE DELETE_GREENSFDES(VAR)
    TYPE(greensfdes) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NIONS_CORE, &
!$ACC&                 VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%NT_GLOBAL, &
!$ACC&                 VAR%NI_GLOBAL, &
!$ACC&                 VAR%POSION, &
!$ACC&                 VAR%NPRO_LMMAX, &
!$ACC&                 VAR%NPRO_LMBASE, &
!$ACC&                 VAR%NPRO_POS, &
!$ACC&                 VAR%NPRO_ENTRIES, &
!$ACC&                 VAR%NLM_LMMAX, &
!$ACC&                 VAR%NLM_POS, &
!$ACC&                 VAR%NLM_ENTRIES, &
!$ACC&                 VAR%COMM, &
!$ACC&                 VAR%LUSEINV, &
!$ACC&                 VAR%NGVECTOR_INV, &
!$ACC&                 VAR%NRPLWV_ROW_DATA_POINTS_NK, &
!$ACC&                 VAR%NRPLWV_COL_DATA_POINTS_NK, &
!$ACC&                 VAR%NRPLWV_COL_MAX_DATA_POINTS_NK, &
!$ACC&                 VAR%MAP_TO_FULL) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_GREENSFDES


!****************** SUBROUTINE DELETE_XC *******************************
!
!***********************************************************************

  SUBROUTINE DELETE_XC(VAR)
    TYPE(xc_info) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ALLOW .OR. .NOT.IS_ON_DEVICE(VAR,1)) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%ID, &
!$ACC&                 VAR%FAMILY, &
!$ACC&                 VAR%KIND, &
!$ACC&                 VAR%ID_LDA_TABLE, &
!$ACC&                 VAR%NPARAM, &
!$ACC&                 VAR%COEFF, &
!$ACC&                 VAR%PARAM, &
!$ACC&                 VAR%NAME) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
    ! safety check
    IF (IS_ON_DEVICE(VAR,1)) &
       CALL vtutor%bug("DELETE_XC: VAR still present on device.", __FILE__, __LINE__)
  END SUBROUTINE DELETE_XC


!****************** SUBROUTINE COPYOUT_WAVESPIN ************************
!
!***********************************************************************

  SUBROUTINE COPYOUT_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.OFFLOAD_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL DELETE_TYPED_VAR(VAR%WDES)
       CALL ACC_DETACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELEN needs to be freed first, the same holds for FERWE and FERTOT.
!$ACC EXIT DATA DELETE( VAR%FERWE, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYOUT_WAVESPIN

  END MODULE mopenacc_deepcopy


!***********************************************************************
!
!> The infrastructure in this module allows subroutines that are executed
!> on the GPU to raise critical exceptions to be acted  upon returning
!> after returning CPU side.
!> It only supports isBug and isError types of exceptions and does not
!> support passing information (arguments) from inside the GPU kernel.
!> Furthermore only one exception can be logged: when ACC_VTUTOR_LOG
!> is called multiple times it will store the instance with the largest
!> IdAndSeverity = Id + IDMAX * SEVERITY (see ACC_VTUTOR_LOG).
!
!***********************************************************************

  MODULE mopenacc_tutor
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ACC_VTUTOR_LOG, ACC_VTUTOR_CHECK

  INTEGER :: IdAndSeverity = -1
  INTEGER :: IDMAX = 1000
!$ACC DECLARE COPYIN(IdAndSeverity,IDMAX)

  CONTAINS

  SUBROUTINE ACC_VTUTOR_LOG(SEVERITY,ID)
!$ACC ROUTINE SEQ
    INTEGER :: SEVERITY,ID
!$ACC ATOMIC UPDATE
    IdAndSeverity = MAX(IdAndSeverity, SEVERITY*IDMAX + ID)
  END SUBROUTINE ACC_VTUTOR_LOG


  SUBROUTINE ACC_VTUTOR_CHECK
    USE moffload_struct_def, ONLY : OFFLOAD_ON
    USE tutor, ONLY : vtutor, DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR
    USE string, ONLY : str
    ! local variables
    CHARACTER(LEN=:), ALLOCATABLE :: FILENAME
    INTEGER :: ID,SEVERITY

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC UPDATE SELF(IdAndSeverity) IF(OFFLOAD_ON)
    IF (IdAndSeverity==-1) RETURN

    SEVERITY = IdAndSeverity/IDMAX
    ID = IdAndSeverity - SEVERITY*IDMAX

    SELECT CASE (ID)
    CASE(DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR)
       FILENAME = 'greens_real_space.F'
    CASE DEFAULT
       FILENAME = 'UnkownFile'
    END SELECT

    CALL vtutor%write(SEVERITY, ID, filename=FILENAME)

    ! currently this point will never be reached because we have
    ! only implemented fatal exceptions ... this might change in
    ! future though
    IdAndSeverity = -1
!$ACC UPDATE DEVICE(IdAndSeverity) IF(OFFLOAD_ON)
    RETURN
  END SUBROUTINE ACC_VTUTOR_CHECK

  END MODULE mopenacc_tutor
#endif // ACC_OFFLOAD

!***********************************************************************
!
! OpenACC related Doxygen documentation page.
!
!***********************************************************************

!> @page openacc OpenACC
!! @tableofcontents
!! @section general All changes
!! All datatypes, subroutines, and functions that have been changed
!! to port VASP to GPU under OpenACC.
!!
!! \li      subrot::add_gamma_from_file
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag
!! \li      subrot::eddiag_exact
!! \li    rmm_diis::eddrmm
!! \li       steep::edstep
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li     locproj::lprj_proall
!! \li            ::potlok
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
!! \li        wave::wvreal
!!
!! @section active Ported
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li      subrot::eddiag
!! \li    rmm_diis::eddrmm
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li        wave::wvreal
!!
!! @section inactive OpenACC deactivated
!! \li      subrot::add_gamma_from_file
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag_exact
!! \li       steep::edstep
!! \li     locproj::lprj_proall
!! \li         pot::potlok
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
