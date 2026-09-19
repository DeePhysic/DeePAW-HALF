#include "symbol.inc"
#ifdef NVCUDA
#ifdef ACC_OFFLOAD
!***********************************************************************
!
!> This module provides MPI related functions that INIT_ACC needs to
!> initialize the OpenACC infrastructure:
!>   - get the total number of devices
!>   - compute the number of ranks per device
!>   - check whether MPI is CUDA-aware.
!
!***********************************************************************
#ifdef MPI
  MODULE mdevice_info

  USE main_mpi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHECK_MPI, GET_NUM_DEVICES_GLOBAL, GET_NUM_RANKS_ON_DEVICE, GET_FREE_MEM_DEVICE

  CONTAINS

!****************** FUNCTION GET_NUM_DEVICES_GLOBAL ********************
!
!***********************************************************************

  FUNCTION GET_NUM_DEVICES_GLOBAL(MY_DEVICE_NUM)

    USE cudafor, ONLY : cudaGetDeviceProperties, cudaDeviceProp

    INTEGER :: GET_NUM_DEVICES_GLOBAL
    INTEGER :: MY_DEVICE_NUM

    TYPE (cudaDeviceProp) :: prop
    CHARACTER :: uuidchar(16)
    INTEGER :: uuidint(16)
    INTEGER, ALLOCATABLE :: uuidints(:,:)
    INTEGER, ALLOCATABLE :: uuidunique(:,:)
    INTEGER :: totalGPUs, I, J
    LOGICAL :: isUnique

    INTEGER :: ISTAT

    ISTAT = cudaGetDeviceProperties(prop, MY_DEVICE_NUM)
    READ(prop%uuid, '(16A1)') uuidchar
    uuidint = ICHAR(uuidchar(:))

    ALLOCATE(uuidints(16, COMM_intra_node_world%NCPU))

    uuidints = 0; uuidints(:, COMM_intra_node_world%NODE_ME) = uuidint(:)
    CALL M_sum_i(COMM_intra_node_world, uuidints(1,1), 16*COMM_intra_node_world%NCPU)

    ALLOCATE(uuidunique(16, COMM_intra_node_world%NCPU))

    uuidunique = 0 ; totalGPUs = 0
    DO I=1, COMM_intra_node_world%NCPU
       isUnique = .TRUE.
       DO J=1,totalGPUs
          IF ( ALL( uuidints(:,I) == uuidunique(:,J) ) ) THEN
             isUnique = .FALSE.
             EXIT
          ENDIF
       ENDDO
       IF ( isUnique ) THEN
          totalGPUs = totalGPUs + 1
          uuidunique(:,totalGPUs) = uuidints(:,I)
       ENDIF
    ENDDO

    DEALLOCATE(uuidints, uuidunique)

    CALL M_sum_i(COMM_inter_node_world, totalGPUs, 1)

    GET_NUM_DEVICES_GLOBAL = totalGPUs

  END FUNCTION GET_NUM_DEVICES_GLOBAL


!****************** FUNCTION GET_NUM_RANKS_ON_DEVICE *******************
!
!***********************************************************************

  FUNCTION GET_NUM_RANKS_ON_DEVICE(MY_DEVICE_NUM)

    USE cudafor, ONLY : cudaGetDeviceProperties, cudaDeviceProp

    INTEGER :: GET_NUM_RANKS_ON_DEVICE
    INTEGER :: MY_DEVICE_NUM

    TYPE (cudaDeviceProp) :: prop
    CHARACTER :: uuidchar(16)
    INTEGER :: uuidint(16)
    INTEGER, ALLOCATABLE :: uuidints(:,:)
    INTEGER :: numranks, I

    INTEGER :: ISTAT

    ISTAT = cudaGetDeviceProperties(prop, MY_DEVICE_NUM)
    READ(prop%uuid, '(16A1)') uuidchar
    uuidint = ICHAR(uuidchar(:))

    ALLOCATE(uuidints(16, COMM_intra_node_world%NCPU))

    uuidints = 0; uuidints(:, COMM_intra_node_world%NODE_ME) = uuidint(:)
    CALL M_sum_i(COMM_intra_node_world, uuidints(1,1), 16*COMM_intra_node_world%NCPU)

    numranks = 0

    DO I=1, COMM_intra_node_world%NCPU
       IF ( ALL( uuidints(:,I) == uuidint(:) ) ) numranks = numranks + 1
    ENDDO

    DEALLOCATE(uuidints)

    GET_NUM_RANKS_ON_DEVICE = numranks

  END FUNCTION GET_NUM_RANKS_ON_DEVICE


!****************** FUNCTION GET_FREE_MEM_DEVICE ***********************
!
!***********************************************************************

  FUNCTION GET_FREE_MEM_DEVICE(MFREE) RESULT(IERR)
    USE openacc, ONLY : ACC_DEVICE_KIND,ACC_GET_DEVICE_TYPE,ACC_DEVICE_NVIDIA,ACC_GET_DEVICE_NUM, &
                        ACC_GET_PROPERTY,ACC_PROPERTY_FREE_MEMORY
    USE iso_c_binding

    REAL(q) :: MFREE
    INTEGER :: IERR

    ! local variables
    INTEGER(KIND=ACC_DEVICE_KIND) :: DEVICE_TYPE
    INTEGER :: IDEVN

    INTEGER(KIND=c_size_t) :: ISIZE

    MFREE=0 ; IERR=0

    DEVICE_TYPE=ACC_GET_DEVICE_TYPE()
    IF (DEVICE_TYPE/=ACC_DEVICE_NVIDIA) THEN
       IERR=1 ; RETURN
    ENDIF

    IDEVN = ACC_GET_DEVICE_NUM(ACC_DEVICE_NVIDIA)
    ISIZE = ACC_GET_PROPERTY(IDEVN,ACC_DEVICE_NVIDIA,ACC_PROPERTY_FREE_MEMORY)

    MFREE = REAL(ISIZE,KIND=q)
  END FUNCTION GET_FREE_MEM_DEVICE


  SUBROUTINE CHECK_MPI()
    USE iso_c_binding
    USE tutor, ONLY : vtutor

#ifndef CRAY_MPICH
    INTERFACE
      INTEGER(c_int) FUNCTION MPIX_Query_cuda_support() BIND(C, name="MPIX_Query_cuda_support")
      END FUNCTION
    END INTERFACE
#endif

    LOGICAL :: CUDA_AWARE_SUPPORT

    CHARACTER(LEN=1) :: ENVVAR_VALUE
    INTEGER :: ENVVAR_STAT

#ifndef CRAY_MPICH
    ! This should tell us if MPI is CUDA-aware
    CUDA_AWARE_SUPPORT = MPIX_Query_cuda_support() == 1
#else
    CUDA_AWARE_SUPPORT = .TRUE.
#endif
    ! However, for OpenMPI some env variables can still deactivate it even though the previous
    ! check was positive
    CALL GET_ENVIRONMENT_VARIABLE("OMPI_MCA_mpi_cuda_support", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
    IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
    CALL GET_ENVIRONMENT_VARIABLE("OMPI_MCA_opal_cuda_support", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
    IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
    ! Just in case we might be non-OpenMPI, and their MPIX_Query_cuda_support behaves similarly
    CALL GET_ENVIRONMENT_VARIABLE("MV2_USE_CUDA", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
    IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
    CALL GET_ENVIRONMENT_VARIABLE("MPICH_RDMA_ENABLED_CUDA", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
    IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
    CALL GET_ENVIRONMENT_VARIABLE("PMPI_GPU_AWARE", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
    IF (ENVVAR_STAT==0) CUDA_AWARE_SUPPORT =(ENVVAR_VALUE == '1')

    IF (.NOT.CUDA_AWARE_SUPPORT)  CALL vtutor%error("CHECK_MPI: your MPI is not CUDA-aware, STOPPING")

    RETURN
  END SUBROUTINE CHECK_MPI

  END MODULE mdevice_info
#endif // MPI


!***********************************************************************
!
! This module contains the infrastructure to create and destroy
! cuFFT plans for grid_3d structures.
!
!***********************************************************************
#define use_cufft_workarea_minimal

  MODULE mcufft
  USE prec
  USE mgrid_struct_def, ONLY : grid_3d

  USE cufft
  USE openacc, ONLY : ACC_ASYNC_SYNC,ACC_ASYNC_NOVAL,acc_get_cuda_stream
  USE moffload_struct_def, ONLY : OFFLOAD_ALLOW,ACC_ASYNC_ASYNC,ACC_ASYNC_Q

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  USE iso_c_binding

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FFT_OFFLOAD_MAKEPLANS,FFT_OFFLOAD_MAKEPLAN,FFT_OFFLOAD_GETPLAN,FFT_OFFLOAD_DESTROYPLAN
  PUBLIC :: fft_offload_plan

  PUBLIC :: FFTBAS_PLAN_OFFLOAD, FFTBRC_PLAN_OFFLOAD
  PUBLIC :: FFTBAS_PLAN_MU_OFFLOAD, FFTBRC_PLAN_MU_OFFLOAD
  PUBLIC :: FFTBAS_PLAN_MPI_OFFLOAD, FFTBRC_PLAN_MPI_OFFLOAD

  TYPE fft_offload_plan
    TYPE(grid_3d), POINTER :: GRID
    INTEGER :: stream
    INTEGER :: nbatch
    INTEGER :: nstride
    INTEGER :: FFT_3D_Z2Z,FFT_3D_D2Z,FFT_3D_Z2D
    INTEGER :: FFT_RC_Z2Z,FFT_IN_Z2Z
    INTEGER :: FFT_RL_Z2Z,FFT_RL_D2Z,FFT_RL_Z2D
    INTEGER :: CUFFT_WORKAREA_POLICY
    CHARACTER(LEN=:), ALLOCATABLE :: tag
#if PGI_BEFORE_XX_X
    LOGICAL :: LFORCE_ACC_WAIT
#endif
  END TYPE fft_offload_plan

  TYPE(fft_offload_plan), ALLOCATABLE :: plans(:)

  INTEGER :: createdPlans = 0

#if PGI_BEFORE_XX_X
  INTEGER :: CUFFT_VER = -1
#endif
#if PGI_BEFORE_XX_X
  INTEGER :: DEVICE_CACHE_SIZE = -1
#endif
#if PGI_BEFORE_XX_X  && defined(use_cufft_workarea_minimal)
  INTEGER :: NV_COMPUTE_CAPABILITY = -1
#endif

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE planAssign
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE planCompare
  END INTERFACE

  INTEGER(c_int), BIND(c,name="CUFFT_WORKAREA_MINIMAL") :: CUFFT_WORKAREA_MINIMAL

  INTERFACE
    INTEGER(c_int) FUNCTION cufftXtSetWorkAreaPolicy(plan, policy, workSize) BIND(c, name="cufftXtSetWorkAreaPolicy")
        IMPORT
        INTEGER(c_int), VALUE :: plan
        INTEGER(c_int), VALUE :: policy
        INTEGER(c_size_t)     :: workSize
    END FUNCTION cufftXtSetWorkAreaPolicy
  END INTERFACE

  CONTAINS

!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLANS *******************
!
!***********************************************************************

  SUBROUTINE FFT_OFFLOAD_MAKEPLANS(GRID,NSTREAMS)
    TYPE(grid_3d), INTENT(IN), TARGET :: GRID
    INTEGER, INTENT(IN), OPTIONAL :: NSTREAMS
    ! local variables
    INTEGER :: N,I

    PROFILING_START('fft_offload_makeplans')

    IF (.NOT.OFFLOAD_ALLOW) THEN
       PROFILING_STOP('fft_offload_makeplans')
       RETURN
    ENDIF

    ! Always create a plan for ACC_ASYNC_SYNC
    CALL FFT_OFFLOAD_MAKEPLAN(GRID,INT(ACC_ASYNC_SYNC))

    ! And for ACC_ASYNC_ASYNC
    CALL FFT_OFFLOAD_MAKEPLAN(GRID,ACC_ASYNC_ASYNC)

    ! Then for additional streams if requested
    N=1; IF (PRESENT(NSTREAMS)) N=NSTREAMS
    DO I=1,N
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,I)
    ENDDO

    ! And last but not least for ACC_ASYNC_NOVAL
    CALL FFT_OFFLOAD_MAKEPLAN(GRID,INT(ACC_ASYNC_NOVAL))

    PROFILING_STOP('fft_offload_makeplans')

  END SUBROUTINE FFT_OFFLOAD_MAKEPLANS


!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLAN ********************
!
!***********************************************************************
  SUBROUTINE FFT_OFFLOAD_MAKEPLAN(GRID,ISTREAM,NBATCH,NSTRIDE,TAG)
    TYPE(grid_3d), INTENT(IN), TARGET      :: GRID
    INTEGER, INTENT(IN)                    :: ISTREAM
    INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
    INTEGER, INTENT(IN), OPTIONAL          :: NSTRIDE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG

    ! local variables
    TYPE(fft_offload_plan) :: NEW

    PROFILING_START('fft_offload_makeplan')

    IF (.NOT.OFFLOAD_ALLOW) THEN
       PROFILING_STOP('fft_offload_makeplan')
       RETURN
    ENDIF

    NEW%GRID   =>GRID
    NEW%STREAM = ISTREAM

    NEW%NBATCH = 1 ; IF (PRESENT(NBATCH)) NEW%NBATCH = MAX(NBATCH,1)
    NEW%NSTRIDE=-1 ; IF (PRESENT(NSTRIDE)) NEW%NSTRIDE = NSTRIDE

    ! when we already have a plan for this GRID+ISTREAM+NBATCH we may exit here
    IF (planFind(NEW)/=0) THEN
       PROFILING_STOP('fft_offload_makeplan')
       RETURN
    ENDIF

    NEW%TAG = "" ; IF (PRESENT(TAG)) NEW%TAG = TAG

    CALL MAKEPLAN(NEW)

    CALL planAdd(NEW)

    PROFILING_STOP('fft_offload_makeplan')

  END SUBROUTINE FFT_OFFLOAD_MAKEPLAN


  SUBROUTINE MAKEPLAN(PLAN)
    IMPLICIT NONE
    TYPE(fft_offload_plan) :: PLAN
    ! local variables
    INTEGER :: NG(3),BATCHES(3),LD,IERR

#ifdef use_cufft_workarea_minimal
    INTEGER(KIND=c_size_t) :: ISIZE
#if PGI_BEFORE_XX_X
    INTEGER :: IDUMMY(1)
#endif
#endif

    IERR=0

#if PGI_BEFORE_XX_X
    IF (CUFFT_VER == -1 ) IERR = cufftGetVersion(CUFFT_VER)
#endif

    NG=PLAN%GRID%NGPTAR

#ifdef use_cufft_workarea_minimal
    ! per default we try to use the minimal workarea policy
    PLAN%CUFFT_WORKAREA_POLICY=CUFFT_WORKAREA_MINIMAL
    ! for R2C FFTs this is not possible
    IF (PLAN%GRID%REAL2CPLX) PLAN%CUFFT_WORKAREA_POLICY=CUFFT_WORKAREA_MINIMAL+1
#if PGI_BEFORE_XX_X
! can be removed once NV HPC SDK no longer ships CUDA 12.0
    IF (NV_COMPUTE_CAPABILITY==-1) NV_COMPUTE_CAPABILITY = GET_DEVICE_COMPUTE_CAPABILITY()
#ifdef __aarch64__
    IF (ANY(NG(2:3)==6) .AND. CUFFT_VER < 11000 .AND. CUFFT_VER > 10502 &
                        .AND. NV_COMPUTE_CAPABILITY >= 70 .AND. NV_COMPUTE_CAPABILITY <= 80) THEN
#else
    IF (ANY(NG(2:3)==6) .AND. CUFFT_VER < 10900 .AND. CUFFT_VER > 10502 &
                        .AND. NV_COMPUTE_CAPABILITY >= 70 .AND. NV_COMPUTE_CAPABILITY <= 80) THEN
#endif
       PLAN%CUFFT_WORKAREA_POLICY=CUFFT_WORKAREA_MINIMAL+1
    ENDIF
#endif
#else
    PLAN%CUFFT_WORKAREA_POLICY=CUFFT_WORKAREA_MINIMAL+1
#endif

    IF (PLAN%GRID%RL%NFAST==1 .AND. PLAN%GRID%RL_FFT%NFAST==1) THEN
       ! serial fft
       IF (PLAN%NBATCH==1) THEN
          ! create single fft plans
          IF (PLAN%GRID%REAL2CPLX) THEN
             ! 3d real-to-complex
             IERR = IERR + cufftPlan3d(PLAN%FFT_3D_D2Z,NG(3),NG(2),NG(1),CUFFT_D2Z)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_D2Z,acc_get_cuda_stream(PLAN%STREAM))
             ! 3d complex-to-real
             IERR = IERR + cufftPlan3d(PLAN%FFT_3D_Z2D,NG(3),NG(2),NG(1),CUFFT_Z2D)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2D,acc_get_cuda_stream(PLAN%STREAM))
          ELSE
#ifdef use_cufft_workarea_minimal
             IF (PLAN%CUFFT_WORKAREA_POLICY==CUFFT_WORKAREA_MINIMAL) THEN
                IERR = IERR + cufftCreate(PLAN%FFT_3D_Z2Z)
                ! write(*,*) 'ierror cufftCreate :',ierr

                IERR = IERR + cufftSetAutoAllocation(PLAN%FFT_3D_Z2Z,0)
                ! write(*,*) 'ierror SetAutoAlloc:',ierr

                IERR = IERR + cufftPlan3d(PLAN%FFT_3D_Z2Z,NG(3),NG(2),NG(1),CUFFT_Z2Z)

                IERR = IERR + cufftXtSetWorkAreaPolicy(PLAN%FFT_3D_Z2Z,CUFFT_WORKAREA_MINIMAL,ISIZE)
#if PGI_BEFORE_XX_X
                IF (CUFFT_VER < 10500) THEN
!$ACC ENTER DATA CREATE(IDUMMY)
!$ACC HOST_DATA USE_DEVICE(IDUMMY)
                   IERR = IERR + cufftSetWorkArea(PLAN%FFT_3D_Z2Z,IDUMMY)
!$ACC END HOST_DATA
!$ACC EXIT DATA DELETE(IDUMMY)
                   ! write(*,*) 'ierror SetWorkArea:',ierr
                ENDIF
#endif
             ELSE
                IERR = IERR + cufftPlan3d(PLAN%FFT_3D_Z2Z,NG(3),NG(2),NG(1),CUFFT_Z2Z)
             ENDIF

             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
#else
             ! 3d complex-to-complex
             IERR = IERR + cufftPlan3d(PLAN%FFT_3D_Z2Z,NG(3),NG(2),NG(1),CUFFT_Z2Z)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
#endif
          ENDIF
       ELSE
          ! size required for in-place FFTs (in complex words)
          LD=PLAN%GRID%MPLWV
          ! allow for custom strided array access
          IF (PLAN%NSTRIDE > 0) LD=PLAN%NSTRIDE

          ! create batched fft plans
          IF (PLAN%GRID%REAL2CPLX) THEN
             ! batched 3d real-to-complex
             IERR = IERR + cufftPlanMany(PLAN%FFT_3D_D2Z,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2), NG(1)+2   ],1,2*LD, &
                                         [NG(3),NG(2),(NG(1)+2)/2],1,  LD, &
                                         CUFFT_D2Z,PLAN%NBATCH)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_D2Z,acc_get_cuda_stream(PLAN%STREAM))
             ! batched 3d complex-to-real
             IERR = IERR + cufftPlanMany(PLAN%FFT_3D_Z2D,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2),(NG(1)+2)/2],1,  LD, &
                                         [NG(3),NG(2), NG(1)+2   ],1,2*LD, &
                                         CUFFT_Z2D,PLAN%NBATCH)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2D,acc_get_cuda_stream(PLAN%STREAM))
          ELSE
             ! batched 3d complex-to-complex
#ifdef use_cufft_workarea_minimal
             IF (PLAN%CUFFT_WORKAREA_POLICY==CUFFT_WORKAREA_MINIMAL) THEN
                IERR = IERR + cufftCreate(PLAN%FFT_3D_Z2Z)
                ! write(*,*) 'ierror cufftCreate :',ierr

                IERR = IERR + cufftSetAutoAllocation(PLAN%FFT_3D_Z2Z,0)
                ! write(*,*) 'ierror SetAutoAlloc:',ierr

                IERR = IERR + cufftMakePlanMany(PLAN%FFT_3D_Z2Z,3,[NG(3),NG(2),NG(1)], &
                                                [NG(3),NG(2),NG(1)],1,LD, &
                                                [NG(3),NG(2),NG(1)],1,LD, &
                                                CUFFT_Z2Z,PLAN%NBATCH,ISIZE)
                ! write(*,*) 'ierror MakePlanMany:',ierr,isize

                IERR = IERR + cufftXtSetWorkAreaPolicy(PLAN%FFT_3D_Z2Z,CUFFT_WORKAREA_MINIMAL,ISIZE)
                ! write(*,*) 'ierror SetWorkAreaPolicy:',ierr,cufft_workarea_minimal
#if PGI_BEFORE_XX_X
                IF (CUFFT_VER < 10500) THEN
!$ACC ENTER DATA CREATE(IDUMMY)
!$ACC HOST_DATA USE_DEVICE(IDUMMY)
                   IERR = IERR + cufftSetWorkArea(PLAN%FFT_3D_Z2Z,IDUMMY)
!$ACC END HOST_DATA
!$ACC EXIT DATA DELETE(IDUMMY)
                   ! write(*,*) 'ierror SetWorkArea:',ierr
                ENDIF
#endif
             ELSE
                IERR = IERR + cufftPlanMany(PLAN%FFT_3D_Z2Z,3,[NG(3),NG(2),NG(1)], &
                                            [NG(3),NG(2),NG(1)],1,LD,[NG(3),NG(2),NG(1)],1,LD,CUFFT_Z2Z,PLAN%NBATCH)
             ENDIF

             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
             ! write(*,*) 'ierror SetStream:',ierr
#else
             IERR = IERR + cufftPlanMany(PLAN%FFT_3D_Z2Z,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2),NG(1)],1,LD,[NG(3),NG(2),NG(1)],1,LD,CUFFT_Z2Z,PLAN%NBATCH)
             IERR = IERR + cufftSetStream(PLAN%FFT_3D_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
#endif
          ENDIF
       ENDIF
    ELSE
       IF (PLAN%NBATCH/=1) THEN
          CALL vtutor%bug("MAKEPLAN: batched parallel ffts not supported yet",__FILE__,__LINE__)
       ENDIF

       ! parallel fft
       BATCHES(1)=PLAN%GRID%RC%NCOL
       BATCHES(2)=PLAN%GRID%IN%NCOL
       BATCHES(3)=PLAN%GRID%RL_FFT%NCOL

       IF (BATCHES(1)>0) THEN
          ! fft on data with reciprocal space layout
          IERR = IERR + cufftPlanMany(PLAN%FFT_RC_Z2Z,1,NG(1), &
                                      NG(1),1,NG(1),NG(1),1,NG(1),CUFFT_Z2Z,BATCHES(1))
          IERR = IERR + cufftSetStream(PLAN%FFT_RC_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
       ENDIF
       IF (BATCHES(2)>0) THEN
          ! fft on data with intermediate layout
          IERR = IERR + cufftPlanMany(PLAN%FFT_IN_Z2Z,1,NG(2), &
                                      NG(2),BATCHES(2),1,NG(2),BATCHES(2),1,CUFFT_Z2Z,BATCHES(2))
          IERR = IERR + cufftSetStream(PLAN%FFT_IN_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
       ENDIF
       IF (BATCHES(3)>0) THEN
          ! fft on data with real space layout
          IF (PLAN%GRID%REAL2CPLX) THEN
             ! real-to-complex
             IERR = IERR + cufftPlanMany(PLAN%FFT_RL_D2Z,1,NG(3), &
                                         NG(3),1,NG(3)+2,NG(3),1,(NG(3)+2)/2,CUFFT_D2Z,BATCHES(3))
             IERR = IERR + cufftSetStream(PLAN%FFT_RL_D2Z,acc_get_cuda_stream(PLAN%STREAM))
             ! complex-to-real
             IERR = IERR + cufftPlanMany(PLAN%FFT_RL_Z2D,1,NG(3), &
                                         NG(3),1,(NG(3)+2)/2,NG(3),1,NG(3)+2,CUFFT_Z2D,BATCHES(3))
             IERR = IERR + cufftSetStream(PLAN%FFT_RL_Z2D,acc_get_cuda_stream(PLAN%STREAM))
          ELSE
             ! complex-to-complex
             IERR = IERR + cufftPlanMany(PLAN%FFT_RL_Z2Z,1,NG(3), &
                                         NG(3),1,NG(3),NG(3),1,NG(3),CUFFT_Z2Z,BATCHES(3))
             IERR = IERR + cufftSetStream(PLAN%FFT_RL_Z2Z,acc_get_cuda_stream(PLAN%STREAM))
          ENDIF
       ENDIF
    ENDIF

#if PGI_BEFORE_XX_X
    CALL L2CACHE_OPT_WAR(PLAN)
#endif

    ! In case something went wrong along the way ...
    IF (IERR/=0) THEN
       CALL vtutor%error("MAKEPLAN: could not create plan")
    ENDIF
  END SUBROUTINE MAKEPLAN


!****************** SUBROUTINE FFT_OFFLOAD_DESTROYPLAN *****************
!
!***********************************************************************

  SUBROUTINE FFT_OFFLOAD_DESTROYPLAN(GRID,ISTREAM,NBATCH,TAG)
    TYPE(grid_3d), INTENT(IN), POINTER     :: GRID
    INTEGER, INTENT(IN), OPTIONAL          :: ISTREAM
    INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG
    ! local variables
    TYPE(fft_offload_plan), ALLOCATABLE :: tempPlans(:)
    INTEGER :: i,j
    LOGICAL :: keepPlan

    IF (.NOT.OFFLOAD_ALLOW) RETURN

    j=0

    ALLOCATE(tempPlans(createdPlans))

    DO i=1,createdPlans
       keepPlan = .FALSE.

       IF (.NOT.ASSOCIATED(GRID,plans(i)%GRID)) keepPLan = .TRUE.

       IF (PRESENT(ISTREAM)) THEN
          IF (plans(i)%stream/=ISTREAM) keepPlan = .TRUE.
       ENDIF

       IF (PRESENT(NBATCH)) THEN
          IF (plans(i)%NBATCH/=MAX(NBATCH,1)) keepPlan = .TRUE.
       ENDIF

       IF (PRESENT(TAG)) keepPlan = (plans(i)%TAG/=TAG)

       IF (keepPlan) THEN
          j=j+1; tempPlans(j)=plans(i)
       ELSE
          CALL DESTROYPLAN(plans(i))
       ENDIF
    ENDDO

    DEALLOCATE(plans)

    IF (j/=0) THEN
       ALLOCATE(plans(j))
       DO i=1,j
          plans(i)=tempPlans(i)
       ENDDO
    ENDIF

    DEALLOCATE(tempPlans)

    createdPlans=j

  END SUBROUTINE FFT_OFFLOAD_DESTROYPLAN


  SUBROUTINE DESTROYPLAN(PLAN)
    IMPLICIT NONE
    TYPE(fft_offload_plan) :: PLAN
    ! local variables
    INTEGER :: ierr
    ierr = 0
    IF (PLAN%GRID%RL%NFAST==1) THEN
       ! serial fft
       IF (PLAN%GRID%REAL2CPLX) THEN
          ierr = ierr + cufftDestroy(PLAN%FFT_3D_D2Z)
          ierr = ierr + cufftDestroy(PLAN%FFT_3D_Z2D)
       ELSE
          ierr = ierr + cufftDestroy(PLAN%FFT_3D_Z2Z)
       ENDIF
    ELSE
       IF (PLAN%GRID%RC%NCOL>0) &
          ierr = ierr + cufftDestroy(PLAN%FFT_RC_Z2Z)
       IF (PLAN%GRID%IN%NCOL>0) &
          ierr = ierr + cufftDestroy(PLAN%FFT_IN_Z2Z)
       IF (PLAN%GRID%RL_FFT%NCOL>0) THEN
          IF (PLAN%GRID%REAL2CPLX) THEN
             ierr = ierr + cufftDestroy(PLAN%FFT_RL_D2Z)
             ierr = ierr + cufftDestroy(PLAN%FFT_RL_Z2D)
          ELSE
             ierr = ierr + cufftDestroy(PLAN%FFT_RL_Z2Z)
          ENDIF
       ENDIF
    ENDIF
    IF (ierr/=0) THEN
       CALL vtutor%error("DESTROYPLAN: could not destroy plans")
    ENDIF
  END SUBROUTINE DESTROYPLAN


!****************** FUNCTION FFT_OFFLOAD_GETPLAN ***********************
!
!***********************************************************************

  FUNCTION FFT_OFFLOAD_GETPLAN(GRID,ISTREAM,NBATCH,RES) RESULT(ISTAT)
    TYPE(grid_3d), INTENT(IN), TARGET :: GRID
    INTEGER, INTENT(IN) :: ISTREAM
    INTEGER, INTENT(IN) :: NBATCH
    TYPE(fft_offload_plan), INTENT(OUT) :: RES
    INTEGER :: ISTAT
    ! local variables
    TYPE(fft_offload_plan) :: q
    INTEGER :: i

    q%GRID   =>GRID
    q%stream = ISTREAM
    q%nbatch = NBATCH

    ISTAT=0

    ! search matching plan
    i=planFind(q)

    IF (i/=0) THEN
       ! return matching plan
       RES=plans(i)
    ELSE
      ! no match found: return in error
      ISTAT=-1
    ENDIF
  END FUNCTION FFT_OFFLOAD_GETPLAN


!****************** SUBROUTINE PLANADD *********************************
!
!***********************************************************************

  SUBROUTINE planAdd(a)
    TYPE(fft_offload_plan), INTENT(IN) :: a
    ! local variables
    INTEGER :: i

    ! Early exit if plan already exists
    IF (planFind(a)/=0) RETURN

    ! Find storage space for a new plan
    i=planAlloc()

    ! Add plan
    plans(i)=a
  END SUBROUTINE planAdd


!****************** FUNCTION PLANALLOC *********************************
!
!***********************************************************************

  FUNCTION planAlloc
    INTEGER :: planAlloc
    ! local variables
    TYPE(fft_offload_plan), ALLOCATABLE :: tempPlans(:)
    INTEGER :: i,oldSize

    ! on first call the plans array is not allocated
    IF (.NOT.ALLOCATED(plans)) THEN
       ALLOCATE(plans(1))
       NULLIFY(plans(1)%grid)
       createdPlans=0

    ! when the plans array is full
    ELSEIF(createdPlans==SIZE(plans)) THEN
       oldSize=SIZE(plans)
       ALLOCATE(tempPlans(oldsize))
       DO i=1,oldSize
          tempPlans(i)=plans(i)
       ENDDO
       DEALLOCATE(plans)
       ALLOCATE(plans(oldsize+1))
       DO i=1,oldSize
          plans(i)=tempPlans(i)
       ENDDO
       DEALLOCATE(tempPlans)
       NULLIFY(plans(oldsize+1)%grid)

    ! error .. this may not happen
    ELSEIF(createdPlans>SIZE(plans)) THEN
       STOP
    ENDIF

    ! advance the number of created plans by one
    createdPlans=createdPlans+1

    ! and pass the index back
    planAlloc=createdPlans
  END FUNCTION planAlloc


!****************** FUNCTION PLANFIND **********************************
!
!***********************************************************************

  FUNCTION planFind(a)
    TYPE(fft_offload_plan), INTENT(IN) :: a
    INTEGER :: planFind
    ! local variables
    INTEGER :: i

    planFind=0
    DO i=1,createdPlans
       IF (plans(i)==a) THEN
          planFind=i; RETURN
       ENDIF
    ENDDO
  END FUNCTION planFind


!****************** SUBROUTINE PLANASSIGN ******************************
!
!***********************************************************************

  SUBROUTINE planAssign(a,b)
    TYPE(fft_offload_plan), INTENT(OUT) :: a
    TYPE(fft_offload_plan), INTENT(IN)  :: b

    a%GRID          =>b%GRID
    a%stream        = b%stream
    a%nbatch        = b%nbatch
    a%nstride       = b%nstride
    a%tag           = b%tag

    a%FFT_3D_Z2Z    = b%FFT_3D_Z2Z
    a%FFT_3D_D2Z    = b%FFT_3D_D2Z
    a%FFT_3D_Z2D    = b%FFT_3D_Z2D

    a%FFT_RC_Z2Z    = b%FFT_RC_Z2Z
    a%FFT_IN_Z2Z    = b%FFT_IN_Z2Z
    a%FFT_RL_Z2Z    = b%FFT_RL_Z2Z
    a%FFT_RL_D2Z    = b%FFT_RL_D2Z
    a%FFT_RL_Z2D    = b%FFT_RL_Z2D

#if PGI_BEFORE_XX_X
    a%LFORCE_ACC_WAIT = b%LFORCE_ACC_WAIT
#endif
  END SUBROUTINE planAssign


!****************** FUNCTION PLANCOMPARE *******************************
!
!***********************************************************************

  FUNCTION planCompare(a,b) RESULT(res)
    TYPE(fft_offload_plan), INTENT(IN) :: a,b
    LOGICAL :: res

    res=.FALSE.

    IF (.NOT.ASSOCIATED(a%GRID,b%GRID)) RETURN
    IF ( a%stream/=b%stream) RETURN
    IF ( a%nbatch/=b%nbatch) RETURN

    res=.TRUE.
  END FUNCTION planCompare

#if PGI_BEFORE_XX_X
  SUBROUTINE L2CACHE_OPT_WAR(PLAN)
    USE moffload_struct_def
    TYPE(fft_offload_plan) :: PLAN
    INTEGER :: IERR
    INTEGER :: dataSize

    PLAN%LFORCE_ACC_WAIT = .FALSE.

    IF ( PLAN%nbatch == 1 ) RETURN

    IF ( CUFFT_VER == -1 ) IERR = cufftGetVersion(CUFFT_VER)
    IF ( CUFFT_VER < 10300 .OR. CUFFT_VER > 10501 ) RETURN

    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(1)) ) RETURN
    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(2)) ) RETURN
    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(3)) ) RETURN

    dataSize = PRODUCT(PLAN%GRID%NGPTAR) * PLAN%nbatch * 16
    IF ( DEVICE_CACHE_SIZE == -1 ) DEVICE_CACHE_SIZE=GET_DEVICE_CACHE_SIZE()
    IF ( DEVICE_CACHE_SIZE > dataSize ) RETURN

    PLAN%LFORCE_ACC_WAIT = .TRUE.

    CONTAINS
      LOGICAL FUNCTION IS_DIM_AFFECTED(X)
        INTEGER :: X
        INTEGER, PARAMETER, DIMENSION(30) :: LIST=(/2, 4, 8, 16, 32, 64, 128, &
                                                    3, 9, 81, 243, 729, 2187, 6561, &
                                                    5, 25, 125, 625, 3125, &
                                                    6, 36, 216, 1296, 7776, &
                                                    7, 49, 343, 2401, &
                                                   11, 121/)
        IS_DIM_AFFECTED = ANY(LIST==X)
        RETURN
      END FUNCTION

      INTEGER FUNCTION GET_DEVICE_CACHE_SIZE()
        USE cudafor, only : cudaGetDeviceProperties, cudaDeviceProp
        USE openacc, only : ACC_DEVICE_NVIDIA, acc_get_device_num
        INTEGER :: IERR
        type (cudaDeviceProp) :: prop

        IERR =cudaGetDeviceProperties(prop, acc_get_device_num(ACC_DEVICE_NVIDIA))
        GET_DEVICE_CACHE_SIZE = prop%l2CacheSize
        IF ( prop%major >= 8 ) GET_DEVICE_CACHE_SIZE=GET_DEVICE_CACHE_SIZE/2
        RETURN
      END FUNCTION GET_DEVICE_CACHE_SIZE

  END SUBROUTINE L2CACHE_OPT_WAR
#endif

#if PGI_BEFORE_XX_X && defined(use_cufft_workarea_minimal)
  FUNCTION GET_DEVICE_COMPUTE_CAPABILITY()

    USE cudafor, only : cudaGetDeviceProperties, cudaDeviceProp
    USE openacc, only : ACC_DEVICE_NVIDIA, acc_get_device_num

    INTEGER :: GET_DEVICE_COMPUTE_CAPABILITY

    TYPE (cudaDeviceProp) :: prop

    INTEGER :: ISTAT

    ISTAT = cudaGetDeviceProperties(prop, acc_get_device_num(ACC_DEVICE_NVIDIA))
    GET_DEVICE_COMPUTE_CAPABILITY = prop%major * 10 + prop%minor

  END FUNCTION GET_DEVICE_COMPUTE_CAPABILITY
#endif

!************************* SUBROUTINE FFTBAS_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      TYPE (grid_3d):: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbas_plan_offload')

      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) THEN
         CALL vtutor%bug("FFTBAS_PLAN_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // str(INT(ACC_ASYNC_Q)), __FILE__, __LINE__)
      ENDIF
!$ACC HOST_DATA USE_DEVICE(C)
      IF (ISIGN<=0) THEN
         IERR=cufftExecZ2Z(PLAN%FFT_3D_Z2Z,C,C,CUFFT_FORWARD)
      ELSE
         IERR=cufftExecZ2Z(PLAN%FFT_3D_Z2Z,C,C,CUFFT_INVERSE)
      ENDIF
!$ACC END HOST_DATA

      PROFILING_STOP('fftbas_plan_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_OFFLOAD


!************************* SUBROUTINE FFTBRC_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      REAL(q), POINTER  :: A(:)
      INTEGER :: IERR

      PROFILING_START('fftbrc_plan_offload')

      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) THEN
         CALL vtutor%bug("FFTBRC_PLAN_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // str(INT(ACC_ASYNC_Q)), __FILE__, __LINE__)
      ENDIF
      CALL c_f_pointer(c_loc(C),A,[2*GRID%MPLWV])
!$ACC HOST_DATA USE_DEVICE(C,A)
      IF (ISIGN<=0) THEN
         IERR=cufftExecD2Z(PLAN%FFT_3D_D2Z,A,C)
      ELSE
         IERR=cufftExecZ2D(PLAN%FFT_3D_Z2D,C,A)
      ENDIF
!$ACC END HOST_DATA

      PROFILING_STOP('fftbrc_plan_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISN)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbas_plan_mu_offload')
!
! Find corresponding plan
!
      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,PLAN)

! Plan not found? Exit in error
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_MU_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // &
            str(INT(ACC_ASYNC_Q)) // " and batch size N= " // str(N),__FILE__,__LINE__)

!
! Execute plan
!
#if PGI_BEFORE_XX_X
!$ACC WAIT(ACC_ASYNC_Q) IF(PLAN%LFORCE_ACC_WAIT)
#endif
!$ACC HOST_DATA USE_DEVICE(C)
      IF (ISN.LE.0) THEN
         IERR=cufftExecZ2Z(PLAN%FFT_3D_Z2Z,C,C,CUFFT_FORWARD)
      ELSE
         IERR=cufftExecZ2Z(PLAN%FFT_3D_Z2Z,C,C,CUFFT_INVERSE)
      ENDIF
!$ACC END HOST_DATA

      PROFILING_STOP('fftbas_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISN)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      REAL(q), POINTER :: R(:)
      INTEGER :: IERR

      PROFILING_START('fftbrc_plan_mu_offload')
!
! Find corresponding plan
!
      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,PLAN)

! Plan not found? Exit in error
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_MU_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // &
            str(INT(ACC_ASYNC_Q)) // " and batch size N= " // str(N),__FILE__,__LINE__)

!
! Connect a real pointer to the data buffer as well
     CALL c_f_pointer(c_loc(C),R,[2*LDC*N])

!
! Execute plan
!

!$ACC HOST_DATA USE_DEVICE(C,R)
      IF (ISN.LE.0) THEN
         IERR=cufftExecD2Z(PLAN%FFT_3D_D2Z,R,C)
      ELSE
         IERR=cufftExecZ2D(PLAN%FFT_3D_Z2D,C,R)
      ENDIF
!$ACC END HOST_DATA

! cleanup
      NULLIFY(R)

      PROFILING_STOP('fftbrc_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISIGN

      ! local variables
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)

      COMPLEX(q), POINTER  :: C(:)
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbas_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))

!$ACC DATA CREATE(RCVBUF,SNDBUF) ASYNC(ACC_ASYNC_Q)
      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) THEN
         CALL vtutor%bug("FFTBAS_PLAN_MPI_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // &
            str(INT(ACC_ASYNC_Q)), __FILE__, __LINE__)
      ENDIF
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
      IF (ISIGN==1) THEN
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RC_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_IN_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RL_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
      ELSE
! transformation along third dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RL_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_IN_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RC_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
      ENDIF
!$ACC END DATA

      PROFILING_STOP('fftbas_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISIGN

      ! local variables
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)

      COMPLEX(q), POINTER  :: C(:)
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbrc_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))

!$ACC DATA CREATE(RCVBUF,SNDBUF) ASYNC(ACC_ASYNC_Q)
      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

      IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) THEN
         CALL vtutor%bug("FFTBAS_PLAN_MPI_OFFLOAD: no FFT plan found for ACC_ASYNC_Q= " // &
            str(INT(ACC_ASYNC_Q)), __FILE__, __LINE__)
      ENDIF
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
      IF (ISIGN==1) THEN
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RC_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_IN_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
!$ACC HOST_DATA USE_DEVICE(C,A)
         IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2D(PLAN%FFT_RL_Z2D,C,A)
!$ACC END HOST_DATA
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
      ELSE
! transformation along third dimension:
!$ACC HOST_DATA USE_DEVICE(C,A)
         IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecD2Z(PLAN%FFT_RL_D2Z,A,C)
!$ACC END HOST_DATA
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_IN_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
         IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(PLAN%FFT_RC_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
      ENDIF
!$ACC END DATA

      PROFILING_STOP('fftbrc_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD

  END MODULE mcufft


!***********************************************************************
!
! This module only holds a few variables related to the use of cublas
! and a subroutine to initialize the cublasHandle ACC_CUBLASHANDLE.
!
! The actual cublas calls are wrapped in ACC_* subroutines that live
! outside of mcublas (see immediately below the aforementioned module).
!
!***********************************************************************

  MODULE mcublas
  USE prec

  USE cublas
  USE openacc,  ONLY : acc_get_cuda_stream

  USE moffload_struct_def, ONLY : OFFLOAD_ALLOW,ACC_ASYNC_Q

  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  TYPE(cublasHandle) :: ACC_CUBLASHANDLE
!$OMP THREADPRIVATE(ACC_CUBLASHANDLE)

  INTEGER(4) :: CUBLAS_VERSION=0

  PRIVATE

  PUBLIC :: DCOPY_OFFLOAD, SCOPY_OFFLOAD, CCOPY_OFFLOAD, ZCOPY_OFFLOAD
  PUBLIC :: DAXPY_OFFLOAD, ZAXPY_OFFLOAD
  PUBLIC :: DSCAL_OFFLOAD, ZSCAL_OFFLOAD, ZDSCAL_OFFLOAD
  PUBLIC :: DDOT_OFFLOAD, ZDOTC_OFFLOAD, CDOTC_OFFLOAD
  PUBLIC :: DGEMV_OFFLOAD, ZGEMV_OFFLOAD
  PUBLIC :: DGEMM_OFFLOAD, ZGEMM_OFFLOAD, SGEMM_OFFLOAD, CGEMM_OFFLOAD
  PUBLIC :: DTRMM_OFFLOAD, ZTRMM_OFFLOAD
  PUBLIC :: DLACPY_OFFLOAD

  PUBLIC :: ACC_CUBLASHANDLE
  PUBLIC :: CUBLAS_INIT, CUBLAS_FINISH
  PUBLIC :: CHAR_TO_FILL, CHAR_TO_DIAG, CHAR_TO_OP, CHAR_TO_SIDE

  CONTAINS

  INTEGER FUNCTION CHAR_TO_OP(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('N','n')
        CHAR_TO_OP = CUBLAS_OP_N
      CASE ('T','t')
        CHAR_TO_OP = CUBLAS_OP_T
      CASE ('C','c')
        CHAR_TO_OP = CUBLAS_OP_C
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_OP: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_OP


  INTEGER FUNCTION CHAR_TO_FILL(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('U','u')
        CHAR_TO_FILL = CUBLAS_FILL_MODE_UPPER
      CASE ('L','l')
        CHAR_TO_FILL = CUBLAS_FILL_MODE_LOWER
      !CASE ('F','f')
      !  CHAR_TO_FILL = CUBLAS_FILL_MODE_FULL
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_FILL: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_FILL


  INTEGER FUNCTION CHAR_TO_SIDE(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('L','l')
        CHAR_TO_SIDE = CUBLAS_SIDE_LEFT
      CASE ('R','r')
        CHAR_TO_SIDE = CUBLAS_SIDE_RIGHT
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_SIDE: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_SIDE


  INTEGER FUNCTION CHAR_TO_DIAG(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('N','n')
        CHAR_TO_DIAG = CUBLAS_DIAG_NON_UNIT
      CASE ('U','u')
        CHAR_TO_DIAG = CUBLAS_DIAG_UNIT
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_DIAG: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_DIAG


  SUBROUTINE CUBLAS_INIT()
    INTEGER(4) :: IERR

    PROFILING_START('cublas_init')

    IF (.NOT.OFFLOAD_ALLOW) THEN
       PROFILING_STOP('cublas_init')
       RETURN
    ENDIF

!$OMP PARALLEL
    IERR = cublasCreate(ACC_CUBLASHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR setting up cublasHandle")
    ENDIF
!$OMP END PARALLEL

    CALL GET_CUBLAS_VERSION()

    PROFILING_STOP('cublas_init')

  END SUBROUTINE CUBLAS_INIT


  SUBROUTINE CUBLAS_FINISH()
    INTEGER(4) :: IERR

    IF (.NOT.OFFLOAD_ALLOW) RETURN

!$OMP PARALLEL
    IERR = cublasDestroy(ACC_CUBLASHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR releasing cublasHandle")
    ENDIF
!$OMP END PARALLEL
  END SUBROUTINE CUBLAS_FINISH

  SUBROUTINE GET_CUBLAS_VERSION()
    USE tutor, ONLY: vtutor
    INTEGER(4) :: IERR

!$OMP PARALLEL
!$OMP MASTER
    IERR = cublasGetVersion(ACC_CUBLASHANDLE,CUBLAS_VERSION)
    IF (IERR/=0) THEN
       CALL vtutor%bug("ERROR in cublasGetVersion", __FILE__, __LINE__)
    ENDIF
!$OMP END MASTER
!$OMP END PARALLEL
  END SUBROUTINE GET_CUBLAS_VERSION


  SUBROUTINE DCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
    REAL(q) :: DX(*),DY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(DX(1:N),DY(1:N)) ASYNC(ACC_ASYNC_Q)
       DY(1:N)=DX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDcopy_v2(ACC_CUBLASHANDLE,N,DX,INCX,DY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE DCOPY_OFFLOAD


  SUBROUTINE SCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
    REAL(qs) :: DX(*),DY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(DX(1:N),DY(1:N)) ASYNC(ACC_ASYNC_Q)
       DY(1:N)=DX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasScopy_v2(ACC_CUBLASHANDLE,N,DX,INCX,DY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE SCOPY_OFFLOAD


  SUBROUTINE CCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
    COMPLEX(qs) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(ZX(1:N),ZY(1:N)) ASYNC(ACC_ASYNC_Q)
       ZY(1:N)=ZX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasCcopy_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE CCOPY_OFFLOAD


  SUBROUTINE ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
    COMPLEX(q) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(ZX(1:N),ZY(1:N)) ASYNC(ACC_ASYNC_Q)
       ZY(1:N)=ZX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZcopy_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ZCOPY_OFFLOAD


  SUBROUTINE DAXPY_OFFLOAD(N,DA,DX,INCX,DY,INCY)
    REAL(q) :: DX(*),DY(*)
    REAL(q) :: DA
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDaxpy_v2(ACC_CUBLASHANDLE,N,DA,DX,INCX,DY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE DAXPY_OFFLOAD


  SUBROUTINE ZAXPY_OFFLOAD(N,ZA,ZX,INCX,ZY,INCY)
    COMPLEX(q) :: ZX(*),ZY(*)
    COMPLEX(q) :: ZA
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZaxpy_v2(ACC_CUBLASHANDLE,N,ZA,ZX,INCX,ZY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ZAXPY_OFFLOAD


  SUBROUTINE DSCAL_OFFLOAD(N,DA,DX,INCX)
    REAL(q) :: DX(*)
    REAL(q) :: DA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX)
    IERR = cublasDscal_v2(ACC_CUBLASHANDLE,N,DA,DX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE DSCAL_OFFLOAD


  SUBROUTINE ZSCAL_OFFLOAD(N,ZA,ZX,INCX)
    COMPLEX(q) :: ZX(*)
    COMPLEX(q) :: ZA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX)
    IERR = cublasZscal_v2(ACC_CUBLASHANDLE,N,ZA,ZX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE ZSCAL_OFFLOAD


  SUBROUTINE ZDSCAL_OFFLOAD(N,DA,ZX,INCX)
    COMPLEX(q) :: ZX(*)
    REAL(q) :: DA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX)
    IERR = cublasZdscal_v2(ACC_CUBLASHANDLE,N,DA,ZX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE ZDSCAL_OFFLOAD


  REAL(q) FUNCTION DDOT_OFFLOAD(N,DX,INCX,DY,INCY)
    REAL(q) :: DX(*),DY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDdot_v2(ACC_CUBLASHANDLE,N,DX,INCX,DY,INCY,DDOT_OFFLOAD)
!$ACC END HOST_DATA
  END FUNCTION DDOT_OFFLOAD


  COMPLEX(q) FUNCTION ZDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
    COMPLEX(q) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY,ZDOTC_OFFLOAD)
!$ACC END HOST_DATA
  END FUNCTION ZDOTC_OFFLOAD


  COMPLEX(qs) FUNCTION CDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
    COMPLEX(qs) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasCdotc_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY,CDOTC_OFFLOAD)
!$ACC END HOST_DATA
  END FUNCTION CDOTC_OFFLOAD


  SUBROUTINE DGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    REAL(q) :: A(LDA,*),X(*),Y(*)
    REAL(q) :: ALPHA,BETA
    INTEGER :: M,N,LDA,INCX,INCY
    CHARACTER(1) :: MODE
    ! local variables
    REAL(q), ALLOCATABLE :: YTMP(:)
    INTEGER :: IDX,IDY,I
    LOGICAL :: LGEMM

    INTEGER(4) :: IERR

    IF (MODE=='N'.OR.MODE=='n') THEN
       IDX=N*INCX
       IDY=M*INCY
    ELSE
       IDX=M*INCX
       IDY=N*INCY
    ENDIF

    IF ( CUBLAS_VERSION < 10100 ) THEN
! Use DGEMM or not? This should be made conditional upon a performance test.
    LGEMM=.FALSE.

    IF (LGEMM.AND.INCX==1.AND.INCY==1) THEN
       IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ENDIF
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (MODE=='N'.OR.MODE=='n') THEN
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,M,1,N,ALPHA,A,LDA,X,N,BETA,Y,IDY)
!$ACC END HOST_DATA
       ELSE
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),CUBLAS_OP_N,N,1,M,ALPHA,A,LDA,X,M,BETA,Y,IDY)
!$ACC END HOST_DATA
       ENDIF
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       RETURN
    ENDIF

    IF (MODE=='T'.OR.MODE=='t') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasDdot_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasDdot_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    ENDIF

    IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
    IERR = cublasDgemv_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DGEMV_OFFLOAD


  SUBROUTINE ZGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    COMPLEX(q) :: A(LDA,*),X(*),Y(*)
    COMPLEX(q) :: ALPHA,BETA
    INTEGER :: M,N,LDA,INCX,INCY
    CHARACTER(1) :: MODE
    ! local variables
    COMPLEX(q), ALLOCATABLE :: YTMP(:)
    INTEGER :: IDX,IDY,I
    LOGICAL :: LGEMM

    INTEGER(4) :: IERR

    IF (MODE=='N'.OR.MODE=='n') THEN
       IDX=N*INCX
       IDY=M*INCY
    ELSE
       IDX=M*INCX
       IDY=N*INCY
    ENDIF

    IF ( CUBLAS_VERSION < 10100 ) THEN
! Use ZGEMM or not? This should be made conditional upon a performance test.
    LGEMM=.FALSE.

    IF (LGEMM.AND.INCX==1.AND.INCY==1) THEN
       IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ENDIF
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (MODE=='N'.OR.MODE=='n') THEN
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,M,1,N,ALPHA,A,LDA,X,N,BETA,Y,IDY)
!$ACC END HOST_DATA
       ELSE
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),CUBLAS_OP_N,N,1,M,ALPHA,A,LDA,X,M,BETA,Y,IDY)
!$ACC END HOST_DATA
       ENDIF
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       RETURN
    ENDIF

    IF (MODE=='C'.OR.MODE=='c') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    IF (MODE=='T'.OR.MODE=='t') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotu_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotu_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    ENDIF

    IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
    IERR = cublasZgemv_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE ZGEMV_OFFLOAD


  SUBROUTINE DGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    REAL(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
    REAL(q) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE DGEMM_OFFLOAD


  SUBROUTINE SGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    REAL(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
    REAL(qs) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasSgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE SGEMM_OFFLOAD


  SUBROUTINE ZGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    COMPLEX(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
    COMPLEX(q) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE ZGEMM_OFFLOAD


  SUBROUTINE CGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    COMPLEX(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
    COMPLEX(qs) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasCgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE CGEMM_OFFLOAD


  SUBROUTINE DTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    REAL(q) :: A(LDA,*),B(LDB,*)
    REAL(q) :: ALPHA
    INTEGER :: LDA,LDB,M,N
    CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasDtrmm_v2(ACC_CUBLASHANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,A,LDA,B,LDB,B,LDB)
!$ACC END HOST_DATA
  END SUBROUTINE DTRMM_OFFLOAD


  SUBROUTINE ZTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    COMPLEX(q) :: A(LDA,*),B(LDB,*)
    COMPLEX(q) :: ALPHA
    INTEGER :: LDA,LDB,M,N
    CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO
    ! local variables
    INTEGER(4) :: IERR

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasZtrmm_v2(ACC_CUBLASHANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,A,LDA,B,LDB,B,LDB)
!$ACC END HOST_DATA
  END SUBROUTINE ZTRMM_OFFLOAD


  SUBROUTINE DLACPY_OFFLOAD(UPLO,M,N,A,LDA,B,LDB)
    REAL(q) :: A(LDA,*),B(LDB,*)
    INTEGER :: LDA,LDB,M,N
    CHARACTER(1) :: UPLO
    ! local variables
    INTEGER(4) :: IERR

    IF (UPLO/='F'.AND.UPLO/='f') THEN
       CALL vtutor%bug("DLACPY_OFFLOAD: only UPLO=F (full matrix) is supported",__FILE__,__LINE__)
    ENDIF
#ifndef use_kernels_construct
    IF (M==LDA.AND.LDB==LDA) THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B)
       IERR = cublasDcopy_v2(ACC_CUBLASHANDLE,M*N,A,1,B,1)
!$ACC END HOST_DATA
       RETURN
    ENDIF
#endif
!$ACC KERNELS PRESENT(A(1:M,1:N),B(1:M,1:N)) ASYNC(ACC_ASYNC_Q)
    B(1:M,1:N)=A(1:M,1:N)
!$ACC END KERNELS
  END SUBROUTINE DLACPY_OFFLOAD

  END MODULE mcublas


!***********************************************************************
!
!> This module holds everything related to the use of cusolverDn.
!
!***********************************************************************

  MODULE mcusolver
  USE prec

  USE cusolverDn, ONLY : cusolverDnHandle, &
                         cusolverDnCreate, &
                         cusolverDnDestroy, &
                         cusolverDnSetStream, &
                         cusolverDnDpotrf_bufferSize, &
                         cusolverDnZpotrf_bufferSize, &
                         cusolverDnDgetrf_bufferSize, &
                         cusolverDnZgetrf_bufferSize, &
                         cusolverDnDsygvd_bufferSize, &
                         cusolverDnZhegvd_bufferSize, &
                         cusolverDnDsyevd_bufferSize, &
                         cusolverDnZheevd_bufferSize, &
                         cusolverDnCheevd_bufferSize, &
                         cusolverDnSsyevd_bufferSize, &
                         cusolverDnDpotrf, &
                         cusolverDnZpotrf, &
                         cusolverDnDgetrf, &
                         cusolverDnDgetrs, &
                         cusolverDnZgetrf, &
                         cusolverDnZgetrs, &
                         cusolverDnDsyevd, &
                         cusolverDnSsyevd, &
                         cusolverDnZheevd, &
                         cusolverDnCheevd, &
                         CUSOLVER_EIG_MODE_NOVECTOR, &
                         CUSOLVER_EIG_MODE_VECTOR, &
                         CUSOLVER_EIG_TYPE_1, &
                         CUSOLVER_EIG_TYPE_2, &
                         CUSOLVER_EIG_TYPE_3

  USE cublas, ONLY :     CUBLAS_FILL_MODE_UPPER, &
                         CUBLAS_FILL_MODE_LOWER, &
                         CUBLAS_OP_N, &
                         CUBLAS_OP_T, &
                         CUBLAS_OP_C, &
                         CUBLAS_SIDE_LEFT, &
                         CUBLAS_DIAG_NON_UNIT, &
                         cublasSetStream, &
                         cublasDtrsm_v2, &
                         cublasZtrsm_v2

  USE openacc,  ONLY : acc_get_cuda_stream,c_devloc

  USE moffload_struct_def, ONLY : OFFLOAD_ALLOW,OFFLOAD_ON,ACC_ASYNC_Q

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  USE mcublas

  IMPLICIT NONE

  TYPE(cusolverDnHandle) :: ACC_CUSOLVERHANDLE
!$OMP THREADPRIVATE(ACC_CUSOLVERHANDLE)

  INTERFACE
  INTEGER(c_int) FUNCTION cusolverDnZhegvd(handle,itype,jobz,uplo,n,A,lda,B,ldb,W,work,lwork,devInfo) BIND(C, name="cusolverDnZhegvd")
    USE ISO_C_BINDING
    USE cusolverDn, ONLY : cusolverDnHandle
    IMPLICIT NONE
    TYPE(cusolverDnHandle),VALUE::handle
    INTEGER(c_int),VALUE::itype
    INTEGER(c_int),VALUE::jobz
    INTEGER(c_int),VALUE::uplo
    INTEGER(c_int),VALUE::n
    !dir$ ignore_tkr(trd) A
    DOUBLE COMPLEX, DIMENSION(*) :: A
    INTEGER(c_int),VALUE::lda
    !dir$ ignore_tkr(trd) B
    DOUBLE COMPLEX, DIMENSION(*) :: B
    INTEGER(c_int),VALUE::ldb
    !dir$ ignore_tkr(trd) W
    DOUBLE PRECISION, DIMENSION(*)::W
    !dir$ ignore_tkr(trd) work
    DOUBLE COMPLEX, DIMENSION(*) :: work
    INTEGER(c_int),VALUE::lwork
    !dir$ ignore_tkr(trd) devInfo
    INTEGER, DIMENSION(1)::devInfo
  END FUNCTION
  END INTERFACE

  INTERFACE
  INTEGER(c_int) FUNCTION cusolverDnDsygvd(handle,itype,jobz,uplo,n,A,lda,B,ldb,W,work,lwork,devInfo) BIND(C, name="cusolverDnDsygvd")
    USE ISO_C_BINDING
    USE cusolverDn, ONLY : cusolverDnHandle
    IMPLICIT NONE
    TYPE(cusolverDnHandle),VALUE::handle
    INTEGER(c_int),VALUE::itype
    INTEGER(c_int),VALUE::jobz
    INTEGER(c_int),VALUE::uplo
    INTEGER(c_int),VALUE::n
    !dir$ ignore_tkr(trd) A
    DOUBLE PRECISION, DIMENSION(*) :: A
    INTEGER(c_int),VALUE::lda
    !dir$ ignore_tkr(trd) B
    DOUBLE PRECISION, DIMENSION(*) :: B
    INTEGER(c_int),VALUE::ldb
    !dir$ ignore_tkr(trd) W
    DOUBLE PRECISION, DIMENSION(*)::W
    !dir$ ignore_tkr(trd) work
    DOUBLE PRECISION, DIMENSION(*) :: work
    INTEGER(c_int),VALUE::lwork
    !dir$ ignore_tkr(trd) devInfo
    INTEGER, DIMENSION(1)::devInfo
  END FUNCTION
  END INTERFACE

  PRIVATE

  PUBLIC :: DPOTRF_OFFLOAD, ZPOTRF_OFFLOAD
  PUBLIC :: DGETRF_OFFLOAD, ZGETRF_OFFLOAD
  PUBLIC :: DGETRS_OFFLOAD, ZGETRS_OFFLOAD
  PUBLIC :: DTRTRI_OFFLOAD, ZTRTRI_OFFLOAD
  PUBLIC :: DSYGV_OFFLOAD, ZHEGV_OFFLOAD
  PUBLIC :: DSYEV_OFFLOAD, ZHEEV_OFFLOAD
  PUBLIC :: SSYEVX_OFFLOAD, DSYEVX_OFFLOAD, ZHEEVX_OFFLOAD, CHEEVX_OFFLOAD

  PUBLIC :: CUSOLVER_INIT, CUSOLVER_FINISH

  CONTAINS

  SUBROUTINE CUSOLVER_INIT()
    INTEGER(4) :: IERR

    PROFILING_START('cusolver_init')

    IF (.NOT.OFFLOAD_ALLOW) THEN
       PROFILING_STOP('cusolver_init')
       RETURN
    ENDIF

!$OMP PARALLEL
    IERR = cusolverDnCreate(ACC_CUSOLVERHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR setting up cusolverHandle")
    ENDIF
!$OMP END PARALLEL

    PROFILING_STOP('cusolver_init')

  END SUBROUTINE CUSOLVER_INIT


  SUBROUTINE CUSOLVER_FINISH()
     INTEGER(4) :: IERR

     IF (.NOT.OFFLOAD_ALLOW) RETURN

!$OMP PARALLEL
     IERR = cusolverDnDestroy(ACC_CUSOLVERHANDLE)
     IF (IERR/=0) THEN
        CALL vtutor%error("ERROR releasing cusolverHandle")
     ENDIF
!$OMP END PARALLEL
  END SUBROUTINE CUSOLVER_FINISH


  SUBROUTINE DPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,UPLO_,INFO_

#ifdef magma_dpotrf
    CALL magmaf_dpotrf_gpu(UPLO,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnDpotrf_bufferSize(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,INFO_)
    IERR = IERR + cusolverDnDpotrf(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,WORK,LWORK,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE DPOTRF_OFFLOAD


  SUBROUTINE ZPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,UPLO_,INFO_

#ifdef magma_zpotrf
    CALL magmaf_zpotrf_gpu(UPLO,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnZpotrf_bufferSize(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,INFO_)
    IERR = IERR + cusolverDnZpotrf(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,WORK,LWORK,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZPOTRF_OFFLOAD


  SUBROUTINE DGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: M,N,LDA,INFO,IPIV(*)
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,INFO_

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnDgetrf_bufferSize(ACC_CUSOLVERHANDLE,M,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,IPIV,INFO_)
    IERR = IERR + cusolverDnDgetrf(ACC_CUSOLVERHANDLE,M,N,A,LDA,WORK,IPIV,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE DGETRF_OFFLOAD


  SUBROUTINE ZGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: M,N,LDA,INFO,IPIV(*)
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,INFO_

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnZgetrf_bufferSize(ACC_CUSOLVERHANDLE,M,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,IPIV,INFO_)
    IERR = IERR + cusolverDnZgetrf(ACC_CUSOLVERHANDLE,M,N,A,LDA,WORK,IPIV,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZGETRF_OFFLOAD


  SUBROUTINE DGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),B(LDB,*)
    INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
    CHARACTER(1) :: TRANS
    ! local variables
    INTEGER(4) :: IERR,TRANS_,INFO_

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (TRANS=='N'.OR.TRANS=='n') THEN
       TRANS_=CUBLAS_OP_N
    ELSE
       TRANS_=CUBLAS_OP_T
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A,IPIV,B,INFO_)
    IERR = IERR + cusolverDnDgetrs(ACC_CUSOLVERHANDLE,TRANS_,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)

    INFO=INFO_*100+IERR
  END SUBROUTINE DGETRS_OFFLOAD


  SUBROUTINE ZGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),B(LDB,*)
    INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
    CHARACTER(1) :: TRANS
    ! local variables
    INTEGER(4) :: IERR,TRANS_,INFO_

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (TRANS=='N'.OR.TRANS=='n') THEN
       TRANS_=CUBLAS_OP_N
    ELSEIF (TRANS=='T'.OR.TRANS=='t') THEN
       TRANS_=CUBLAS_OP_T
    ELSE
       TRANS_=CUBLAS_OP_C
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A,IPIV,B,INFO_)
    IERR = IERR + cusolverDnZgetrs(ACC_CUSOLVERHANDLE,TRANS_,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZGETRS_OFFLOAD


  SUBROUTINE DTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
    USE mcublas
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO,DIAG
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: I,J
    INTEGER(4) :: IERR

#ifdef magma_dtrtri
    CALL magmaf_dtrtri_gpu(UPLO,DIAG,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       ALLOCATE(WORK(N*N))
!$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(WORK) ASYNC(ACC_ASYNC_Q)
       WORK=0
       DO I=1,N
          WORK(I+(I-1)*N)=1
       ENDDO
!$ACC END KERNELS

!$ACC HOST_DATA USE_DEVICE(A,WORK)
       IERR = IERR + cublasDtrsm_v2(ACC_CUBLASHANDLE,CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, &
                                    CUBLAS_OP_N,CUBLAS_DIAG_NON_UNIT,N,N,done,A,LDA,WORK,N)
!$ACC END HOST_DATA

!$ACC KERNELS PRESENT(A,WORK) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          DO J=1,N
             A(J,I)=WORK(J+(I-1)*N)
          ENDDO
       ENDDO
!$ACC END KERNELS

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
       DEALLOCATE(WORK)
    ENDIF

    INFO=IERR
  END SUBROUTINE DTRTRI_OFFLOAD


  SUBROUTINE ZTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
    USE mcublas
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO,DIAG
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: I,J
    INTEGER(4) :: IERR

#ifdef magma_ztrtri
    CALL magmaf_ztrtri_gpu(UPLO,DIAG,N,c_devloc(A),LDA,INFO)
    RETURN
#endif

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       ALLOCATE(WORK(N*N))
!$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(WORK) ASYNC(ACC_ASYNC_Q)
       WORK=0
       DO I=1,N
          WORK(I+(I-1)*N)=1
       ENDDO
!$ACC END KERNELS

!$ACC HOST_DATA USE_DEVICE(A,WORK)
       IERR = IERR + cublasZtrsm_v2(ACC_CUBLASHANDLE,CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, &
                                    CUBLAS_OP_N,CUBLAS_DIAG_NON_UNIT,N,N,zone,A,LDA,WORK,N)
!$ACC END HOST_DATA

!$ACC KERNELS PRESENT(A,WORK) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          DO J=1,N
             A(J,I)=WORK(J+(I-1)*N)
          ENDDO
       ENDDO
!$ACC END KERNELS

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
       DEALLOCATE(WORK)
    ENDIF

    INFO=IERR
  END SUBROUTINE ZTRTRI_OFFLOAD


  SUBROUTINE DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
    INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: ITYPE_,JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

#ifdef magma_dsygv
!   CALL magmaf_dsygv_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (ITYPE==1) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_1
    ELSEIF (ITYPE==2) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_2
    ELSE
       ITYPE_=CUSOLVER_EIG_TYPE_3
    ENDIF

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,B,W)
    IERR = IERR + cusolverDnDsygvd_bufferSize(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsygvd(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE DSYGV_OFFLOAD


  SUBROUTINE ZHEGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),B(LDB,*),WORK(*)
    REAL(q) :: W(*),RWORK(*)
    INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: ITYPE_,JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

#ifdef magma_zhegv
!   CALL magmaf_zhegv_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (ITYPE==1) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_1
    ELSEIF (ITYPE==2) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_2
    ELSE
       ITYPE_=CUSOLVER_EIG_TYPE_3
    ENDIF

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,B,W)
       IERR = IERR + cusolverDnZhegvd_bufferSize(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,LWORK_)
!$ACC END HOST_DATA

       ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZhegvd(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZHEGV_OFFLOAD


  SUBROUTINE DSYEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),W(*),WORK(*)
    INTEGER :: N,LDA,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

#ifdef magma_dsyev
!   CALL magmaf_dsyev_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnDsyevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsyevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE DSYEV_OFFLOAD


  SUBROUTINE ZHEEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),WORK(*)
    REAL(q) :: W(*),RWORK(*)
    INTEGER :: N,LDA,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

#ifdef magma_zheev
!   CALL magmaf_zheev_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnZheevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZheevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZHEEV_OFFLOAD


  SUBROUTINE DSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*),VL,VU,ABSTOL
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
!$ACC DECLARE DEVICE_RESIDENT(WORK_)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_

#ifdef magma_dsyevx
!   CALL magmaf_dsyevx_gpu
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("DSYEVX_OFFLOAD: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnDsyevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsyevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE DSYEVX_OFFLOAD


  SUBROUTINE ZHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),Z(LDZ,*),WORK(*)
    REAL(q) :: W(*),RWORK(*),VL,VU,ABSTOL
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    COMPLEX(q), POINTER :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_
#ifdef magma_zheevx
    COMPLEX(q) :: WA(N,N),WZ(N,M)
#endif

#ifdef magma_zheevx
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_zheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK,     -1,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    LWORK_=WORK(1); ALLOCATE(WORK_(LWORK_))
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_zheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK_,LWORK_,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    DEALLOCATE(WORK_)
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("ZHEEVX_OFFLOAD: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnZheevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA

    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_,WORK_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZheevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ZHEEVX_OFFLOAD


  SUBROUTINE CHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
    IMPLICIT NONE
    COMPLEX(qs) :: A(LDA,*),Z(LDZ,*),WORK(*)
    REAL(qs) :: W(*),RWORK(*),VL,VU,ABSTOL
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    COMPLEX(qs), POINTER :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_
#ifdef magma_cheevx
    COMPLEX(qs) :: WA(N,N),WZ(N,M)
#endif

#ifdef magma_cheevx
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_cheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK,     -1,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    LWORK_=WORK(1); ALLOCATE(WORK_(LWORK_))
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_cheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK_,LWORK_,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    DEALLOCATE(WORK_)
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("CHEEVX_OFFLOAD: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnCheevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA

    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_,WORK_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnCheevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE CHEEVX_OFFLOAD


  SUBROUTINE SSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
    IMPLICIT NONE
    REAL(qs) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
    REAL(qs) :: ABSTOL,VL,VU
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    REAL(qs), ALLOCATABLE :: WORK_(:)
!$ACC DECLARE DEVICE_RESIDENT(WORK_)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_

#ifdef magma_dsyevx
!   CALL magmaf_dsyevx_gpu
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("DSYEVX_OFFLOAD: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnSsyevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnSsyevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE SSYEVX_OFFLOAD

  END MODULE mcusolver

!***********************************************************************
!
!> This module holds everything related to the use of cal/nccl_comm
!
!***********************************************************************
  MODULE mmp_comm
#if PGI_BEFORE_25_9
#ifdef CUSOLVERMP
  USE cusolvermp, ONLY : cal_comm_create_mpi, cal_comm_destroy, mpcomm => cal_comm
#endif
#ifdef CUBLASMP
  USE cublasmp, ONLY : cal_comm_create_mpi, cal_comm_destroy, mpcomm => cal_comm
#endif
#else
  USE nccl, mpcomm => ncclComm
#endif
  USE tutor, ONLY: vtutor
  USE string, ONLY: str

  IMPLICIT NONE
  PRIVATE

#if defined(CUSOLVERMP) || defined(CUBLASMP)
  PUBLIC :: mpcomm, GET_MP_COMM,  MP_COMM_FREE

  TYPE dt_mp_comm_struct
    INTEGER(4)     :: MPICOMM
    INTEGER        :: MYRANK
    INTEGER        :: MYRANKS
    TYPE(mpcomm)   :: ACC_MP_COMM
    INTEGER,ALLOCATABLE :: RANKS(:)
  END TYPE

  INTEGER, PARAMETER        :: MAX_COM = 1024
  TYPE(dt_mp_comm_struct)   :: COMMS(MAX_COM)
  INTEGER                   :: NCOMM = 0

  CONTAINS

  SUBROUTINE GET_RANKS(RANKS, mpicomm, SZ)
    USE main_mpi
    INTEGER :: RANKS(SZ)
    INTEGER :: MPICOMM
    INTEGER :: IERR, SZ
    INTEGER :: R(SZ)
    INTEGER :: GRP1, GRP2
    INTEGER :: I

    CALL MPI_COMM_GROUP(COMM_WORLD%MPI_COMM, GRP1, IERR)
    CALL MPI_COMM_GROUP(            MPICOMM, GRP2, IERR)
    CALL MPI_COMM_SIZE(MPICOMM, SZ, IERR)

    DO I=0, SZ-1
       R(I+1) = I
    ENDDO

    CALL MPI_GROUP_TRANSLATE_RANKS(GRP2, SZ, R, GRP1, RANKS(1), IERR)
    CALL MPI_GROUP_FREE(GRP1, IERR)
    CALL MPI_GROUP_FREE(GRP2, IERR)
  END SUBROUTINE

  TYPE(mpcomm) FUNCTION GET_MP_COMM(CTXT)
    use openacc
    use mopenacc_struct_def
    USE main_mpi
    IMPLICIT NONE
    INTEGER :: CTXT
    INTEGER :: BCTXT, MPICOMM, MYRANK, MYRANKS
    INTEGER :: I, RES, IERR
    INTEGER, EXTERNAL :: BLACS2SYS_HANDLE
    TYPE(mpcomm) :: ACC_MP_COMM
#if ! PGI_BEFORE_25_9
    TYPE(communic) :: COMM_TMP
#endif
    INTEGER, ALLOCATABLE :: RANKS(:)

    PROFILING_START('get_mp_comm')
    CALL BLACS_GET(CTXT, 10, BCTXT)
    MPICOMM = BLACS2SYS_HANDLE(BCTXT)
    CALL MPI_COMM_RANK(MPICOMM, MYRANK, IERR)
    CALL MPI_COMM_SIZE(MPICOMM, MYRANKS, IERR)

    ALLOCATE(RANKS(MYRANKS))
    CALL GET_RANKS(RANKS, mpicomm, MYRANKS)

    IF(NCOMM /= 0) THEN
        DO I = 1,NCOMM
           IF(COMMS(I)%MYRANKS /= MYRANKS .OR. &
              COMMS(I)%MYRANK  /= MYRANK  .OR. &
              ANY(COMMS(I)%RANKS /= RANKS))  CYCLE

           GET_MP_COMM = COMMS(I)%ACC_MP_COMM
           DEALLOCATE(RANKS)
           PROFILING_STOP('get_mp_comm')
           RETURN
        ENDDO
    ENDIF

    IF(NCOMM == MAX_COM) THEN
       CALL vtutor%error("GET_MP_COMM: error: the number of mp comms are limited by 1024" )
    ENDIF
#if PGI_BEFORE_25_9
    IERR = cal_comm_create_mpi(MPICOMM, MYRANK, MYRANKS, acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_MP_COMM)
    IF ( IERR /= 0 ) CALL vtutor%error("GET_MP_COMM: error in cal_comm_create_mpi: " // str(IERR) )
#else
#ifndef USENCCL
#error "cuSolverMp and cuBlasMp libs required NCCL to be enabled"
#endif
    IF ( .NOT. LUSENCCL ) CALL vtutor%warning('Using NCCL is disabled')

    COMM_TMP%NODE_ME = MYRANK+1
    COMM_TMP%IONODE  = 1
    COMM_TMP%NCPU    = MYRANKS
    COMM_TMP%MPI_COMM = MPICOMM
    COMM_TMP%LUSENCCL = .TRUE.
    CALL nccl_init_comm(COMM_TMP)

    IF ( .NOT. COMM_TMP%LUSENCCL) THEN
        CALL vtutor%bug('Failed to create NCCL communicator ', __FILE__, __LINE__)
    ENDIF
    ACC_MP_COMM = COMM_TMP%NCCL_COMM
#endif
    NCOMM = NCOMM + 1
    COMMS(NCOMM)%MPICOMM = MPICOMM
    COMMS(NCOMM)%MYRANK  = MYRANK
    COMMS(NCOMM)%MYRANKS = MYRANKS
    COMMS(NCOMM)%ACC_MP_COMM = ACC_MP_COMM
    ALLOCATE(COMMS(NCOMM)%RANKS(MYRANKS))
    COMMS(NCOMM)%RANKS = RANKS
    DEALLOCATE(RANKS)
    GET_MP_COMM = ACC_MP_COMM
    PROFILING_STOP('get_mp_comm')
  END FUNCTION


  SUBROUTINE MP_COMM_FREE()
      USE main_mpi
      IMPLICIT NONE
      INTEGER :: I
#if PGI_BEFORE_25_9
      INTEGER :: IERR
#else
      TYPE(ncclResult)  :: ncclRes
#endif

      IF(NCOMM == 0) THEN
         RETURN
      ENDIF

      DO I = 1, NCOMM
         DEALLOCATE(COMMS(I)%RANKS)
#if PGI_BEFORE_25_9
         !IERR = cal_comm_destroy(COMMS(I)%ACC_MP_COMM)
         !IF ( IERR /= 0 ) CALL vtutor%error("MP_COMM_FREE: error in cal_comm_destroy: " // str(IERR) )
#else
         ncclRes = ncclCommDestroy(COMMS(I)%ACC_MP_COMM)
         IF ( ncclRes /= ncclSuccess ) CALL vtutor%error("MP_COMM_FREE: error in ncclCommDestroy" )
#endif
      ENDDO
      NCOMM = 0
  END SUBROUTINE
#endif
  END MODULE mmp_comm

#define CHK_CUBLASMP_ERR(err, msg)  IF ( err%stat /= 0 ) CALL vtutor%error(msg // str(err%stat) )

  MODULE mmp_grid
#ifdef CUBLASMP
  USE cublasmp
#endif

#if defined(CUSOLVERMP) || defined(CUBLASMP)
  PUBLIC :: GET_CUBLASMP_GRID, MP_GRID_FREE

  TYPE dt_mp_grid_struct
    INTEGER(4)     :: MPICOMM
    INTEGER        :: MYRANK
    INTEGER        :: MYRANKS
    INTEGER        :: NPROW
    INTEGER        :: NPCOL
#ifdef CUBLASMP
    TYPE(cublasMpGridLayout) :: MAPPING
    TYPE(cublasMpGrid) :: CUBLASMP_GRID
#endif
    INTEGER, ALLOCATABLE :: RANKS(:)
  END TYPE

  INTEGER, PARAMETER        :: MAX_GRID = 100
  TYPE(dt_mp_grid_struct)   :: GRIDS(MAX_GRID)
  INTEGER                   :: NGRID = 0

  CONTAINS

#ifdef CUBLASMP
  TYPE(cublasMpGrid) FUNCTION GET_CUBLASMP_GRID(CTXT, CUBLASMP_HANDLER, MAPPING)
    USE mmp_comm
    USE main_mpi
    USE cublasmp
    use tutor, ONLY: vtutor
    use string, ONLY: str

    IMPLICIT NONE
    INTEGER                  :: CTXT
    TYPE(cublasMpHandle)     :: CUBLASMP_HANDLER
    TYPE(cublasMpGridLayout) :: MAPPING
    TYPE(cublasMpStatus)     :: IERR
    TYPE(mpcomm) :: ACC_MP_COMM
    INTEGER :: BCTXT, MPICOMM, MYRANK, MYRANKS
    INTEGER :: MYROW, MYCOL, NPROW, NPCOL
    INTEGER :: I,RES
    TYPE(cublasMpGrid) :: CUBLASMP_GRID
    INTEGER, EXTERNAL :: BLACS2SYS_HANDLE
    INTEGER, ALLOCATABLE :: RANKS(:)

    GET_CUBLASMP_GRID = cublasMpGrid(C_NULL_PTR)
    IF(CTXT == -1 ) RETURN

    PROFILING_START('GET_CUBLASMP_GRID')
    CALL BLACS_GRIDINFO(CTXT, NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK
    CALL BLACS_GET(CTXT, 10, BCTXT)
    MPICOMM = BLACS2SYS_HANDLE(BCTXT)
    CALL MPI_COMM_RANK(MPICOMM, MYRANK, I)
    CALL MPI_COMM_SIZE(MPICOMM, MYRANKS, I)
    ALLOCATE(RANKS(MYRANKS))
    CALL MPI_Allgather(COMM_WORLD%NODE_ME-1, 1, MPI_INTEGER, RANKS(1), 1, MPI_INTEGER, MPICOMM, IERR)

    IF(NGRID > 0) THEN
       DO I=1, NGRID
           IF(GRIDS(I)%MYRANKS /= MYRANKS .OR. &
              GRIDS(I)%MYRANK  /= MYRANK  .OR. &
              GRIDS(I)%NPROW   /= NPROW   .OR. &
              GRIDS(I)%NPCOL   /= NPCOL   .OR. &
              GRIDS(I)%MAPPING /= MAPPING .OR. &
              ANY(GRIDS(I)%RANKS /= RANKS))  CYCLE

           GET_CUBLASMP_GRID = GRIDS(I)%CUBLASMP_GRID
           DEALLOCATE(RANKS)
           PROFILING_STOP('GET_CUBLASMP_GRID')
           RETURN
       ENDDO
    ENDIF

    IF(NGRID == MAX_GRID) THEN
        CALL vtutor%error("GET_CUBLASMP_GRID: error: the number of mp grids are limited by 100" )
    ENDIF

    ! make sure we have an appropriate communicator
    ACC_MP_COMM = GET_MP_COMM(CTXT)

    IERR = cublasMpGridCreate(NPROW, NPCOL, MAPPING, ACC_MP_COMM, CUBLASMP_GRID)

    CHK_CUBLASMP_ERR(IERR, "GET_CUBLASMP_GRID: error in cublasMpGridCreate: ")

    NGRID = NGRID + 1

    GRIDS(NGRID)%MYRANKS = MYRANKS
    GRIDS(NGRID)%MYRANK  = MYRANK
    GRIDS(NGRID)%NPROW   = NPROW
    GRIDS(NGRID)%NPCOL   = NPCOL
    GRIDS(NGRID)%MAPPING = MAPPING
    GRIDS(NGRID)%MPICOMM = MPICOMM
    GRIDS(NGRID)%CUBLASMP_GRID = CUBLASMP_GRID
    ALLOCATE(GRIDS(NGRID)%RANKS(MYRANKS))
    GRIDS(NGRID)%RANKS = RANKS
    GET_CUBLASMP_GRID = GRIDS(NGRID)%CUBLASMP_GRID
    DEALLOCATE(RANKS)
    PROFILING_STOP('GET_CUBLASMP_GRID')
  END FUNCTION GET_CUBLASMP_GRID
#endif

  SUBROUTINE MP_GRID_FREE()
    use tutor, ONLY: vtutor
    use string, ONLY: str
      IMPLICIT NONE
      INTEGER :: I
      TYPE(cublasMpStatus)  :: IERR

      IF(NGRID == 0) THEN
         RETURN
      ENDIF

      DO I = 1, NGRID
         DEALLOCATE(GRIDS(I)%RANKS)
         IERR = cublasMpGridDestroy(GRIDS(I)%CUBLASMP_GRID)
         CHK_CUBLASMP_ERR(IERR, "MP_GRID_FREE: error in cublasMpGridDestroy: ")
      ENDDO
      NGRID = 0
  END SUBROUTINE

#endif ! defined(CUSOLVERMP) || defined(CUBLASMP)
  END MODULE

!***********************************************************************
!
!> This module holds everything related to the use of cublasmp.
!
!***********************************************************************

  MODULE mcublasmp

  USE prec
  USE mopenacc_struct_def
  USE scala_struct_def
  USE mcublas, ONLY : CHAR_TO_FILL, CHAR_TO_DIAG
  USE openacc, ONLY : ACC_GET_DEVICE_NUM, ACC_DEVICE_NVIDIA, ACC_GET_CUDA_STREAM
  USE cusolverDn, ONLY : CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR
#ifdef CUBLASMP
  USE mmp_comm, ONLY : mpcomm, GET_MP_COMM, MP_COMM_FREE
  USE mmp_grid
  USE cublasmp
#endif

  USE string, ONLY : str
  USE tutor, ONLY : vtutor

  PRIVATE

  PUBLIC :: PDTRTRI_OFFLOAD, &
            PZTRTRI_OFFLOAD, &
            PSTRMR2D_OFFLOAD, &
            PDTRMR2D_OFFLOAD, &
            PCTRMR2D_OFFLOAD, &
            PZTRMR2D_OFFLOAD, &
            PSGEMR2D_OFFLOAD, &
            PDGEMR2D_OFFLOAD, &
            PCGEMR2D_OFFLOAD, &
            PZGEMR2D_OFFLOAD, &
            CUBLASMP_FREE

#ifdef CUBLASMP
   TYPE dt_cublasmp_struct
      INTEGER              :: acc_stream
      TYPE(cublasMpHandle) :: ACC_CUBLASMP_HANDLE
   END TYPE

   INTEGER, PARAMETER        :: MAX_CUBLASMP_HANDLERS = 128
   TYPE(dt_cublasmp_struct)  :: CUBLASMP_HANDLERS(MAX_CUBLASMP_HANDLERS)
   INTEGER                   :: NCUBLASMPHANDLERS = 0
#endif

  CONTAINS

#ifdef CUBLASMP
  TYPE(cublasMpHandle) FUNCTION GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)
    use cublasmp
    use tutor, ONLY: vtutor
    use string, ONLY: str
    use openacc
    IMPLICIT NONE
    INTEGER              :: I, ACC_ASYNC_Q
    TYPE(cublasMpHandle) :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpStatus) :: IERR

    PROFILING_START('GET_CUBLASMP_HANDLER')
    IF(NCUBLASMPHANDLERS > 0) THEN
    DO I=1,NCUBLASMPHANDLERS
       IF(CUBLASMP_HANDLERS(I)%acc_stream == ACC_ASYNC_Q) THEN
          GET_CUBLASMP_HANDLER = CUBLASMP_HANDLERS(I)%ACC_CUBLASMP_HANDLE
          PROFILING_STOP('GET_CUBLASMP_HANDLER')
          return
       ENDIF
    ENDDO
    ENDIF
    IF(NCUBLASMPHANDLERS == MAX_CUBLASMP_HANDLERS) THEN
       CALL vtutor%error("GET_CUBLASMP_HANDLER: error: the number of cublasmp handlers are limited by " // str(MAX_CUBLASMP_HANDLERS) )
    ENDIF
    IERR = cublasMpCreate(ACC_CUBLASMP_HANDLE, acc_get_cuda_stream(ACC_ASYNC_Q))
    CHK_CUBLASMP_ERR(IERR, "GET_CUBLASMP_HANDLER: error in cublasMpCreate: ")

    NCUBLASMPHANDLERS = NCUBLASMPHANDLERS + 1
    CUBLASMP_HANDLERS(NCUBLASMPHANDLERS)%acc_stream = ACC_ASYNC_Q
    CUBLASMP_HANDLERS(NCUBLASMPHANDLERS)%ACC_CUBLASMP_HANDLE = ACC_CUBLASMP_HANDLE
    GET_CUBLASMP_HANDLER = ACC_CUBLASMP_HANDLE
    PROFILING_STOP('GET_CUBLASMP_HANDLER')
  END FUNCTION

  SUBROUTINE CUBLASMP_FREE()
    use tutor, ONLY: vtutor
    use string, ONLY: str
    use cublasmp
    IMPLICIT NONE
    INTEGER :: I
    TYPE(cublasMpStatus) :: IERR

    IF(NCUBLASMPHANDLERS == 0) THEN
       RETURN
    ENDIF

    DO I = 1, NCUBLASMPHANDLERS
       IERR = cublasMpDestroy(CUBLASMP_HANDLERS(I)%ACC_CUBLASMP_HANDLE)
       CHK_CUBLASMP_ERR(IERR, "CUBLASMP_FREE: error in cublasMpDestroy: ")
    ENDDO
    NCUBLASMPHANDLERS = 0
  END SUBROUTINE
#endif


  SUBROUTINE PDTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
    USE scala, ONLY: SET_DIAGONALE_Z, SET_DIAGONALE_D

    IMPLICIT NONE

    REAL(q) :: A(*)
    INTEGER :: DESCA(*)
    INTEGER :: IA, INFO, JA, N
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cublasMpGridLayout) :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)     :: IERR
    TYPE(cublasMpHandle)     :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)       :: ACC_CUBLASMP_GRID_A
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA
    TYPE(mpcomm)             :: ACC_MP_COMM_A
    TYPE(cudaDataType)       :: data_type = cudaDataType(CUDA_R_64F)
    TYPE(cublasComputeType)  :: compute_type = cublasComputeType(CUBLAS_COMPUTE_64F)
    INTEGER(8)               :: devworksize, hostworksize
    INTEGER(1), ALLOCATABLE  :: devwork(:), hostwork(:)
    REAL(q), ALLOCATABLE  :: B(:,:)
    REAL(q), ALLOCATABLE     :: AD(:)
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC

    PROFILING_START('PDTRTRI_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ! prepare Descriptor for matrix A
    ACC_MP_COMM_A = GET_MP_COMM(DESCA(CTXT_))
    CALL BLACS_GRIDINFO( DESCA(CTXT_), NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK

    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PDTRTRI_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    NP = MAX(1,NUMROC(N,DESCA(MB_),MYROW,0,NPROW))
    NQ = MAX(1,NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL))

    ALLOCATE(B(NP,NQ))
    ALLOCATE(AD(N))

    ! set diagonal to one
!$ACC ENTER DATA CREATE(B,AD) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(B,AD) ASYNC(ACC_ASYNC_Q)
    B=dzero
    AD=done
!$ACC END KERNELS
    CALL SET_DIAGONALE_D(N, B, AD, DESCA)
!$ACC EXIT DATA DELETE(AD) ASYNC(ACC_ASYNC_Q)

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrsm_bufferSize(ACC_CUBLASMP_HANDLE, CUBLAS_SIDE_LEFT, CHAR_TO_FILL(UPLO), &
                                   CUBLAS_OP_N,CHAR_TO_DIAG(DIAG), N, N, &
                                   done, A,  IA, JA, ACC_CUBLASMP_DESCA, B, IA, JA, ACC_CUBLASMP_DESCA, &
                                   compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDTRTRI_OFFLOAD: error in cublasMpTrsm_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrsm(ACC_CUBLASMP_HANDLE, CUBLAS_SIDE_LEFT,CHAR_TO_FILL(UPLO), &
                                     CUBLAS_OP_N,CHAR_TO_DIAG(DIAG), &
                                     N,N, done, A,  IA, JA, ACC_CUBLASMP_DESCA, B, IA, JA, ACC_CUBLASMP_DESCA, &
                                     compute_type, devwork, devworksize, hostwork, hostworksize)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDTRTRI_OFFLOAD: error in cublasMptrsm: ")

    CALL WDCOPY(NQ*NP,B,1,A,1)

    ! deallocate memory
!$ACC EXIT DATA DELETE(B,devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(B,devwork,hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PDTRTRI_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PDTRTRI_OFFLOAD')
#else
    CALL vtutor%bug("PDTRTRI_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDTRTRI_OFFLOAD


  SUBROUTINE PZTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
    USE scala, ONLY: SET_DIAGONALE_Z, SET_DIAGONALE_D

    IMPLICIT NONE

    COMPLEX(q) :: A(*)
    INTEGER :: DESCA(*)
    INTEGER :: IA, INFO, JA, N
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cublasMpGridLayout) :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)     :: IERR
    TYPE(cublasMpHandle)     :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)       :: ACC_CUBLASMP_GRID_A
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA
    TYPE(mpcomm)             :: ACC_MP_COMM_A
    TYPE(cudaDataType)       :: data_type = cudaDataType(CUDA_C_64F)
    TYPE(cublasComputeType)  :: compute_type = cublasComputeType(CUBLAS_COMPUTE_64F)
    INTEGER(8)               :: devworksize, hostworksize
    INTEGER(1), ALLOCATABLE  :: devwork(:), hostwork(:)
    COMPLEX(q), ALLOCATABLE  :: B(:,:)
    REAL(q), ALLOCATABLE     :: AD(:)
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC

    PROFILING_START('PZTRTRI_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ! prepare Descriptor for matrix A
    ACC_MP_COMM_A = GET_MP_COMM(DESCA(CTXT_))
    CALL BLACS_GRIDINFO( DESCA(CTXT_), NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK


    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PZTRTRI_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    NP = MAX(1,NUMROC(N,DESCA(MB_),MYROW,0,NPROW))
    NQ = MAX(1,NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL))

    ALLOCATE(B(NP,NQ))
    ALLOCATE(AD(N))

    ! set diagonal to one
!$ACC ENTER DATA CREATE(B,AD) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(B,AD) ASYNC(ACC_ASYNC_Q)
    B=zzero
    AD=zone
!$ACC END KERNELS
    CALL SET_DIAGONALE_Z(N, B, AD, DESCA)
!$ACC EXIT DATA DELETE(AD) ASYNC(ACC_ASYNC_Q)

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrsm_bufferSize(ACC_CUBLASMP_HANDLE, CUBLAS_SIDE_LEFT, CHAR_TO_FILL(UPLO), &
                                   CUBLAS_OP_N,CHAR_TO_DIAG(DIAG), N, N, &
                                   zone, A,  IA, JA, ACC_CUBLASMP_DESCA, B, IA, JA, ACC_CUBLASMP_DESCA, &
                                   compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZTRTRI_OFFLOAD: error in cublasMpTrsm_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrsm(ACC_CUBLASMP_HANDLE, CUBLAS_SIDE_LEFT,CHAR_TO_FILL(UPLO), &
                                     CUBLAS_OP_N,CHAR_TO_DIAG(DIAG), &
                                     N,N, zone, A,  IA, JA, ACC_CUBLASMP_DESCA, B, IA, JA, ACC_CUBLASMP_DESCA, &
                                     compute_type,devwork, devworksize, hostwork, hostworksize)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZTRTRI_OFFLOAD: error in cublasMptrsm: ")

    CALL WZCOPY(NQ*NP,B,1,A,1)

    ! deallocate memory
!$ACC EXIT DATA DELETE(B,devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(B,devwork,hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PZTRTRI_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PZTRTRI_OFFLOAD')
#else
    CALL vtutor%bug("PZTRTRI_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZTRTRI_OFFLOAD


  SUBROUTINE PSTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    REAL(qs) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_R_32F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: UPLO_, DIAG_
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PSTRMR2D_OFFLOAD')

    IF (DIAG=='N'.OR.DIAG=='n') THEN
       DIAG_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       DIAG_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocate scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrmr2d_bufferSize(ACC_CUBLASMP_HANDLE,UPLO_,DIAG_, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpTrmr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrmr2d(ACC_CUBLASMP_HANDLE, UPLO_, DIAG_, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PSTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PSTRMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PSTRMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PSTRMR2D_OFFLOAD


  SUBROUTINE PDTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    REAL(q) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: UPLO_, DIAG_
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PDTRMR2D_OFFLOAD')

    IF (DIAG=='N'.OR.DIAG=='n') THEN
       DIAG_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       DIAG_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocate scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrmr2d_bufferSize(ACC_CUBLASMP_HANDLE,UPLO_,DIAG_, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpTrmr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrmr2d(ACC_CUBLASMP_HANDLE, UPLO_, DIAG_, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PDTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PDTRMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PDTRMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDTRMR2D_OFFLOAD


  SUBROUTINE PCTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    COMPLEX(qs) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_C_32F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: UPLO_, DIAG_
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PCTRMR2D_OFFLOAD')

    IF (DIAG=='N'.OR.DIAG=='n') THEN
       DIAG_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       DIAG_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocate scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrmr2d_bufferSize(ACC_CUBLASMP_HANDLE,UPLO_,DIAG_, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpTrmr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrmr2d(ACC_CUBLASMP_HANDLE, UPLO_, DIAG_, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PCTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PCTRMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PCTRMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PCTRMR2D_OFFLOAD


  SUBROUTINE PZTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    COMPLEX(q) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT
    CHARACTER(1) :: UPLO, DIAG

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: UPLO_, DIAG_
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PZTRMR2D_OFFLOAD')

    IF (DIAG=='N'.OR.DIAG=='n') THEN
       DIAG_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       DIAG_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocate scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpTrmr2d_bufferSize(ACC_CUBLASMP_HANDLE,UPLO_,DIAG_, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpTrmr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpTrmr2d(ACC_CUBLASMP_HANDLE, UPLO_, DIAG_, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PZTRMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PZTRMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PZTRMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZTRMR2D_OFFLOAD


  SUBROUTINE PSGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    REAL(qs) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_R_32F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PSGEMR2D_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpGemr2d_bufferSize(ACC_CUBLASMP_HANDLE, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpGemr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpGemr2d(ACC_CUBLASMP_HANDLE, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PSGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PSGEMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PSGEMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PSGEMR2D_OFFLOAD


  SUBROUTINE PDGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    REAL(q) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PDGEMR2D_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpGemr2d_bufferSize(ACC_CUBLASMP_HANDLE, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpGemr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpGemr2d(ACC_CUBLASMP_HANDLE, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PDGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PDGEMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PDGEMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDGEMR2D_OFFLOAD


  SUBROUTINE PZGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    COMPLEX(q) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PZGEMR2D_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpGemr2d_bufferSize(ACC_CUBLASMP_HANDLE, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpGemr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpGemr2d(ACC_CUBLASMP_HANDLE, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PZGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PZGEMR2D_OFFLOAD')
#else
    CALL vtutor%bug("PZGEMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZGEMR2D_OFFLOAD


  SUBROUTINE PCGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

    IMPLICIT NONE

    COMPLEX(qs) :: A(*), B(*)
    INTEGER :: DESCA(*), DESCB(*)
    INTEGER :: IA, JA, IB, JB, M, N, ICTXT

#ifdef CUBLASMP
    TYPE(cudaDataType)             :: data_type = cudaDataType(CUDA_C_32F)
    INTEGER(1), ALLOCATABLE        :: devwork(:), hostwork(:)
    TYPE(cublasMpGridLayout)       :: MAPPING = CUBLASMP_GRID_LAYOUT_ROW_MAJOR
    TYPE(cublasMpStatus)           :: IERR
    TYPE(cublasMpHandle)           :: ACC_CUBLASMP_HANDLE
    TYPE(cublasMpGrid)             :: ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_GRID_B
    TYPE(cublasMpMatrixDescriptor) :: ACC_CUBLASMP_DESCA, ACC_CUBLASMP_DESCB
    TYPE(mpcomm)                   :: ACC_MP_COMM_A, ACC_MP_COMM_B, ACC_GLOBAL_COMM
    INTEGER(8)                     :: devworksize, hostworksize
    INTEGER                        :: NPROW, NPCOL, MYPROW, MYPCOL

    PROFILING_START('PCGEMR2D_OFFLOAD')

    ACC_CUBLASMP_HANDLE = GET_CUBLASMP_HANDLER(ACC_ASYNC_Q)

    ACC_GLOBAL_COMM = GET_MP_COMM(ICTXT)

    ! prepare Descriptor for matrix A
    ACC_CUBLASMP_GRID_A = GET_CUBLASMP_GRID(DESCA(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCA(M_), DESCA(N_), DESCA(MB_), DESCA(NB_), &
                                          DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_A, ACC_CUBLASMP_DESCA)

    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! prepare Descriptor for matrix B
    ACC_CUBLASMP_GRID_B = GET_CUBLASMP_GRID(DESCB(CTXT_), ACC_CUBLASMP_HANDLE, MAPPING)

    IERR = cublasMpMatrixDescriptorCreate(DESCB(M_), DESCB(N_), DESCB(MB_), DESCB(NB_), &
                                          DESCB(RSRC_), DESCB(CSRC_), DESCB(LLD_),      &
                                          data_type, ACC_CUBLASMP_GRID_B, ACC_CUBLASMP_DESCB)

    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorCreate: ")

    ! calculate and allocat scratch space
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasMpGemr2d_bufferSize(ACC_CUBLASMP_HANDLE, M, N,     &
                                     A, IA, JA, ACC_CUBLASMP_DESCA, &
                                     B, IB, JB, ACC_CUBLASMP_DESCB, &
                                     devworksize, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpGemr2d_bufferSize: ")

    ! allocate memory
    ALLOCATE(hostwork(hostworksize))
    ALLOCATE(devwork(devworksize))
!$ACC ENTER DATA CREATE(devwork) ASYNC(ACC_ASYNC_Q)

    ! perform redistibution
!$ACC HOST_DATA USE_DEVICE(A,B,devwork)
    IERR = cublasMpGemr2d(ACC_CUBLASMP_HANDLE, M, N,     &
                           A, IA, JA, ACC_CUBLASMP_DESCA, &
                           B, IB, JB, ACC_CUBLASMP_DESCB, &
                           devwork, devworksize, hostwork, hostworksize, ACC_GLOBAL_COMM)
!$ACC END HOST_DATA
    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpGemr2d: ")
    ! deallocate memory
!$ACC EXIT DATA DELETE(devwork) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(devwork)
    DEALLOCATE(hostwork)

    !destroy descriptors
    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCA)
    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    IERR = cublasMpMatrixDescriptorDestroy(ACC_CUBLASMP_DESCB)
    CHK_CUBLASMP_ERR(IERR, "PCGEMR2D_OFFLOAD: error in cublasMpMatrixDescriptorDestroy: ")

    PROFILING_STOP('PCGEMR2D_OFFLOAD')

#else
    CALL vtutor%bug("PCGEMR2D_OFFLOAD: CUBLASMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PCGEMR2D_OFFLOAD


  END MODULE mcublasmp


!***********************************************************************
!
!> This module holds everything related to the use of cusolvermp.
!
!***********************************************************************

  MODULE mcusolvermp

    USE prec
    USE mopenacc_struct_def
    USE scala_struct_def
    USE mcublas, ONLY : CHAR_TO_FILL, CHAR_TO_DIAG
    USE openacc, ONLY : ACC_GET_DEVICE_NUM, ACC_DEVICE_NVIDIA, ACC_GET_CUDA_STREAM
    USE cusolverDn, ONLY : CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR
#ifdef CUSOLVERMP
    USE mmp_comm, ONLY : mpcomm, GET_MP_COMM, MP_COMM_FREE
    USE cusolvermp
#endif

    USE string, ONLY : str
    USE tutor, ONLY : vtutor

    PRIVATE

    PUBLIC :: PDPOTRF_OFFLOAD, &
              PZPOTRF_OFFLOAD, &
              PSSYEVX_OFFLOAD, &
              PDSYEVX_OFFLOAD, &
              PCHEEVX_OFFLOAD, &
              PZHEEVX_OFFLOAD, &
              PSSYEVD_OFFLOAD, &
              PDSYEVD_OFFLOAD, &
              PCHEEVD_OFFLOAD, &
              PZHEEVD_OFFLOAD, &
              PSSYEV_OFFLOAD, &
              PDSYEV_OFFLOAD, &
              PCHEEV_OFFLOAD, &
              PZHEEV_OFFLOAD, &
              CUSOLVERMP_FREE

#ifdef CUSOLVERMP
    TYPE dt_cusolvermp_struct
       INTEGER                :: acc_stream
       INTEGER                :: device_num
       TYPE(cusolverMpHandle) :: ACC_CUSOLVERMP_HANDLE
    END TYPE

    INTEGER, PARAMETER        :: MAX_CUSOLVERMP_HANDLERS = 128
    TYPE(dt_cusolvermp_struct):: CUSOLVERMP_HANDLERS(MAX_CUSOLVERMP_HANDLERS)
    INTEGER                   :: NCUSOLVERMPHANDLERS = 0
#endif

  CONTAINS

#ifdef CUSOLVERMP
  TYPE(cusolverMpHandle) FUNCTION GET_CUSOLVERMP_HANDLER(DEV_NUM, ACC_ASYNC_Q)
      use cublasmp
      use tutor, ONLY: vtutor
      use string, ONLY: str
      use openacc
      IMPLICIT NONE
      INTEGER                :: I, IERR, ACC_ASYNC_Q, DEV_NUM
      TYPE(cusolverMpHandle) :: ACC_CUSOLVERMP_HANDLE

      PROFILING_START('GET_CUSOLVERMP_HANDLER')
      IF(NCUSOLVERMPHANDLERS > 0) THEN
      DO I=1,NCUSOLVERMPHANDLERS
         IF(CUSOLVERMP_HANDLERS(I)%acc_stream == ACC_ASYNC_Q .AND. &
            CUSOLVERMP_HANDLERS(I)%device_num == DEV_NUM ) THEN
              GET_CUSOLVERMP_HANDLER = CUSOLVERMP_HANDLERS(I)%ACC_CUSOLVERMP_HANDLE
              PROFILING_STOP('GET_CUSOLVERMP_HANDLER')
              return
         ENDIF
      ENDDO
      ENDIF
      IF(NCUSOLVERMPHANDLERS == MAX_CUSOLVERMP_HANDLERS) THEN
          CALL vtutor%error("GET_CUSOLVERMP_HANDLER: error: the number of cublasmp handlers are limited by " // str(MAX_CUSOLVERMP_HANDLERS) )
      ENDIF
      IERR = cusolverMpCreate(ACC_CUSOLVERMP_HANDLE, DEV_NUM, acc_get_cuda_stream(ACC_ASYNC_Q))
      IF ( IERR /= 0 ) CALL vtutor%error("GET_CUSOLVERMP_HANDLER: error in cusolverMpCreate: " // str(IERR) )

      NCUSOLVERMPHANDLERS = NCUSOLVERMPHANDLERS + 1
      CUSOLVERMP_HANDLERS(NCUSOLVERMPHANDLERS)%acc_stream = ACC_ASYNC_Q
      CUSOLVERMP_HANDLERS(NCUSOLVERMPHANDLERS)%device_num = DEV_NUM
      CUSOLVERMP_HANDLERS(NCUSOLVERMPHANDLERS)%ACC_CUSOLVERMP_HANDLE = ACC_CUSOLVERMP_HANDLE
      GET_CUSOLVERMP_HANDLER = ACC_CUSOLVERMP_HANDLE
      PROFILING_STOP('GET_CUSOLVERMP_HANDLER')
  END FUNCTION


  SUBROUTINE CUSOLVERMP_FREE()
      use tutor, ONLY: vtutor
      use string, ONLY: str
      use cusolvermp
      IMPLICIT NONE
      INTEGER :: I, IERR

      IF(NCUSOLVERMPHANDLERS == 0) THEN
         RETURN
      ENDIF

      DO I = 1, NCUSOLVERMPHANDLERS
         IERR = cusolverMpDestroy(CUSOLVERMP_HANDLERS(I)%ACC_CUSOLVERMP_HANDLE)
         IF ( IERR /= 0 ) CALL vtutor%error("CUSOLVERMP_FREE: error in cusolverMpDestroy: " // str(IERR) )
      ENDDO
      NCUSOLVERMPHANDLERS = 0
  END SUBROUTINE
#endif

  SUBROUTINE PDPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)

    IMPLICIT NONE

    REAL(q) :: A(*)
    INTEGER :: DESCA(*)
    INTEGER :: IA, INFO, JA, N
    CHARACTER(1) :: UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESC
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESC, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = cusolverMpPotrf_buffersize(ACC_CUSOLVERMP_HANDLE, CHAR_TO_FILL(UPLO), N, A, IA, JA, ACC_CUSOLVERMP_MATDESC, &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpPotrf_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A,devwork,devinfo)
    IERR = cusolverMpPotrf(ACC_CUSOLVERMP_HANDLE, CHAR_TO_FILL(UPLO), N, A, IA, JA, ACC_CUSOLVERMP_MATDESC, &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpPotrf: " // str(IERR) )

#if PGI_BEFORE_25_9
    IERR = cal_stream_sync(ACC_MP_COMM,acc_get_cuda_stream(ACC_ASYNC_Q))
#endif

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpPotrf devInfo: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESC)
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PDPOTRF_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PDPOTRF_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDPOTRF_OFFLOAD


  SUBROUTINE PZPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)

    IMPLICIT NONE

    COMPLEX(q) :: A(*)
    INTEGER :: DESCA(*)
    INTEGER :: IA, INFO, JA, N
    CHARACTER(1) :: UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESC
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESC, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = cusolverMpPotrf_buffersize(ACC_CUSOLVERMP_HANDLE, CHAR_TO_FILL(UPLO), N, A, IA, JA, ACC_CUSOLVERMP_MATDESC, &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpPotrf_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A,devwork,devinfo)
    IERR = cusolverMpPotrf(ACC_CUSOLVERMP_HANDLE, CHAR_TO_FILL(UPLO), N, A, IA, JA, ACC_CUSOLVERMP_MATDESC, &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpPotrf: " // str(IERR) )

#if PGI_BEFORE_25_9
    IERR = cal_stream_sync(ACC_MP_COMM,acc_get_cuda_stream(ACC_ASYNC_Q))
#endif

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpPotrf devInfo: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESC)
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PZPOTRF_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PZPOTRF_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZPOTRF_OFFLOAD


  SUBROUTINE PSSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL, IU, ABSTOL, M, NZ, W,&
                ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

    IMPLICIT NONE

    REAL(qs) :: A(*), GAP(*), W(*), WORK(*), Z(*)
    REAL(qs) :: ABSTOL, ORFAC, VL, VU
    INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
    INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
    CHARACTER(1) :: JOBZ, RANGE, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LIWORK == -1 ) THEN
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVX_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PSSYEVX_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PSSYEVX_OFFLOAD


  SUBROUTINE PDSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL, IU, ABSTOL, M, NZ, W,&
                ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

    IMPLICIT NONE

    REAL(q) :: A(*), GAP(*), W(*), WORK(*), Z(*)
    REAL(q) :: ABSTOL, ORFAC, VL, VU
    INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
    INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
    CHARACTER(1) :: JOBZ, RANGE, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LIWORK == -1 ) THEN
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVX_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PDSYEVX_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDSYEVX_OFFLOAD


  SUBROUTINE PCHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL, IU, ABSTOL, M, NZ, W,&
                ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

    IMPLICIT NONE

    COMPLEX(qs) :: A(*), WORK(*), Z(*)
    REAL(qs) :: GAP(*), RWORK(*), W(*)
    REAL(qs) :: ABSTOL, ORFAC, VL, VU
    INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
    INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
    CHARACTER(1) :: JOBZ, RANGE, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1 .OR. LIWORK == -1 ) THEN
       RWORK(1) = 1
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVX_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PCHEEVX_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PCHEEVX_OFFLOAD


  SUBROUTINE PZHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL, IU, ABSTOL, M, NZ, W,&
                ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

    IMPLICIT NONE

    COMPLEX(q) :: A(*), WORK(*), Z(*)
    REAL(q) :: GAP(*), RWORK(*), W(*)
    REAL(q) :: ABSTOL, ORFAC, VL, VU
    INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
    INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
    CHARACTER(1) :: JOBZ, RANGE, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1 .OR. LIWORK == -1 ) THEN
       RWORK(1) = 1
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVX_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PZHEEVX_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZHEEVX_OFFLOAD


  SUBROUTINE PSSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK,  INFO)

    IMPLICIT NONE

    REAL(qs) :: A(*), W(*), WORK(*), Z(*)
    INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LIWORK == -1 ) THEN
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEVD_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PSSYEVD_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PSSYEVD_OFFLOAD


  SUBROUTINE PDSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK,  INFO)

    IMPLICIT NONE

    REAL(q) :: A(*), W(*), WORK(*), Z(*)
    INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LIWORK == -1 ) THEN
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEVD_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PDSYEVD_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDSYEVD_OFFLOAD


  SUBROUTINE PCHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK,  INFO)

    IMPLICIT NONE

    COMPLEX(qs) :: A(*), WORK(*), Z(*)
    REAL(qs) :: RWORK(*), W(*)
    INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1 .OR. LIWORK == -1 ) THEN
       RWORK(1) = 1
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEVD_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PCHEEVD_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PCHEEVD_OFFLOAD


  SUBROUTINE PZHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK,  INFO)

    IMPLICIT NONE

    COMPLEX(q) :: A(*), WORK(*), Z(*)
    REAL(q) :: RWORK(*), W(*)
    INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1 .OR. LIWORK == -1 ) THEN
       RWORK(1) = 1
       IWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEVD_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PZHEEVD_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZHEEVD_OFFLOAD


  SUBROUTINE PSSYEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, INFO)

    IMPLICIT NONE

    REAL(qs) :: A(*), WORK(*), Z(*)
    REAL(qs) :: W(*)
    INTEGER :: DESCA(*), DESCZ(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, N, LWORK
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1) THEN
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PSSYEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PSSYEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PSSYEV_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PSSYEV_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PSSYEV_OFFLOAD


  SUBROUTINE PDSYEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, INFO)

    IMPLICIT NONE

    REAL(q) :: A(*), WORK(*), Z(*)
    REAL(q) :: W(*)
    INTEGER :: DESCA(*), DESCZ(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, N, LWORK
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_R_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_R_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1) THEN
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PDSYEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PDSYEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PDSYEV_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PDSYEV_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PDSYEV_OFFLOAD


  SUBROUTINE PCHEEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO)

    IMPLICIT NONE

    COMPLEX(qs) :: A(*), WORK(*), Z(*)
    REAL(qs) :: RWORK(*), W(*)
    INTEGER :: DESCA(*), DESCZ(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_32F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_32F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1) THEN
       RWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PCHEEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PCHEEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PCHEEV_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PCHEEV_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PCHEEV_OFFLOAD


  SUBROUTINE PZHEEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO)

    IMPLICIT NONE

    COMPLEX(q) :: A(*), WORK(*), Z(*)
    REAL(q) :: RWORK(*), W(*)
    INTEGER :: DESCA(*), DESCZ(*)
    INTEGER :: IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
    CHARACTER(1) :: JOBZ, UPLO

#ifdef CUSOLVERMP
    ! local variables
    TYPE(cudaDataType)               :: compute_type = cudaDataType(CUDA_C_64F)
    TYPE(cudaDataType)               :: data_type = cudaDataType(CUDA_C_64F)
    INTEGER(1), ALLOCATABLE          :: devwork(:), hostwork(:)
    INTEGER(8)                       :: devworksize, hostworksize
    INTEGER(4), ALLOCATABLE          :: devinfo(:)
    INTEGER(4)                       :: IERR
    INTEGER                          :: NPROW, NPCOL, MYPROW, MYPCOL
    TYPE(cusolverMpHandle)           :: ACC_CUSOLVERMP_HANDLE
    TYPE(cusolverMpGrid)             :: ACC_CUSOLVERMP_GRID
    TYPE(mpcomm)                     :: ACC_MP_COMM
    TYPE(cusolverMpMatrixDescriptor) :: ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_MATDESCZ
    INTEGER(4)                       :: MAPPING = CUSOLVERMP_GRID_MAPPING_ROW_MAJOR

    IF ( LWORK == -1 .OR. LRWORK == -1) THEN
       RWORK(1) = 1
       WORK(1) = 1
       INFO = 0
       RETURN
    ENDIF

    ACC_MP_COMM = GET_MP_COMM(DESCA(CTXT_))

    ACC_CUSOLVERMP_HANDLE = GET_CUSOLVERMP_HANDLER(acc_get_device_num(ACC_DEVICE_NVIDIA), ACC_ASYNC_Q)

    CALL BLACS_GRIDINFO(DESCA(CTXT_), NPROW, NPCOL, MYPROW, MYPCOL)

!note(sm): VASP actually uses ROW MAJOR for ACFDT, but it's always a 1D grid.
    IERR = cusolverMpCreateDeviceGrid(ACC_CUSOLVERMP_HANDLE, ACC_CUSOLVERMP_GRID, ACC_MP_COMM, &
              NPROW, NPCOL, MAPPING)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpCreateDeviceGrid: " // str(IERR) )

    IF (DESCA(RSRC_) /= 0 .OR. DESCA(CSRC_) /=0 ) CALL vtutor%error("PZHEEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCA(RSRC_)) // str(DESCA(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCA, ACC_CUSOLVERMP_GRID, data_type, DESCA(M_), DESCA(N_), &
              DESCA(MB_), DESCA(NB_), DESCA(RSRC_), DESCA(CSRC_), DESCA(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

    IF (DESCZ(RSRC_) /= 0 .OR. DESCZ(CSRC_) /=0 ) CALL vtutor%error("PZHEEV_OFFLOAD: rsrc or csrc not 0: " // &
                                                     str(DESCZ(RSRC_)) // str(DESCZ(CSRC_)) )

    IERR = cusolverMpCreateMatrixDesc(ACC_CUSOLVERMP_MATDESCZ, ACC_CUSOLVERMP_GRID, data_type, DESCZ(M_), DESCZ(N_), &
              DESCZ(MB_), DESCZ(NB_), DESCZ(RSRC_), DESCZ(CSRC_), DESCZ(LLD_))
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpCreateMatrixDesc: " // str(IERR) )

!$ACC HOST_DATA USE_DEVICE(A, W, Z)
    IERR = cusolverMpSyevd_buffersize(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devworksize, hostworksize)
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpSyevd_buffersize: " // str(IERR) )

    ALLOCATE(devwork(devworksize), hostwork(hostworksize), devinfo(1))
    devinfo(1)=0
!$ACC ENTER DATA CREATE(devwork) COPYIN(devinfo)

!$ACC HOST_DATA USE_DEVICE(A, W, Z, devwork, devinfo)
    IERR = cusolverMpSyevd(ACC_CUSOLVERMP_HANDLE, JOBZ, CHAR_TO_FILL(UPLO), N, A, IA, JA, &
              ACC_CUSOLVERMP_MATDESCA, W, Z, IZ, JZ, ACC_CUSOLVERMP_MATDESCZ,  &
              compute_type, devwork, devworksize, hostwork, hostworksize, devinfo(1))
!$ACC END HOST_DATA
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpSyevd: " // str(IERR) )

!$ACC EXIT DATA DELETE(devwork) COPYOUT(devinfo) WAIT
    IF ( devinfo(1) /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpSyevd devInfo: " // str(IERR) )

    DEALLOCATE(devwork, hostwork, devinfo)

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCA)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyMatrixDesc(ACC_CUSOLVERMP_MATDESCZ)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpDestroyMatrixDesc: " // str(IERR) )

    IERR = cusolverMpDestroyGrid(ACC_CUSOLVERMP_GRID)
    IF ( IERR /= 0 ) CALL vtutor%error("PZHEEV_OFFLOAD: error in cusolverMpDeviceGrid: " // str(IERR) )

    INFO=0
#else
    CALL vtutor%bug("PZHEEV_OFFLOAD: CUSOLVERMP not defined", __FILE__, __LINE__)
#endif

  END SUBROUTINE PZHEEV_OFFLOAD

  END MODULE mcusolvermp


#ifdef CUDA_PINNED_MEMORY
  MODULE mcuda_pinned_memory
  USE prec
  USE iso_c_binding, ONLY : C_PTR,C_SIZE_T,C_LOC

  USE cudafor, ONLY : cudaHostRegister, cudaHostRegisterPortable, &
                      cudaHostUnregister, cudaPointerGetAttributes, cudaPointerAttributes

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PIN_MEMORY,UNPIN_MEMORY,C_LOC,C_SIZE_T

  CONTAINS

  SUBROUTINE PIN_MEMORY(CPTR,ISIZE,ISTATUS)
    TYPE(C_PTR), VALUE :: CPTR
    INTEGER(KIND=C_SIZE_T) :: ISIZE
    INTEGER, OPTIONAL :: ISTATUS
    ! local variables
    INTEGER :: IERR

    IERR = cudaHostRegister(CPTR,ISIZE,cudaHostRegisterPortable)

    IF (PRESENT(ISTATUS)) ISTATUS=IERR
    RETURN
  END SUBROUTINE PIN_MEMORY


  SUBROUTINE UNPIN_MEMORY(CPTR,ISTATUS)
    TYPE(C_PTR), VALUE :: CPTR
    INTEGER, OPTIONAL :: ISTATUS
    ! local variables
    TYPE(cudaPointerAttributes) :: ATTR
    INTEGER :: IERR

    IERR = cudaPointerGetAttributes(ATTR,CPTR)

    IF (IERR==0) IERR = cudaHostUnregister(CPTR)

    IF (PRESENT(ISTATUS)) ISTATUS=IERR
    RETURN
  END SUBROUTINE UNPIN_MEMORY

  END MODULE mcuda_pinned_memory
#endif // CUDA_PINNED_MEMORY

#ifdef USENCCL
  SUBROUTINE nccl_init_comm(COMM)

    USE moffload_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclResult, ncclUniqueId, ncclInvalidUsage, ncclSuccess, &
                     ncclGetUniqueId, ncclCommInitRank
    USE tutor, ONLY: vtutor

    IMPLICIT NONE

    TYPE(communic) COMM

    ! local variables
    TYPE(ncclResult)   ncclRes
    TYPE(ncclUniqueId) ncclUID
    INTEGER ierror,i

    IF ( .NOT. LUSENCCL ) RETURN

    IF (COMM%NODE_ME==1) ncclRes = ncclGetUniqueId(ncclUID)

    CALL M_bcast_c_char(COMM, ncclUID%internal, 128)

    ncclRes = ncclCommInitRank( COMM%NCCL_COMM, COMM%NCPU, ncclUID, COMM%NODE_ME-1)

    IF ( ncclRes == ncclInvalidUsage ) THEN
       COMM%LUSENCCL = .FALSE.
    ELSE
       IF ( ncclRes /= ncclSuccess ) &
          CALL vtutor%bug('nccl_init_comm: Error in ncclCommInitRank', __FILE__, __LINE__)
       COMM%LUSENCCL = .TRUE.
    ENDIF

    CALL M_barrier( COMM )

  END SUBROUTINE nccl_init_comm


  SUBROUTINE nccl_group_start()
    USE nccl, ONLY : ncclResult, ncclGroupStart
    ! local variables
    TYPE(ncclResult) :: ncclRes
    ncclRes = ncclGroupStart()
  END SUBROUTINE nccl_group_start


  SUBROUTINE nccl_group_end()
    USE nccl, ONLY : ncclResult, ncclGroupEnd
    ! local variables
    TYPE(ncclResult) :: ncclRes
    ncclRes = ncclGroupEnd()
  END SUBROUTINE nccl_group_end


  SUBROUTINE nccl_bcast_s(vec, counts, root, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclBcast, ncclFloat
    TYPE(communic) :: comm
    REAL(qs) :: vec(counts)
    INTEGER :: counts, root, stream
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec)
    ncclRes = ncclBcast(c_devloc(vec), counts, ncclFloat, root, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_bcast_s


  SUBROUTINE nccl_bcast_d(vec, counts, root, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclBcast, ncclDouble
    TYPE(communic) :: comm
    REAL(q) :: vec(counts)
    INTEGER :: counts, root, stream
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec)
    ncclRes = ncclBcast(c_devloc(vec), counts, ncclDouble, root, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_bcast_d


  SUBROUTINE nccl_reduce_d(vec1, vec2, counts, root, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclReduce, ncclDouble, ncclSum
    TYPE(communic) :: comm
    REAL(q) :: vec1(counts), vec2(counts)
    INTEGER :: counts, root, stream
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclDouble, ncclSum, root, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_reduce_d


  SUBROUTINE nccl_send_d(vec, counts, peer, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclSend, ncclDouble
    TYPE(communic) :: comm
    INTEGER :: counts, peer, stream
    REAL(q) :: vec(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec)
    ncclRes = ncclSend(c_devloc(vec), counts, ncclDouble, peer, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_send_d


  SUBROUTINE nccl_recv_d(vec, counts, peer, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclRecv, ncclDouble
    TYPE(communic) :: comm
    INTEGER :: counts, peer, stream
    REAL(q) :: vec(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec)
    ncclRes = ncclRecv(c_devloc(vec), counts, ncclDouble, peer, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_recv_d


  SUBROUTINE nccl_allgather_d(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllGather, ncclDouble
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    REAL(q) :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllGather(c_devloc(vec1), c_devloc(vec2), counts, ncclDouble, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE nccl_allgather_d


  SUBROUTINE nccl_min_i(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclInt, ncclMin
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    INTEGER :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclInt, ncclMin, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_min_i


  SUBROUTINE nccl_min_i8(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclInt64, ncclMin
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    INTEGER(qi8) :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclInt64, ncclMin, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_min_i8


  SUBROUTINE nccl_max_i(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclInt, ncclMax
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    INTEGER :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclInt, ncclMax, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_max_i


  SUBROUTINE nccl_max_i8(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclInt64, ncclMax
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    INTEGER(qi8) :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclInt64, ncclMax, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_max_i8


  SUBROUTINE nccl_sum_s(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclFloat, ncclSum
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    REAL(qs) :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclFloat, ncclSum, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_sum_s


  SUBROUTINE nccl_sum_d(vec1, vec2, counts, comm, stream)
    USE prec
    USE mopenacc_struct_def
    USE mpimy, ONLY : communic
    USE nccl, ONLY : ncclComm, ncclResult, ncclAllReduce, ncclDouble, ncclSum
    TYPE(communic) :: comm
    INTEGER :: counts, stream
    REAL(q) :: vec1(counts), vec2(counts)
    ! local variables
    TYPE(ncclResult) :: ncclRes
!$ACC HOST_DATA USE_DEVICE(vec1,vec2)
    ncclRes = ncclAllReduce(c_devloc(vec1), c_devloc(vec2), counts, ncclDouble, ncclSum, comm%nccl_comm, acc_get_cuda_stream(stream))
!$ACC END HOST_DATA
  END SUBROUTINE  nccl_sum_d

#endif //USENCCL

#endif // ACC_OFFLOAD
#endif // NVCUDA
