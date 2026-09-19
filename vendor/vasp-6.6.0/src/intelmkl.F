#include "symbol.inc"
#define divide_and_conquer
#ifdef INTELMKL
#ifdef OMP_OFFLOAD
!***********************************************************************
!
!***********************************************************************
#ifdef MPI
    MODULE mdevice_info

      USE main_mpi
      USE iso_c_binding

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: GET_NUM_DEVICES_GLOBAL, GET_NUM_RANKS_ON_DEVICE, GET_FREE_MEM_DEVICE

      INTERFACE
         SUBROUTINE ze_device_get_uuid(id, device_num, info) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr),    VALUE :: id
            INTEGER(c_int), VALUE :: device_num
            INTEGER(c_int)        :: info
         END SUBROUTINE ze_device_get_uuid
      END INTERFACE

    CONTAINS

!****************** FUNCTION GET_NUM_DEVICES_GLOBAL ********************
!
!***********************************************************************

      INTEGER FUNCTION GET_NUM_DEVICES_GLOBAL(MY_DEVICE_NUM)
      INTEGER :: MY_DEVICE_NUM
      ! local variables
      CHARACTER(KIND=c_char), TARGET :: uuidchar(16)
      INTEGER :: uuidint(16)
      INTEGER, ALLOCATABLE :: uuidints(:,:)
      INTEGER, ALLOCATABLE :: uuidunique(:,:)
      INTEGER :: totalGPUs, I, J, INFO
      LOGICAL :: isUnique

      CALL ze_device_get_uuid(c_loc(uuidchar), MY_DEVICE_NUM, INFO)

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

      INTEGER FUNCTION GET_NUM_RANKS_ON_DEVICE(MY_DEVICE_NUM)
      INTEGER :: MY_DEVICE_NUM
      ! local variables
      CHARACTER(KIND=c_char), TARGET :: uuidchar(16)
      INTEGER :: uuidint(16)
      INTEGER, ALLOCATABLE :: uuidints(:,:)
      INTEGER :: numranks, I, INFO

      CALL ze_device_get_uuid(c_loc(uuidchar), MY_DEVICE_NUM, INFO)
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

      REAL(q) :: MFREE
      INTEGER :: IERR

      IERR = -1
      MFREE = 0

      END FUNCTION GET_FREE_MEM_DEVICE

    END MODULE mdevice_info
#endif // MPI


!***********************************************************************
!
!***********************************************************************
    include "mkl_dfti_omp_offload.f90"

    MODULE mmklfft

      USE prec
      USE moffload_struct_def, ONLY : OFFLOAD_ALLOW
      USE mgrid_struct_def, ONLY : grid_3d

      USE tutor
      USE string

      USE iso_c_binding

      USE mkl_dfti_omp_offload

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: FFT_OFFLOAD_MAKEPLANS, FFT_OFFLOAD_MAKEPLAN, FFT_OFFLOAD_DESTROYPLAN, FFT_OFFLOAD_GETPLAN
      PUBLIC :: fft_offload_plan

      PUBLIC :: FFTBAS_PLAN_OFFLOAD, FFTBRC_PLAN_OFFLOAD
      PUBLIC :: FFTBAS_PLAN_MU_OFFLOAD, FFTBRC_PLAN_MU_OFFLOAD
      PUBLIC :: FFTBAS_PLAN_MPI_OFFLOAD, FFTBRC_PLAN_MPI_OFFLOAD

      TYPE plan_type
         TYPE(dfti_descriptor), POINTER :: desc_fwd
         TYPE(dfti_descriptor), POINTER :: desc_bck
         INTEGER :: fwd_stride(4), bkw_stride(4)
         INTEGER :: fwd_dist, bkw_dist
         INTEGER :: Asize,Bsize
      END TYPE plan_type

      TYPE fft_offload_plan
         TYPE(grid_3d), POINTER :: GRID
         INTEGER :: stream
         INTEGER :: nbatch
         INTEGER :: nstride
         CHARACTER(LEN=:), ALLOCATABLE :: tag
         ! plans for serial FFTs (3D)
         TYPE(plan_type) :: FFT_3D_Z2Z
         TYPE(plan_type) :: FFT_3D_D2Z
         ! plans for MPI distributed FFTs (1D)
         TYPE(plan_type) :: FFT_RC_Z2Z
         TYPE(plan_type) :: FFT_IN_Z2Z
         TYPE(plan_type) :: FFT_RL_Z2Z
         TYPE(plan_type) :: FFT_RL_D2Z
      END TYPE fft_offload_plan

      TYPE(fft_offload_plan), ALLOCATABLE :: plans(:)

      INTEGER :: createdPlans = 0

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE planAssign
      END INTERFACE

      INTERFACE OPERATOR(==)
         MODULE PROCEDURE planCompare
      END INTERFACE

    CONTAINS

!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLANS *******************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_MAKEPLANS(GRID, NSTREAMS)
      TYPE(grid_3d), INTENT(IN), TARGET :: GRID
      INTEGER, INTENT(IN), OPTIONAL :: NSTREAMS

      PROFILING_START('fft_offload_makeplans')

      IF (.NOT.OFFLOAD_ALLOW) THEN
         PROFILING_STOP('fft_offload_makeplans')
         RETURN
      ENDIF

      CALL FFT_OFFLOAD_MAKEPLAN(GRID,-1)

      PROFILING_STOP('fft_offload_makeplans')

      END SUBROUTINE FFT_OFFLOAD_MAKEPLANS


!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLAN ********************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_MAKEPLAN(GRID, ISTREAM, NBATCH, NSTRIDE, TAG)
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

      NEW%GRID => GRID

      NEW%STREAM = -1

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
      TYPE(fft_offload_plan), TARGET :: PLAN
      ! local variables
      TYPE(plan_type), POINTER :: ptype
      INTEGER :: NG(3), BATCHES(3), IERR, IERR_, LD

      NG = PLAN%GRID%NGPTAR

      IERR = 0

      ! serial FFT (3D)
      IF (PLAN%GRID%RL%NFAST==1 .AND. PLAN%GRID%RL_FFT%NFAST==1) THEN

         ! single FFT
         IF (PLAN%NBATCH==1) THEN
            ! real-2-complex
            IF (PLAN%GRID%REAL2CPLX) THEN
               ! FFT_3D_D2Z: 3d real-2-complex/complex-2-real
               ptype => PLAN%FFT_3D_D2Z

               ptype%fwd_stride = [0, 1, (NG(1)+2), NG(2)*(NG(1)+2)]
               ptype%bkw_stride = [0, 1, (NG(1)+2)/2, NG(2)*(NG(1)+2)/2]
               ptype%fwd_dist = 1
               ptype%bkw_dist = 1
               ptype%Asize = (NG(1)+2)*NG(2)*NG(3)
               ptype%Bsize = (NG(1)+2)/2*NG(2)*NG(3)

               IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_REAL, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, ptype%fwd_stride);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, ptype%bkw_stride);
               IERR = IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

               !IERR =IERR + DftiSetValue(ptype%desc, DFTI_INPUT_STRIDES, [0, 1, NG(1), NG(1)*NG(2)]);
               !IERR =IERR + DftiSetValue(ptype%desc, DFTI_OUTPUT_STRIDES, [0, 1, INT(NG(1)/2.0)+1, NG(2)*(INT(NG(1)/2.0)+1)] );

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               IERR = IERR + DftiCreateDescriptor(ptype%desc_bck, DFTI_DOUBLE, DFTI_REAL, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_INPUT_STRIDES, ptype%bkw_stride);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_OUTPUT_STRIDES, ptype%fwd_stride);
               IERR = IERR + DftiSetValue( ptype%desc_bck, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_bck)
               IERR = IERR + IERR_

            ! complex-2-complex
            ELSE
               ! FFT_3D_Z2Z: 3d complex-2-complex
               ptype => PLAN%FFT_3D_Z2Z

               ptype%fwd_stride = [0, 1, NG(1), NG(1)*NG(2)]
               ptype%bkw_stride = [0, 1, NG(1), NG(1)*NG(2)]
               ptype%fwd_dist = 1
               ptype%bkw_dist = 1
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = NG(1)*NG(2)*NG(3)

               IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_COMPLEX, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)
               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               !IERR =IERR + DftiCreateDescriptor(ptype%desc_bck, DFTI_DOUBLE, DFTI_COMPLEX, 3, [NG(1),NG(2),NG(3)]);
               !IERR =IERR + DftiSetValue( ptype%desc_bck, DFTI_PLACEMENT, DFTI_INPLACE)

               !!$OMP DISPATCH
               !IERR = IERR + DftiCommitDescriptor(ptype%desc_bck)
               ptype%desc_bck => ptype%desc_fwd
            ENDIF
         ! batched FFT
         ELSE
!!            CALL vtutor%bug("MAKEPLAN: batched serial FFTs not implemented yet", __FILE__, __LINE__)
            ! write(*,*) 'MAKEPLAN: batched serial FFTs starting, NBATCH=',PLAN%NBATCH
            LD = PLAN%GRID%MPLWV; IF (PLAN%NSTRIDE > 0) LD = PLAN%NSTRIDE

            IF (PLAN%GRID%REAL2CPLX) THEN
               ! write(*,*) 'MAKEPLAN: batched serial FFTs real-2-complex'
               ptype => PLAN%FFT_3D_D2Z
               ptype%fwd_stride = [0, 1, 2*(NG(1)/2+1), NG(2)*2*(NG(1)/2+1)]
               ptype%bkw_stride = [0, 1, (NG(1)/2+1), NG(2)*(NG(1)/2+1)]
               ptype%fwd_dist = 2*LD
               ptype%bkw_dist = LD
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = (NG(1)+2)/2*NG(2)*NG(3)

               IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_REAL, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, PLAN%NBATCH );
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, ptype%fwd_stride);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, ptype%bkw_stride);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, ptype%fwd_dist);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, ptype%bkw_dist);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               IERR = IERR + DftiCreateDescriptor(ptype%desc_bck, DFTI_DOUBLE, DFTI_REAL, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_NUMBER_OF_TRANSFORMS, PLAN%NBATCH );
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_INPUT_STRIDES, ptype%bkw_stride);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_OUTPUT_STRIDES, ptype%fwd_stride);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_INPUT_DISTANCE, ptype%bkw_dist);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_OUTPUT_DISTANCE, ptype%fwd_dist);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_bck)
               IERR = IERR + IERR_
            ELSE
               ! write(*,*) 'MAKEPLAN: batched serial FFTs complex-2-complex'
               ptype => PLAN%FFT_3D_Z2Z
               ptype%fwd_stride = [0, 1, NG(1), NG(1)*NG(2)]
               ptype%bkw_stride = [0, 1, NG(1), NG(1)*NG(2)]
               ptype%fwd_dist = LD
               ptype%bkw_dist = LD
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = NG(1)*NG(2)*NG(3)

               IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_COMPLEX, 3, [NG(1),NG(2),NG(3)]);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, PLAN%NBATCH );
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, ptype%fwd_stride);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, ptype%bkw_stride);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, ptype%fwd_dist);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, ptype%bkw_dist);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, NG(1)*NG(2)*NG(3));
               IERR = IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               ptype%desc_bck => ptype%desc_fwd
            ENDIF
         ENDIF

      ! parallel FFT
      ELSE
         IF (PLAN%NBATCH/=1) THEN
            CALL vtutor%bug("MAKEPLAN: batched parallel ffts not supported yet",__FILE__,__LINE__)
         ENDIF

         BATCHES(1)=PLAN%GRID%RC%NCOL
         BATCHES(2)=PLAN%GRID%IN%NCOL
         BATCHES(3)=PLAN%GRID%RL_FFT%NCOL

         ! FFT_RC_Z2Z: 1d complex-2-complex FFTs on data with reciprocal space layout
         IF (BATCHES(1)>0) THEN
            ptype => PLAN%FFT_RC_Z2Z

            ptype%fwd_stride = [0, 1, 1, 1]
            ptype%bkw_stride = [0, 1, 1, 1]
            ptype%fwd_dist = NG(1)
            ptype%bkw_dist = NG(1)
            ptype%Asize = NG(1)*BATCHES(1)
            ptype%Bsize = NG(1)*BATCHES(1)

            IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_COMPLEX, 1, NG(1));
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, BATCHES(1) );
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, [0,1]);
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, [0,1]);
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, NG(1));
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, NG(1));
            IERR = IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

            !$OMP DISPATCH
            IERR_= DftiCommitDescriptor(ptype%desc_fwd)
            IERR = IERR + IERR_

            ptype%desc_bck => ptype%desc_fwd
         ENDIF
         ! FFT_IN_Z2Z: 1d complex-2-complex FFTs on data with intermediate layout
         IF (BATCHES(2)>0) THEN
            ptype => PLAN%FFT_IN_Z2Z

            ptype%fwd_stride = [0, BATCHES(2), 1, 1]
            ptype%bkw_stride = [0, BATCHES(2), 1, 1]
            ptype%fwd_dist = 1
            ptype%bkw_dist = 1
            ptype%Asize = NG(2)*BATCHES(2)
            ptype%Bsize = NG(2)*BATCHES(2)

            IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_COMPLEX, 1, NG(2));
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, BATCHES(2) );
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, [0, BATCHES(2)]);
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, [0, BATCHES(2)]);
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, 1)
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, 1)
            IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

            !$OMP DISPATCH
            IERR_= DftiCommitDescriptor(ptype%desc_fwd)
            IERR = IERR + IERR_

            ptype%desc_bck => ptype%desc_fwd
         ENDIF
         ! 1d FFTs on data with real space layout
         IF (BATCHES(3)>0) THEN
            ! FFT_RL_D2Z: real-2-complex/complex-2-real
            IF (PLAN%GRID%REAL2CPLX) THEN
               ptype => PLAN%FFT_RL_D2Z

               ptype%fwd_stride = [0, 1, 1, 1]
               ptype%bkw_stride = [0, 1, 1, 1]
               ptype%fwd_dist = NG(3)+2
               ptype%bkw_dist = (NG(3)+2)/2
               ptype%Asize = (NG(3)+2)*BATCHES(3)
               ptype%Bsize = (NG(3)+2)/2*BATCHES(3)

               IERR = IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_REAL, 1, NG(3));
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, BATCHES(3) );
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, [0, 1]);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, [0, 1]);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, NG(3)+2);
               IERR = IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, (NG(3)+2)/2);
               IERR = IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               IERR = IERR + DftiCreateDescriptor(ptype%desc_bck, DFTI_DOUBLE, DFTI_REAL, 1, NG(3));
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_NUMBER_OF_TRANSFORMS, BATCHES(3) );
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_INPUT_STRIDES, [0, 1]);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_OUTPUT_STRIDES, [0, 1]);
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_INPUT_DISTANCE, (NG(3)+2)/2 );
               IERR = IERR + DftiSetValue(ptype%desc_bck, DFTI_OUTPUT_DISTANCE, NG(3)+2 );
               IERR = IERR + DftiSetValue( ptype%desc_bck, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_bck)
               IERR = IERR + IERR_

            ! FFT_RL_Z2Z: complex-2-complex
            ELSE
               ptype => PLAN%FFT_RL_Z2Z

               ptype%fwd_stride = [0, 1, 1, 1]
               ptype%bkw_stride = [0, 1, 1, 1]
               ptype%fwd_dist = NG(3)
               ptype%bkw_dist = NG(3)
               ptype%Asize = NG(3)*batches(3)
               ptype%Bsize = NG(3)*batches(3)

               IERR =IERR + DftiCreateDescriptor(ptype%desc_fwd, DFTI_DOUBLE, DFTI_COMPLEX, 1, NG(3));
               IERR =IERR + DftiSetValue(ptype%desc_fwd, DFTI_NUMBER_OF_TRANSFORMS, BATCHES(3) );
               IERR =IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_STRIDES, [0, 1]);
               IERR =IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_STRIDES, [0, 1]);
               IERR =IERR + DftiSetValue(ptype%desc_fwd, DFTI_INPUT_DISTANCE, NG(3));
               IERR =IERR + DftiSetValue(ptype%desc_fwd, DFTI_OUTPUT_DISTANCE, NG(3));
               IERR =IERR + DftiSetValue( ptype%desc_fwd, DFTI_PLACEMENT, DFTI_INPLACE)

               !$OMP DISPATCH
               IERR_= DftiCommitDescriptor(ptype%desc_fwd)
               IERR = IERR + IERR_

               ptype%desc_bck => ptype%desc_fwd
            ENDIF
         ENDIF
      ENDIF

      ! In case something went wrong along the way ...
      IF (IERR/=0) THEN
         CALL vtutor%error("MAKEPLAN: could not create plan")
      ENDIF
      END SUBROUTINE MAKEPLAN


!****************** SUBROUTINE FFT_OFFLOAD_DESTROYPLAN *****************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_DESTROYPLAN(GRID,ISTREAM,NBATCH,TAG)
      TYPE(grid_3d), INTENT(IN), TARGET     :: GRID
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

         IF (.NOT.ASSOCIATED(plans(i)%GRID,GRID)) keepPLan = .TRUE.

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
      INTEGER :: IERR

      IERR = 0

      IF (PLAN%GRID%RL%NFAST==1) THEN
         ! serial fft
         IF (PLAN%GRID%REAL2CPLX) THEN
            IERR = IERR + DftiFreeDescriptor(PLAN%FFT_3D_D2Z%desc_fwd)
            IERR = IERR + DftiFreeDescriptor(PLAN%FFT_3D_D2Z%desc_bck)
         ELSE
            IERR = IERR + DftiFreeDescriptor(PLAN%FFT_3D_Z2Z%desc_fwd)
         ENDIF
      ELSE
         IF (PLAN%GRID%RC%NCOL>0) IERR = IERR + DftiFreeDescriptor(PLAN%FFT_RC_Z2Z%desc_fwd)
         IF (PLAN%GRID%IN%NCOL>0) IERR = IERR + DftiFreeDescriptor(PLAN%FFT_IN_Z2Z%desc_fwd)
         IF (PLAN%GRID%RL_FFT%NCOL>0) THEN
            IF (PLAN%GRID%REAL2CPLX) THEN
               IERR = IERR + DftiFreeDescriptor(PLAN%FFT_RL_D2Z%desc_fwd)
               IERR = IERR + DftiFreeDescriptor(PLAN%FFT_RL_D2Z%desc_bck)
            ELSE
               IERR = IERR + DftiFreeDescriptor(PLAN%FFT_RL_Z2Z%desc_fwd)
            ENDIF
         ENDIF
      ENDIF
      ! In case something went wrong along the way ...
      IF (IERR/=0) THEN
         CALL vtutor%error("DESTROYPLAN: could not destroy plan")
      ENDIF
      END SUBROUTINE DESTROYPLAN


!****************** SUBROUTINE FFT_OFFLOAD_GETPLAN *********************
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

      FUNCTION planAlloc()
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

      a%FFT_RC_Z2Z    = b%FFT_RC_Z2Z
      a%FFT_IN_Z2Z    = b%FFT_IN_Z2Z
      a%FFT_RL_Z2Z    = b%FFT_RL_Z2Z
      a%FFT_RL_D2Z    = b%FFT_RL_D2Z
      END SUBROUTINE planAssign


!****************** FUNCTION PLANCOMPARE *******************************
!
!***********************************************************************

      FUNCTION planCompare(a,b) RESULT(res)
      TYPE(fft_offload_plan), INTENT(IN) :: a,b
      LOGICAL :: res

      res=.FALSE.

      IF (.NOT.ASSOCIATED(a%GRID,b%GRID)) RETURN
      IF (a%stream/=b%stream) RETURN
      IF (a%nbatch/=b%nbatch) RETURN

      res=.TRUE.
      END FUNCTION planCompare


!************************* SUBROUTINE FFTBAS_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbas_plan_offload')

      IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISIGN<=0) THEN
         !$OMP DISPATCH IS_DEVICE_PTR(C)
         IERR = DftiComputeForward(PLAN%FFT_3D_Z2Z%desc_fwd,C(1))
      ELSE
         !$OMP DISPATCH IS_DEVICE_PTR(C)
         IERR = DftiComputeBackward(PLAN%FFT_3D_Z2Z%desc_bck,C(1))
      ENDIF

      PROFILING_STOP('fftbas_plan_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_OFFLOAD


!************************* SUBROUTINE FFTBRC_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q), TARGET :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      REAL(q), POINTER :: A(:)

      PROFILING_START('fftbrc_plan_offload')

      CALL c_f_pointer(c_loc(C),A,[2*GRID%MPLWV])

      IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISIGN<=0) THEN
         !$OMP DISPATCH IS_DEVICE_PTR(A)
         IERR = DftiComputeForward(PLAN%FFT_3D_D2Z%desc_fwd,A)
      ELSE
         !$OMP DISPATCH IS_DEVICE_PTR(A)
         IERR = DftiComputeBackward(PLAN%FFT_3D_D2Z%desc_bck,A)
      ENDIF

      PROFILING_STOP('fftbrc_plan_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

!!      CALL vtutor%bug("FFTBRC_PLAN_MU_OFFLOAD: not implemented yet", __FILE__, __LINE__)

      PROFILING_START('fftbas_plan_mu_offload')

      IERR=FFT_OFFLOAD_GETPLAN(GRID,-1,N,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_MU_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISN<=0) THEN
         !$OMP DISPATCH IS_DEVICE_PTR(C)
         IERR = DftiComputeForward(PLAN%FFT_3D_Z2Z%desc_fwd,C(1))
      ELSE
         !$OMP DISPATCH IS_DEVICE_PTR(C)
         IERR = DftiComputeBackward(PLAN%FFT_3D_Z2Z%desc_bck,C(1))
      ENDIF

      PROFILING_STOP('fftbas_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q), TARGET :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      REAL(q), POINTER :: A(:)

!!      CALL vtutor%bug("FFTBRC_PLAN_OFFLOAD: not implemented yet", __FILE__, __LINE__)

      PROFILING_START('fftbrc_plan_mu_offload')

      IERR=FFT_OFFLOAD_GETPLAN(GRID,-1,N,PLAN)

      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_MU_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(C),A,[2*LDC*N])

      IF (ISN<=0) THEN
         !$OMP DISPATCH IS_DEVICE_PTR(A)
         IERR = DftiComputeForward(PLAN%FFT_3D_D2Z%desc_fwd,A)
      ELSE
         !$OMP DISPATCH IS_DEVICE_PTR(A)
         IERR = DftiComputeBackward(PLAN%FFT_3D_D2Z%desc_bck,A)
      ENDIF

      PROFILING_STOP('fftbrc_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE tutor, ONLY : vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q), TARGET :: A(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR, IERR_

      COMPLEX(q), POINTER :: C(:)
      COMPLEX(q), POINTER, SAVE :: RCVBUF(:),SNDBUF(:)

      PROFILING_START('fftbas_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))
!$OMP TARGET ENTER DATA MAP(ALLOC:SNDBUF,RCVBUF)

      IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_MPI_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

!=======================================================================
! transformation (q->r)
!=======================================================================
      IF (ISIGN == 1) THEN
         IF (GRID%RC%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeBackward(PLAN%FFT_RC_Z2Z%desc_bck,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%IN%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeBackward(PLAN%FFT_IN_Z2Z%desc_bck,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%RL_FFT%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeBackward(PLAN%FFT_RL_Z2Z%desc_bck,C)
            IERR = IERR + IERR_
         ENDIF

!=======================================================================
! transformation (r->q)
!=======================================================================
      ELSE
         IF (GRID%RL_FFT%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeForward(PLAN%FFT_RL_Z2Z%desc_fwd,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%IN%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeForward(PLAN%FFT_IN_Z2Z%desc_fwd,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%RC%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeForward(PLAN%FFT_RC_Z2Z%desc_fwd,C)
            IERR = IERR + IERR_
         ENDIF
      ENDIF

!$OMP TARGET EXIT DATA MAP(RELEASE:SNDBUF,RCVBUF)
      PROFILING_STOP('fftbas_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE tutor, ONLY : vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q), TARGET :: A(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR, IERR_

      COMPLEX(q), POINTER :: C(:)
      COMPLEX(q), POINTER, SAVE :: RCVBUF(:),SNDBUF(:)

      PROFILING_START('fftbrc_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))
!$OMP TARGET ENTER DATA MAP(ALLOC:SNDBUF,RCVBUF)

      IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_MPI_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

!=======================================================================
! transformation (q->r)
!=======================================================================
      IF (ISIGN == 1) THEN
         IF (GRID%RC%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeBackward(PLAN%FFT_RC_Z2Z%desc_bck,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%IN%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeBackward(PLAN%FFT_IN_Z2Z%desc_bck,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%RL_FFT%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(A)
            IERR_= DftiComputeBackward(PLAN%FFT_RL_D2Z%desc_bck,A(1))
            IERR = IERR + IERR_
         ENDIF

!=======================================================================
! transformation (r->q)
!=======================================================================
      ELSE
         IF (GRID%RL_FFT%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(A)
            IERR_= DftiComputeForward(PLAN%FFT_RL_D2Z%desc_fwd,A(1))
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%IN%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeForward(PLAN%FFT_IN_Z2Z%desc_fwd,C)
            IERR = IERR + IERR_
         ENDIF

         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%RC%NCOL>0) THEN
            !$OMP DISPATCH IS_DEVICE_PTR(C)
            IERR_= DftiComputeForward(PLAN%FFT_RC_Z2Z%desc_fwd,C)
            IERR = IERR + IERR_
         ENDIF
      ENDIF

!$OMP TARGET EXIT DATA MAP(RELEASE:SNDBUF,RCVBUF)
      PROFILING_STOP('fftbrc_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD

    END MODULE mmklfft


!***********************************************************************
!
!***********************************************************************
    include "mkl_blas_omp_offload_lp64.f90"

    MODULE mmklblas

      USE prec
      USE moffload_struct_def, ONLY : OFFLOAD_ALLOW
      USE onemkl_blas_omp_offload_lp64, ONLY:   DCOPY, &
                                                SCOPY, &
                                                CCOPY, &
                                                ZCOPY, &
                                                DAXPY, &
                                                ZAXPY, &
                                                DSCAL, &
                                                ZSCAL, &
                                                ZDSCAL, &
                                                DDOT, &
                                                ZDOTC, &
                                                CDOTC, &
                                                DGEMV, &
                                                ZGEMV, &
                                                DGEMM, &
                                                ZGEMM, &
                                                SGEMM, &
                                                CGEMM, &
                                                DTRMM, &
                                                ZTRMM

      USE tutor

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: MKLBLAS_INIT, MKLBLAS_FINISH
      PUBLIC :: DCOPY_OFFLOAD, SCOPY_OFFLOAD, CCOPY_OFFLOAD, ZCOPY_OFFLOAD
      PUBLIC :: DAXPY_OFFLOAD, ZAXPY_OFFLOAD
      PUBLIC :: DSCAL_OFFLOAD, ZSCAL_OFFLOAD, ZDSCAL_OFFLOAD
      PUBLIC :: DDOT_OFFLOAD, ZDOTC_OFFLOAD, CDOTC_OFFLOAD
      PUBLIC :: DGEMV_OFFLOAD, ZGEMV_OFFLOAD
      PUBLIC :: DGEMM_OFFLOAD, ZGEMM_OFFLOAD, SGEMM_OFFLOAD, CGEMM_OFFLOAD
      PUBLIC :: DTRMM_OFFLOAD, ZTRMM_OFFLOAD
      PUBLIC :: DLACPY_OFFLOAD

    CONTAINS

      SUBROUTINE MKLBLAS_INIT()
         IF (.NOT.OFFLOAD_ALLOW) RETURN
      END SUBROUTINE MKLBLAS_INIT


      SUBROUTINE MKLBLAS_FINISH()
         IF (.NOT.OFFLOAD_ALLOW) RETURN
      END SUBROUTINE MKLBLAS_FINISH


      SUBROUTINE DCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         REAL(q) :: DX(*),DY(*)

         !$omp dispatch
         CALL DCOPY(N,DX,INCX,DY,INCY)
      END SUBROUTINE DCOPY_OFFLOAD


      SUBROUTINE SCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         REAL(qs) :: DX(*),DY(*)

         !$omp dispatch
         CALL SCOPY(N,DX,INCX,DY,INCY)
      END SUBROUTINE SCOPY_OFFLOAD


      SUBROUTINE CCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         COMPLEX(qs) :: ZX(*),ZY(*)

         !$omp dispatch
         CALL CCOPY(N,ZX,INCX,ZY,INCY)
      END SUBROUTINE CCOPY_OFFLOAD


      SUBROUTINE ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         COMPLEX(q) :: ZX(*),ZY(*)

         !$omp dispatch
         CALL ZCOPY(N,ZX,INCX,ZY,INCY)
      END SUBROUTINE ZCOPY_OFFLOAD


      SUBROUTINE DAXPY_OFFLOAD(N,DA,DX,INCX,DY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         REAL(q) :: DX(*),DY(*)
         REAL(q) :: DA

         !$omp dispatch
         CALL DAXPY(N,DA,DX,INCX,DY,INCY)
      END SUBROUTINE DAXPY_OFFLOAD


      SUBROUTINE ZAXPY_OFFLOAD(N,ZA,ZX,INCX,ZY,INCY)
         IMPLICIT NONE
         INTEGER :: N,INCX,INCY
         COMPLEX(q) :: ZX(*),ZY(*)
         COMPLEX(q) :: ZA

         !$omp dispatch
         CALL ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
      END SUBROUTINE ZAXPY_OFFLOAD


      SUBROUTINE DSCAL_OFFLOAD(N,DA,DX,INCX)
         IMPLICIT NONE
         INTEGER :: N,INCX
         REAL(q) :: DX(*)
         REAL(q) :: DA

         !$omp dispatch
         CALL DSCAL(N,DA,DX,INCX)
      END SUBROUTINE DSCAL_OFFLOAD


      SUBROUTINE ZSCAL_OFFLOAD(N,ZA,ZX,INCX)
         IMPLICIT NONE
         INTEGER :: N,INCX
         COMPLEX(q) :: ZX(*)
         COMPLEX(q) :: ZA

         !$omp dispatch
         CALL ZSCAL(N,ZA,ZX,INCX)
      END SUBROUTINE ZSCAL_OFFLOAD


      SUBROUTINE ZDSCAL_OFFLOAD(N,DA,ZX,INCX)
         IMPLICIT NONE
         INTEGER :: N,INCX
         COMPLEX(q) :: ZX(*)
         REAL(q) :: DA

         !$omp dispatch
         CALL ZDSCAL(N,DA,ZX,INCX)
      END SUBROUTINE ZDSCAL_OFFLOAD


      REAL(q) FUNCTION DDOT_OFFLOAD(N,DX,INCX,DY,INCY)
         IMPLICIT NONE
         REAL(q) :: DX(*),DY(*)
         INTEGER :: N, INCX, INCY

         !$omp dispatch
         DDOT_OFFLOAD = DDOT(N,DX,INCX,DY,INCY)
      END FUNCTION DDOT_OFFLOAD


      COMPLEX(q) FUNCTION ZDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
         IMPLICIT NONE
         COMPLEX(q) :: ZX(*),ZY(*)
         INTEGER :: N, INCX, INCY

         !$omp dispatch
         ZDOTC_OFFLOAD = ZDOTC(N,ZX,INCX,ZY,INCY)
      END FUNCTION ZDOTC_OFFLOAD


      COMPLEX(qs) FUNCTION CDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
         IMPLICIT NONE
         COMPLEX(qs) :: ZX(*),ZY(*)
         INTEGER :: N, INCX, INCY

         !$omp dispatch
         CDOTC_OFFLOAD = CDOTC(N,ZX,INCX,ZY,INCY)
      END FUNCTION CDOTC_OFFLOAD


      SUBROUTINE DGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),X(*),Y(*)
         REAL(q) :: ALPHA,BETA
         INTEGER :: M,N,LDA,INCX,INCY,SIZE_Y
         CHARACTER(1) :: MODE
         !local variables
         INTEGER :: IDY

         IF (MODE=='N'.OR.MODE=='n') THEN
            IDY=M*INCY
         ELSE
            IDY=N*INCY
         ENDIF


         IF (ABS(BETA) > 0.0) THEN
            !$omp target enter data map(to:Y(1:IDY))
         ELSE
            !$omp target enter data map(alloc:Y(1:IDY))
         ENDIF

         !$omp dispatch
         CALL DGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

         !$omp target exit data map(from:Y(1:IDY))
      END SUBROUTINE DGEMV_OFFLOAD


      SUBROUTINE ZGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*),X(*),Y(*)
         COMPLEX(q) :: ALPHA,BETA
         INTEGER :: M,N,LDA,INCX,INCY,SIZE_Y
         CHARACTER(1) :: MODE
         !local variables
         INTEGER :: IDY

         IF (MODE=='N'.OR.MODE=='n') THEN
            IDY=M*INCY
         ELSE
            IDY=N*INCY
         ENDIF

         IF (ABS(BETA) > 0.0) THEN
            !$omp target enter data map(to:Y(1:IDY))
         ELSE
            !$omp target enter data map(alloc:Y(1:IDY))
         ENDIF

         !$omp dispatch
         CALL ZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

         !$omp target exit data map(from:Y(1:IDY))
      END SUBROUTINE ZGEMV_OFFLOAD


      SUBROUTINE DGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
         REAL(q) :: ALPHA,BETA
         INTEGER :: M,N,K,LDA,LDB,LDC,KA,KB
         CHARACTER(1) :: TRANSA,TRANSB

         IF (N < 1) RETURN

         !$omp dispatch
         CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      END SUBROUTINE DGEMM_OFFLOAD


      SUBROUTINE ZGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
         COMPLEX(q) :: ALPHA,BETA
         INTEGER :: M,N,K,LDA,LDB,LDC,KA,KB
         CHARACTER(1) :: TRANSA,TRANSB

         IF (N < 1) RETURN

         !$omp dispatch
         CALL ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      END SUBROUTINE ZGEMM_OFFLOAD


      SUBROUTINE SGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         IMPLICIT NONE
         REAL(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
         REAL(qs) :: ALPHA,BETA
         INTEGER :: M,N,K,LDA,LDB,LDC,KA,KB
         CHARACTER(1) :: TRANSA,TRANSB

         IF (N < 1) RETURN

         !$omp dispatch
         CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      END SUBROUTINE SGEMM_OFFLOAD


      SUBROUTINE CGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         IMPLICIT NONE
         COMPLEX(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
         COMPLEX(qs) :: ALPHA,BETA
         INTEGER :: M,N,K,LDA,LDB,LDC,KA,KB
         CHARACTER(1) :: TRANSA,TRANSB

         IF (N < 1) RETURN

         !$omp dispatch
         CALL CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      END SUBROUTINE CGEMM_OFFLOAD


      SUBROUTINE DTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),B(LDB,*)
         REAL(q) :: ALPHA
         INTEGER :: M,N,LDA,LDB,KA
         CHARACTER(1) :: SIDE,UPLO,TRANSA,DIAG

         !$omp dispatch
         CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      END SUBROUTINE DTRMM_OFFLOAD


      SUBROUTINE ZTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*),B(LDB,*)
         COMPLEX(q) :: ALPHA
         INTEGER :: M,N,LDA,LDB
         CHARACTER(1) :: SIDE,UPLO,TRANSA,DIAG

         !$omp dispatch
         CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      END SUBROUTINE ZTRMM_OFFLOAD


      SUBROUTINE DLACPY_OFFLOAD(UPLO,M,N,A,LDA,B,LDB)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),B(LDB,*)
         INTEGER :: LDA,LDB,M,N
         CHARACTER(1) :: UPLO

         IF (UPLO/='F'.AND.UPLO/='f') THEN
            CALL vtutor%bug("DLACPY_OFFLOAD: only UPLO=F (full matrix) is supported",__FILE__,__LINE__)
         ENDIF

#ifndef use_kernels_construct
         IF (M==LDA.AND.LDB==LDA) THEN
            !$omp dispatch
            CALL DCOPY(M*N,A,1,B,1)
            RETURN
         ENDIF
#endif
         !$omp target
         B(1:M,1:N)=A(1:M,1:N)
         !$omp end target
      END SUBROUTINE DLACPY_OFFLOAD

    END MODULE mmklblas


!***********************************************************************
!
!***********************************************************************
    include "mkl_lapack_omp_offload_lp64.f90"

    MODULE mmkllapack
      USE prec
      USE tutor
      USE moffload_struct_def, ONLY : OFFLOAD_ALLOW
      ! the information about which routines are fully GPU, hybrid or CPU only is taken from:
      ! https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-0/openmp-offload-for-onemkl.html
      USE onemkl_lapack_omp_offload_lp64, ONLY: DPOTRF, & ! fully GPU
                                                ZPOTRF, & ! fully GPU
                                                DGETRF, & ! hybrid,  (all computations on the CPU for n <= 256)
                                                ZGETRF, & ! hybrid,  (all computations on the CPU for n <= 256)
                                                DGETRS, & ! fully GPU
                                                ZGETRS, & ! fully GPU
                                                DTRTRI, & ! fully GPU
                                                ZTRTRI, & ! fully GPU
                                              ! DSYGV, & ! not present in onemkl_lapack_omp_offload_lp64 yet
                                                DSYGVD, & ! hybrid, partly on CPU (used instead of DSYGV)
                                                ZHEGVD, & ! ? not quite clear from documentation probably interface only
                                                DSYEV, &  ! interface only, all computations on CPU
                                                DSYEVD, & ! hybrid, partly on CPU
                                                DSYEVX, & ! hybrid, partly on CPU
                                                SSYEVX, & ! hybrid, partly on CPU
                                                ZHEEV, &  ! interface only, all computations on CPU
                                                ZHEEVD, & ! hybrid, partly on CPU
                                                ZHEEVX, & ! hybrid, partly on CPU
                                                CHEEVX    ! hybrid, partly on CPU

      PRIVATE

      PUBLIC :: MKLLAPACK_INIT, MKLLAPACK_FINISH
      PUBLIC :: DPOTRF_OFFLOAD, ZPOTRF_OFFLOAD
      PUBLIC :: DGETRF_OFFLOAD, ZGETRF_OFFLOAD
      PUBLIC :: DGETRS_OFFLOAD, ZGETRS_OFFLOAD
      PUBLIC :: DTRTRI_OFFLOAD, ZTRTRI_OFFLOAD
      PUBLIC :: DSYGV_OFFLOAD, ZHEGV_OFFLOAD
      PUBLIC :: DSYEV_OFFLOAD, ZHEEV_OFFLOAD
      PUBLIC :: SSYEVX_OFFLOAD, DSYEVX_OFFLOAD, ZHEEVX_OFFLOAD, CHEEVX_OFFLOAD

    CONTAINS

      SUBROUTINE MKLLAPACK_INIT()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      END SUBROUTINE MKLLAPACK_INIT


      SUBROUTINE MKLLAPACK_FINISH()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      END SUBROUTINE MKLLAPACK_FINISH


      SUBROUTINE DPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*)
         INTEGER :: INFO
         INTEGER :: N,LDA
         CHARACTER(1) :: UPLO
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL DPOTRF(UPLO,N,A,LDA,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE DPOTRF_OFFLOAD


      SUBROUTINE ZPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*)
         INTEGER :: INFO
         INTEGER :: N,LDA
         CHARACTER(1) :: UPLO
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL ZPOTRF(UPLO,N,A,LDA,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE ZPOTRF_OFFLOAD


      SUBROUTINE DGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*)
         INTEGER :: IPIV(*),INFO
         INTEGER :: M,N,LDA
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL DGETRF(M,N,A,LDA,IPIV,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE DGETRF_OFFLOAD


      SUBROUTINE ZGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*)
         INTEGER :: IPIV(*),INFO
         INTEGER :: M,N,LDA
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL ZGETRF(M,N,A,LDA,IPIV,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE ZGETRF_OFFLOAD


      SUBROUTINE DGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),B(LDB,*)
         INTEGER :: IPIV(*)
         INTEGER :: N,NRHS,LDA,LDB,INFO
         CHARACTER(1) :: TRANS
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE DGETRS_OFFLOAD


      SUBROUTINE ZGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*),B(LDB,*)
         INTEGER :: IPIV(*)
         INTEGER :: N,NRHS,LDA,LDB,INFO
         CHARACTER(1) :: TRANS
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE ZGETRS_OFFLOAD


      SUBROUTINE DTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*)
         INTEGER :: N,LDA,INFO
         CHARACTER(1) :: UPLO,DIAG
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL DTRTRI(UPLO,DIAG,N,A,LDA,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE DTRTRI_OFFLOAD


      SUBROUTINE ZTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*)
         INTEGER :: N,LDA,INFO
         CHARACTER(1) :: UPLO,DIAG
         !local variables
         INTEGER :: INFO_

         !$omp target enter data map(alloc: INFO_)
         !$omp dispatch
         CALL ZTRTRI(UPLO,DIAG,N,A,LDA,INFO_)

         !$omp target exit data map(from: INFO_)
         INFO = INFO_ * 100
      END SUBROUTINE ZTRTRI_OFFLOAD


      SUBROUTINE DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
         INTEGER :: INFO
         INTEGER :: ITYPE,N,LDA,LDB,LWORK
         CHARACTER(1) :: JOBZ,UPLO
         !local variables
         REAL(q), ALLOCATABLE :: WORK_(:)
         REAL(q) :: WORK_Q(1)
         INTEGER, ALLOCATABLE :: IWORK_(:)
         INTEGER :: IWORK_Q(1)
         INTEGER :: INFO_, LWORK_, LIWORK_

         IF ((JOBZ .NE. 'N') .AND. (JOBZ .NE. 'V')) THEN
            CALL vtutor%bug("DSYGV_OFFLOAD: invalid JOBZ",__FILE__,__LINE__)
         ENDIF

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: W(1:N), INFO_, WORK_Q, IWORK_Q)
         !$omp dispatch
         CALL DSYGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK_Q,-1,IWORK_Q,-1,INFO_)
         !$omp target exit data map(from: WORK_Q, IWORK_Q)
         LWORK_ = WORK_Q(1)
         LIWORK_ = IWORK_Q(1)
         ALLOCATE(WORK_(LWORK_), IWORK_(LIWORK_))

         !$omp target enter data map(alloc:WORK_, IWORK_)
         !$omp dispatch
         CALL DSYGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK_,LWORK_,IWORK_,LIWORK_,INFO_)
         !$omp target exit data map(from: W(1:N), INFO_)
         !$omp target exit data map(release: IWORK_, WORK_)
         DEALLOCATE(WORK_, IWORK_)
         INFO = INFO_ * 100
      END SUBROUTINE DSYGV_OFFLOAD

      ! SUBROUTINE DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
      !    IMPLICIT NONE
      !    REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
      !    INTEGER :: INFO
      !    INTEGER :: ITYPE,N,LDA,LDB,LWORK
      !    CHARACTER(1) :: JOBZ,UPLO
      !
      !    !$omp dispatch
      !    CALL DSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
      ! END SUBROUTINE DSYGV_OFFLOAD

      SUBROUTINE ZHEGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
        IMPLICIT NONE
        COMPLEX(q) :: A(LDA,*),B(LDB,*)
        REAL(q) :: W(*), RWORK(*)
        INTEGER :: INFO
        COMPLEX(q) :: WORK(*)
        INTEGER :: ITYPE,N,LDA,LDB,LWORK,LRWORK
        CHARACTER(1) :: JOBZ,UPLO
        !local variables
        COMPLEX(q), ALLOCATABLE :: WORK_(:)
        COMPLEX(q) :: WORK_Q(1)
        REAL(q), ALLOCATABLE :: RWORK_(:)
        REAL(q) :: RWORK_Q(1)
        INTEGER, ALLOCATABLE :: IWORK_(:)
        INTEGER :: IWORK_Q(1)
        INTEGER :: INFO_, LWORK_, LRWORK_, LIWORK_

         IF ((JOBZ .NE. 'N') .AND. (JOBZ .NE. 'V')) THEN
            CALL vtutor%bug("ZHEGV_OFFLOAD: invalid JOBZ",__FILE__,__LINE__)
         ENDIF

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc:IWORK_Q, RWORK_Q, WORK_Q, W(1:N), INFO_)
         !$omp dispatch
         CALL ZHEGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK_Q,-1,RWORK_Q,-1,IWORK_Q,-1,INFO_)
         !$omp target exit data map(from: IWORK_Q, RWORK_Q, WORK_Q)
         LIWORK_ = IWORK_Q(1)
         LRWORK_ = RWORK_Q(1)
         LWORK_ = WORK_Q(1)
         ALLOCATE(IWORK_(LIWORK_), RWORK_(LRWORK_), WORK_(LWORK_))

         !$omp target enter data map(alloc: IWORK_, RWORK_, WORK_)
         !$omp target enter data map(to: LIWORK_, LRWORK_, LWORK_)
         !$omp dispatch
         CALL ZHEGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK_,LWORK_,RWORK_,LRWORK_,IWORK_,LIWORK_,INFO_)
         !$omp target exit data map(from:W(1:N), INFO_)
         !$omp target exit data map(release:IWORK_, RWORK_, WORK_, LIWORK_, LRWORK_, LWORK_)
         DEALLOCATE(IWORK_, RWORK_, WORK_)
         INFO = INFO_ * 100
      END SUBROUTINE ZHEGV_OFFLOAD


      SUBROUTINE DSYEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),W(*),WORK(*)
         INTEGER :: INFO
         INTEGER :: N,LDA,LWORK
         CHARACTER(1) :: JOBZ,UPLO
         ! local variables
         REAL(q), ALLOCATABLE :: WORK_(:)
         REAL(q) :: WORK_Q(1)
         INTEGER :: LWORK_, INFO_
#ifdef divide_and_conquer
         INTEGER :: LIWORK_
         INTEGER, ALLOCATABLE :: IWORK_(:)
         INTEGER :: IWORK_Q(1)

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: W(1:N), INFO_, WORK_Q, IWORK_Q)
         !$omp dispatch
         CALL DSYEVD(JOBZ,UPLO,N,A,LDA,W,WORK_Q,-1,IWORK_Q,-1,INFO_)
         !$omp target exit data map(from: WORK_Q, IWORK_Q)
         LWORK_ = WORK_Q(1)
         LIWORK_ = IWORK_Q(1)
         ALLOCATE(WORK_(LWORK_), IWORK_(LIWORK_))

         !$omp target enter data map(alloc: WORK_, IWORK_)
         !$omp dispatch
         CALL DSYEVD(JOBZ,UPLO,N,A,LDA,W,WORK_,LWORK_,IWORK_,LIWORK_,INFO_)
         !$omp target exit data map(from: W(1:N), INFO_)
         !$omp target exit data map(release: WORK_, IWORK_)
         DEALLOCATE(WORK_, IWORK_)
         INFO = INFO_ * 100
#else         
         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: WORK_Q, W(1:N), INFO_)
         !$omp dispatch
         CALL DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK_Q,-1,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK_,LWORK_,INFO_)
         !$omp target exit data map(from: W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)
         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
#endif
      END SUBROUTINE DSYEV_OFFLOAD


      SUBROUTINE ZHEEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*)
         REAL(q) :: W(*),RWORK(*)
         INTEGER :: INFO
         COMPLEX(q) :: WORK(*)
         INTEGER :: N,LDA,LWORK
         CHARACTER(1) :: JOBZ,UPLO
         ! local variables
         COMPLEX(q), ALLOCATABLE :: WORK_(:)
         COMPLEX(q) :: WORK_Q(1)
         INTEGER :: INFO_, LWORK_
#ifdef divide_and_conquer
         INTEGER ::LRWORK_, LIWORK_
         INTEGER, ALLOCATABLE :: IWORK_(:)
         REAL(q), ALLOCATABLE :: RWORK_(:)
         REAL(q) :: RWORK_Q(1)
         INTEGER :: IWORK_Q(1)

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: W(1:N), INFO_, WORK_Q, RWORK_Q, IWORK_Q)
         !$omp dispatch
         CALL ZHEEVD(JOBZ,UPLO,N,A,LDA,W,WORK_Q,-1,RWORK_Q,-1,IWORK_Q,-1,INFO_)
         !$omp target exit data map(from: WORK_Q, RWORK_Q, IWORK_Q)
         LWORK_ = WORK_Q(1)
         LRWORK_ = RWORK_Q(1)
         LIWORK_ = IWORK_Q(1)
         ALLOCATE(WORK_(LWORK_), RWORK_(LRWORK_), IWORK_(LIWORK_))

         !$omp target enter data map(alloc: WORK_, RWORK_, IWORK_)
         !$omp dispatch
         CALL ZHEEVD(JOBZ,UPLO,N,A,LDA,W,WORK_,LWORK_,RWORK_,LRWORK_,IWORK_,LIWORK_,INFO_)
         !$omp target exit data map(from: W(1:N), INFO_)
         !$omp target exit data map(release: WORK_, RWORK_, IWORK_)
         DEALLOCATE(WORK_, RWORK_, IWORK_)
         INFO = INFO_ * 100
#else
         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: W(1:N), INFO_, WORK_Q)
         !$omp dispatch
         CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK_Q,-1,RWORK,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK_,LWORK_,RWORK,INFO_)
         !$omp target exit data map(from: W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)
         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
#endif
      END SUBROUTINE ZHEEV_OFFLOAD


      SUBROUTINE SSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
         IMPLICIT NONE
         REAL(qs) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
         INTEGER :: IFAIL(*),INFO
         REAL(qs) :: ABSTOL,VL,VU
         INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
         CHARACTER(1) :: JOBZ,RANGE,UPLO
         ! local variables
         REAL(qs), ALLOCATABLE :: WORK_(:)
         REAL(qs) :: WORK_Q(1)
         INTEGER :: LWORK_, INFO_

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: M, W(1:N), INFO_, WORK_Q)
         !$omp dispatch
         CALL SSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_Q,-1,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL SSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_,LWORK_,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: M, W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)

         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
      END SUBROUTINE SSYEVX_OFFLOAD


      SUBROUTINE DSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
         IMPLICIT NONE
         REAL(q) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
         INTEGER :: IFAIL(*),INFO
         REAL(q) :: VL,VU,ABSTOL
         INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
         CHARACTER(1) :: JOBZ,RANGE,UPLO
         ! local variables
         REAL(q), ALLOCATABLE :: WORK_(:)
         REAL(q) :: WORK_Q(1)
         INTEGER :: LWORK_, INFO_

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: M, W(1:N), INFO_, WORK_Q)
         !$omp dispatch
         CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_Q,-1,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_,LWORK_,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: M, W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)

         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
      END SUBROUTINE DSYEVX_OFFLOAD


      SUBROUTINE ZHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
         IMPLICIT NONE
         COMPLEX(q) :: A(LDA,*),Z(LDZ,*),WORK(*)
         REAL(q) :: W(*)
         INTEGER :: IFAIL(*),INFO
         REAL(q) :: VL,VU,ABSTOL,RWORK(*)
         INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
         CHARACTER(1) :: JOBZ,RANGE,UPLO
         ! local variables
         COMPLEX(q), ALLOCATABLE :: WORK_(:)
         COMPLEX(q) :: WORK_Q(1)
         INTEGER :: LWORK_, INFO_

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: M, W(1:N), INFO_, WORK_Q)
         !$omp dispatch
         CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_Q,-1,RWORK,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_,LWORK_,RWORK,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: M, W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)
         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
      END SUBROUTINE ZHEEVX_OFFLOAD


      SUBROUTINE CHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
         IMPLICIT NONE
         COMPLEX(qs) :: A(LDA,*),Z(LDZ,*),WORK(*)
         REAL(qs) :: W(*)
         INTEGER :: IFAIL(*),INFO
         REAL(qs) :: VL,VU,ABSTOL,RWORK(*)
         INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
         CHARACTER(1) :: JOBZ,RANGE,UPLO
         ! local variables
         COMPLEX(qs), ALLOCATABLE :: WORK_(:)
         COMPLEX(qs) :: WORK_Q(1)
         INTEGER :: LWORK_, INFO_

         ! Determine the optimal workspace size
         !$omp target enter data map(alloc: M, W(1:N), INFO_, WORK_Q)
         !$omp dispatch
         CALL CHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_Q,-1,RWORK,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: WORK_Q)
         LWORK_ = WORK_Q(1)
         ALLOCATE(WORK_(LWORK_))

         !$omp target enter data map(alloc: WORK_)
         !$omp dispatch
         CALL CHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK_,LWORK_,RWORK,IWORK,IFAIL,INFO_)
         !$omp target exit data map(from: M, W(1:N), INFO_)
         !$omp target exit data map(release: WORK_)
         DEALLOCATE(WORK_)
         INFO = INFO_ * 100
      END SUBROUTINE CHEEVX_OFFLOAD

    END MODULE mmkllapack
#endif // OMP_OFFLOAD
#endif // INTELMKL
