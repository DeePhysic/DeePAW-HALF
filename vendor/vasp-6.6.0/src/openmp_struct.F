#include "symbol.inc"
    MODULE mopenmp_struct_def
      ! overall threading strategy

      !> @ref openmp :
      !> number of threads to be used in \"general\":
      !> by default this is set to the maximum number of available threads,
      !> mopenmp::omp_nthreads=omp_get_max_threads()
      INTEGER :: omp_nthreads = 1

      ! specific routines

      !> @ref openmp :
      !> number of threads used in ::m_alltoall_d_omp
      !> (experimental and unused for now).
      INTEGER :: omp_nthreads_alltoall = 1

      !> @ref openmp :
      !> number of threads used in the domain decomposition of the
      !> real space projection operators.
      !> Default: mopenmp::omp_nthreads_nonlr_rspace=mopenmp::omp_nthreads
      INTEGER :: omp_nthreads_nonlr_rspace = 1

      !> @ref openmp :
      !> number of threads to be used in connection with OpenACC:
      !> in most instances this involves a distribution over orbitals
      INTEGER :: omp_nthreads_acc = 1

      !> @ref openmp :
      !> if set to true the domain of the real space projection operators
      !> is decomposed over (x,y)-planes, in a round-robin fashion, if set
      !> to false the projectors are distributed over (z,y)-columns.
      LOGICAL :: omp_nonlr_planewise = .TRUE.

      !> @ref openmp :
      !> if set to true ::fftmakeplan or ::fftmakeplan_mpi will call
      !> dfftw_init_threads and set mopenmp::omp_dfftw_init_threads=.FALSE.
      LOGICAL :: omp_dfftw_init_threads = .TRUE.

#ifdef OMP_OFFLOAD
      !> asynchronous execution stream (default ACC_ASYNC_Q = 0: synchronous execution)
      INTEGER :: ACC_ASYNC_Q = 0

      !> The device number to which this MPI-rank connects itself
      !> by means of OMP_SET_DEFAULT_DEVICE in INIT_OFFLOAD.
      INTEGER :: OMP_MY_DEVICE_NUM

      !> wrapper around OMP_TARGET_IS_PRESENT
      LOGICAL, EXTERNAL :: IS_ON_DEVICE

      INTERFACE
         FUNCTION OMP_TARGET_IS_PRESENT(PTR, DEVICE_NUM) BIND(C)
            USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr, c_int
            TYPE(c_ptr), VALUE :: PTR
            INTEGER(c_int), VALUE :: DEVICE_NUM
            INTEGER(c_int) :: OMP_TARGET_IS_PRESENT
         END FUNCTION OMP_TARGET_IS_PRESENT
      END INTERFACE
#endif
    END MODULE mopenmp_struct_def

#ifdef OMP_OFFLOAD
!****************** FUNCTION IS_ON_DEVICE ******************************
!
!***********************************************************************

      LOGICAL FUNCTION IS_ON_DEVICE(VAR)
      USE mopenmp_struct_def, ONLY : OMP_MY_DEVICE_NUM, OMP_TARGET_IS_PRESENT
      USE iso_c_binding, ONLY : c_loc
      INTEGER, TARGET :: VAR(*)
      IS_ON_DEVICE = (OMP_TARGET_IS_PRESENT(c_loc(VAR), OMP_MY_DEVICE_NUM) == 1)
      END FUNCTION IS_ON_DEVICE
#endif
