#include "symbol.inc"
#ifdef OFFLOADING
!***********************************************************************
!
!***********************************************************************

    MODULE moffload_deepcopy
#ifdef ACC_OFFLOAD
      USE mopenacc_deepcopy
#endif // ACC_OFFLOAD
#ifdef OMP_OFFLOAD
      USE mopenmp_deepcopy
#endif // OMP_OFFLOAD
    END MODULE moffload_deepcopy


!***********************************************************************
!
!***********************************************************************

    MODULE moffload_fft
#ifdef NVCUDA
      USE mcufft
#endif // NVCUDA
#ifdef CRAYHIP
      USE mrocfft
#endif // CRAYHIP
#ifdef INTELMKL
      USE mmklfft
#endif // INTELMKL
    END MODULE moffload_fft


!***********************************************************************
!
!***********************************************************************

    MODULE moffload_blas
#ifdef NVCUDA
      USE mcublas, BLAS_OFFLOAD_INIT   => CUBLAS_INIT, &
                   BLAS_OFFLOAD_FINISH => CUBLAS_FINISH
#endif // NVCUDA
#ifdef CRAYHIP
      USE mrocblas, BLAS_OFFLOAD_INIT   => ROCBLAS_INIT, &
                    BLAS_OFFLOAD_FINISH => ROCBLAS_FINISH
#endif // CRAYHIP
#ifdef INTELMKL
      USE mmklblas, BLAS_OFFLOAD_INIT   => MKLBLAS_INIT, &
                    BLAS_OFFLOAD_FINISH => MKLBLAS_FINISH
#endif // INTELMKL
    END MODULE moffload_blas


!***********************************************************************
!
!***********************************************************************

    MODULE moffload_lapack
#ifdef NVCUDA
      USE mcusolver,  LAPACK_OFFLOAD_INIT   => CUSOLVER_INIT, &
                      LAPACK_OFFLOAD_FINISH => CUSOLVER_FINISH
#endif // NVCUDA
#ifdef CRAYHIP
      USE mrocsolver, LAPACK_OFFLOAD_INIT   => ROCSOLVER_INIT, &
                      LAPACK_OFFLOAD_FINISH => ROCSOLVER_FINISH
#endif // CRAYHIP
#ifdef INTELMKL
      USE mmkllapack, LAPACK_OFFLOAD_INIT   => MKLLAPACK_INIT, &
                      LAPACK_OFFLOAD_FINISH => MKLLAPACK_FINISH
#endif // INTELMKL
    END MODULE moffload_lapack


!***********************************************************************
!
!***********************************************************************

    MODULE moffload_scalapack
#ifdef NVCUDA
      USE mcusolvermp
      USE mcublasmp
#endif // NVCUDA
    END MODULE moffload_scalapack


!***********************************************************************
!
!***********************************************************************

    MODULE moffload_tutor
#ifdef ACC_OFFLOAD
      USE mopenacc_tutor, ONLY : VTUTOR_OFFLOAD_LOG   => ACC_VTUTOR_LOG, &
                                 VTUTOR_OFFLOAD_CHECK => ACC_VTUTOR_CHECK
#endif // ACC_OFFLOAD
    END MODULE moffload_tutor


!***********************************************************************
!
!***********************************************************************

    MODULE moffload

      USE moffload_struct_def
      USE moffload_deepcopy
      USE moffload_fft, ONLY : FFT_OFFLOAD_MAKEPLANS, &
                               FFT_OFFLOAD_MAKEPLAN, &
                               FFT_OFFLOAD_DESTROYPLAN
      USE moffload_blas, ONLY : BLAS_OFFLOAD_INIT, &
                                BLAS_OFFLOAD_FINISH
      USE moffload_lapack, ONLY : LAPACK_OFFLOAD_INIT, &
                                  LAPACK_OFFLOAD_FINISH
      USE moffload_tutor
      USE moffload_control

#ifdef ACC_OFFLOAD
      USE mopenacc, ONLY : INIT_OFFLOAD, &
                           ACC_SET_ASYNC_Q,ACC_SYNC_ASYNC_Q, &
                           GET_FREE_MEM_DEVICE,ACC_GET_BATCH_SIZE, &
                           SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX
#endif // ACC_OFFLOAF

#ifdef OMP_OFFLOAD
      USE mopenmp, ONLY : INIT_OFFLOAD, &
                          SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX
#ifndef IFX_BEFORE_20241008
      USE mopenmp, ONLY : OMPX_DUMP_MAPPING_TABLES
#endif
#endif // OMP_OFFLOAD

    END MODULE moffload

#endif // OFFLOADING
