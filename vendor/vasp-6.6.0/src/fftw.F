#include "symbol.inc"

!************************* SUBROUTINE FFTMAKEPLAN **********************
!
!> @details @ref openmp :
!> if mopenmp::omp_dfftw_init_threads=.TRUE. this subroutine calls
!> dfftw_init_threads, and sets mopenmp::omp_dfftw_init_threads=.FALSE. to
!> make sure dfftw_init_threads is called only once.
!
!***********************************************************************

      SUBROUTINE FFTMAKEPLAN(c,grid)

      USE prec
      USE mgrid_struct_def
      USE mopenmp_struct_def, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) grid
      COMPLEX(q) c(*), cdummy

      ! local variables
      INTEGER nx,ny,nz
      TYPE(c_ptr) plan
      INTEGER(c_int) EFFORT

!$    INTEGER IERROR, THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

      nx=grid%ngx
      ny=grid%ngy
      nz=grid%ngz

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (.NOT.grid%real2cplx) THEN
         CALL dfftw_plan_dft_3d(plan,nx,ny,nz, &
              c, c, &
              FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(plan)
         CALL dfftw_plan_dft_3d(plan,nx,ny,nz,&
              c, c, &
              FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(plan)
      ELSE
         CALL dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,&
              c,c, &
              EFFORT)
         CALL dfftw_destroy_plan(plan)
         CALL dfftw_plan_dft_c2r_3d(plan,nx,ny,nz,&
              c,c, &
              EFFORT)
         CALL dfftw_destroy_plan(plan)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
      END SUBROUTINE FFTMAKEPLAN


!************************* SUBROUTINE FFTBAS_PLAN **********************
!
!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN(C,GRID,ISIGN)
      USE prec
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID

      COMPLEX(q) :: C(*)
      INTEGER    :: ISIGN

      ! local variables
      TYPE(c_ptr) PLAN
      INTEGER(c_int) EFFORT

#if defined(fftw_cache_plans) || defined(__NEC_TUNE__)
      TYPE(c_ptr), SAVE :: PLAN_FWD = c_null_ptr, PLAN_BWD = c_null_ptr
      INTEGER, SAVE :: NFFT_FWD, NFFT_BWD
#endif

!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbas_plan')

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
#if !defined(fftw_cache_plans) && !defined(__NEC_TUNE__)
      IF (ISIGN.LE.0) THEN
         CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_FORWARD, EFFORT)
      ELSE
         CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_BACKWARD,EFFORT)
      ENDIF
#else // fftw_cache_plans or __NEC_TUNE__
      IF (ISIGN.LE.0) THEN
!!         write(*,*) 'forward:',plan_fwd,nfft_fwd,grid%nplwv,grid%ngx,grid%ngy,grid%ngz
         IF (c_associated(PLAN_FWD) .AND. NFFT_FWD == GRID%NPLWV) THEN
!!            write(*,*) 'using plan_fwd:',plan_fwd
            PLAN=PLAN_FWD
         ELSE
            IF (c_associated(PLAN_FWD)) CALL dfftw_destroy_plan(PLAN_FWD)
            CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_FORWARD,EFFORT)
            PLAN_FWD=PLAN
            NFFT_FWD=GRID%NPLWV
!!            write(*,*) 'generated plan_fwd:',plan_fwd,nfft_fwd
         ENDIF
      ELSE
!!         write(*,*) 'bckward:',plan_bwd,nfft_bwd,grid%nplwv,grid%ngx,grid%ngy,grid%ngz
         IF (c_associated(PLAN_BWD) .AND. NFFT_BWD == GRID%NPLWV) THEN
!!            write(*,*) 'using plan_bck:',plan_bwd
            PLAN=PLAN_BWD
         ELSE
            IF (c_associated(PLAN_BWD)) CALL dfftw_destroy_plan(PLAN_BWD)
            CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_BACKWARD,EFFORT)
            PLAN_BWD=PLAN
            NFFT_BWD=GRID%NPLWV
!!            write(*,*) 'generated plan_bck:',plan_bwd,nfft_bwd
         ENDIF
      ENDIF
#endif // fftw_cache_plans or __NEC_TUNE__
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      CALL dfftw_execute_dft(PLAN,C,C)

#if !defined(fftw_cache_plans) && !defined(__NEC_TUNE__)
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
#endif // fftw_cache_plans or __NEC_TUNE__

      PROFILING_STOP('fftbas_plan')

      RETURN
      END SUBROUTINE FFTBAS_PLAN


!************************* SUBROUTINE FFTBRC_PLAN **********************
!
!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN(C,GRID,ISIGN)
      USE prec
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      COMPLEX(q) C(*)
      INTEGER ISIGN

      ! local variables
      TYPE(c_ptr) PLAN
      INTEGER(c_int) EFFORT

!$    INTEGER THREADS, IERROR
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbrc_plan')

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (ISIGN.LE.0) THEN
         CALL dfftw_plan_dft_r2c_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,EFFORT)
      ELSE
         CALL dfftw_plan_dft_c2r_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,EFFORT)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      IF (ISIGN.LE.0) THEN
         CALL dfftw_execute_dft_r2c(PLAN,C,C)
      ELSE
         CALL dfftw_execute_dft_c2r(PLAN,C,C)
      ENDIF
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbrc_plan')

      RETURN
      END SUBROUTINE FFTBRC_PLAN


!************************* SUBROUTINE FFT1D_C2C_PLAN *******************
!
!> One dimensional FFT
!
!***********************************************************************

    SUBROUTINE FFT1D_C2C_PLAN(C,NFFT,ISIGN)
      USE prec
      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      COMPLEX(q) :: C(*)
      INTEGER    :: NFFT
      INTEGER    :: ISIGN

      ! local variables
      TYPE(c_ptr) PLAN

      IF (ISIGN<=0) THEN
         CALL dfftw_plan_dft_1d(PLAN,NFFT,C,C,FFTW_FORWARD, FFTW_ESTIMATE)
      ELSE
         CALL dfftw_plan_dft_1d(PLAN,NFFT,C,C,FFTW_BACKWARD,FFTW_ESTIMATE)
      ENDIF

      CALL dfftw_execute(PLAN,C,C)

      CALL dfftw_destroy_plan(PLAN)

    END SUBROUTINE FFT1D_C2C_PLAN


!****************** SUBROUTINE FFTMAKEPLAN_MU **************************
!
!***********************************************************************

    SUBROUTINE FFTMAKEPLAN_MU(N,C,LDC,GRID)
      USE prec
      USE mgrid_struct_def
      USE mopenmp_struct_def, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC

      ! local variables
      TYPE(c_ptr) :: PLAN
!!      INTEGER, POINTER :: INULLPTR => NULL()

      INTEGER(c_int) :: EFFORT

!$    INTEGER :: THREADS, ierror
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftmakeplan_mu')

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      IF (.NOT.GRID%REAL2CPLX) THEN
         CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(PLAN)

         CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(PLAN)
      ELSE
         CALL dfftw_plan_many_dft_r2c(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      EFFORT)
         CALL dfftw_destroy_plan(PLAN)

         CALL dfftw_plan_many_dft_c2r(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      EFFORT)
         CALL dfftw_destroy_plan(PLAN)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftmakeplan_mu')

      RETURN
    END SUBROUTINE FFTMAKEPLAN_MU


!****************** SUBROUTINE FFTBAS_PLAN_MU **************************
!
!***********************************************************************

    SUBROUTINE FFTBAS_PLAN_MU(N,C,LDC,GRID,ISN)
      USE prec
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(c_ptr) :: PLAN

      INTEGER(c_int) :: DIRECTION,EFFORT

!$    INTEGER :: THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbas_plan_mu')

!=======================================================================
! Call FFTW
!=======================================================================
      IF (ISN.LE.0) THEN
         DIRECTION=FFTW_FORWARD
      ELSE
         DIRECTION=FFTW_BACKWARD
      ENDIF

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif

!
! How many OpenMP threads are available?
!
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!
! Create batched 3d fft plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                               C(1), GRID%NGPTAR, 1, LDC, &
                               C(1), GRID%NGPTAR, 1, LDC, &
                               DIRECTION, EFFORT)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

!
! Execute plan
!
      CALL dfftw_execute_dft(PLAN,C(1),C(1))

!
! Destroy plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbas_plan_mu')

      RETURN
    END SUBROUTINE FFTBAS_PLAN_MU


!****************** SUBROUTINE FFTBRC_PLAN_MU **************************
!
!***********************************************************************

    SUBROUTINE FFTBRC_PLAN_MU(N,C,LDC,GRID,ISN)
      USE prec
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(c_ptr) :: PLAN

      INTEGER(c_int) :: EFFORT

!$    INTEGER :: THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbrc_plan_mu')

!=======================================================================
! Call FFTW
!=======================================================================
      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif

!
! How many OpenMP threads are available?
!
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!
! Create and execute batched 3d r2c or c2r fft plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      IF (ISN.LE.0) THEN
         CALL dfftw_plan_many_dft_r2c(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      EFFORT)
         CALL dfftw_execute_dft_r2c(PLAN,C(1),C(1))
      ELSE
         CALL dfftw_plan_many_dft_c2r(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      EFFORT)
         CALL dfftw_execute_dft_c2r(PLAN,C(1),C(1))
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

!
! Destroy plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbrc_plan_mu')

      RETURN
    END SUBROUTINE FFTBRC_PLAN_MU

#ifdef MPI
!************************* SUBROUTINE FFTMAKEPLAN_MPI ******************
!
!  FFTW requires a plan therefore this calling interface is provided
!  which calls the FFTBAS and FFTBRC routine for generating these
!  plans
!
!> @details @ref openmp :
!> if mopenmp::omp_dfftw_init_threads=.TRUE. this subroutine calls
!> dfftw_init_threads, and sets mopenmp::omp_dfftw_init_threads=.FALSE.
!> to make sure dfftw_init_threads is called only once.
!
!***********************************************************************

    SUBROUTINE FFTMAKEPLAN_MPI(A,GRID)
      USE prec
      USE mgrid_struct_def
      USE mopenmp_struct_def, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)

      ! local variables
      TYPE(c_ptr) :: planx, plany, planz
      INTEGER :: NX,NY,NZ,IDX,IDY,IDZ

      INTEGER(c_int) :: EFFORT

!$    INTEGER IERROR,THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IDX=NX
      IDY=NY
      IDZ=NZ

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (GRID%RC%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(planx)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_c2r(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             A(1), NZ, 1, IDZ+2 , &
                             EFFORT)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_BACKWARD, EFFORT)
         ENDIF
         CALL dfftw_destroy_plan(planz)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_r2c(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ+2 , &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             EFFORT)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_FORWARD, EFFORT)
         ENDIF
         CALL dfftw_destroy_plan(planz)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RC%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(planx)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
    END SUBROUTINE


!****************** SUBROUTINE FFTBAS_PLAN_MPI *************************
!
!   3-d parallel complex to complex fast fourier transformation
!   written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
!***********************************************************************

    SUBROUTINE FFTBAS_PLAN_MPI(A,GRID,ISIGN)
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN          !  direction of fft
     ! local variables
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
      TYPE(c_ptr) :: planx, plany, planz
      INTEGER :: NX,NY,NZ,IDX,IDY,IDZ

      INTEGER(c_int) EFFORT  !  planning effort

!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbas_plan_mpi')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))

      IDX=NX
      IDY=NY
      IDZ=NZ

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
!=======================================================================
! create plans
!=======================================================================
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (ISIGN==1) THEN
         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%RL_FFT%NCOL > 0) THEN
!            WRITE(*,*) 'complex to complex'
            CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                              A(1), NZ, 1, IDZ , &
                              A(1), NZ, 1, IDZ , &
                              FFTW_BACKWARD, EFFORT)
         ENDIF
      ELSE
         IF (GRID%RL_FFT%NCOL > 0) THEN
!            WRITE(*,*) 'detected inverse complex to complex'
            CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                              A(1), NZ, 1, IDZ , &
                              A(1), NZ, 1, IDZ , &
                              FFTW_FORWARD, EFFORT)
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_FORWARD, EFFORT)
         ENDIF

         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_FORWARD, EFFORT)
         ENDIF
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
      IF (ISIGN==1) THEN
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            CALL dfftw_execute_dft(planz,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
      ELSE
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            CALL dfftw_execute_dft(planz,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL  dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
      ENDIF

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      IF (GRID%RC%NCOL > 0) THEN
         call dfftw_destroy_plan(planx)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         call dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         call dfftw_destroy_plan(planz)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbas_plan_mpi')

      RETURN
    END SUBROUTINE FFTBAS_PLAN_MPI


!****************** SUBROUTINE FFTBRC_PLAN_MPI *************************
!
!   3-d parallel real to complex fast fourier transformation
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
!
!***********************************************************************

    SUBROUTINE FFTBRC_PLAN_MPI(A,GRID,ISIGN)
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid_struct_def

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN          !  direction of fft
     ! local variables
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
      TYPE(c_ptr) :: planx, plany, planz
      INTEGER :: NX,NY,NZ,IDX,IDY,IDZ

      INTEGER(c_int) EFFORT  !  planning effort

!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftbrc_plan_mpi')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,MAX(GRID%MPLWV,1))
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,MAX(GRID%MPLWV,1))

      IDX=NX
      IDY=NY
      IDZ=NZ

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
!=======================================================================
! create plans
!=======================================================================
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (ISIGN==1) THEN
         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%RL_FFT%NCOL > 0) THEN
!            WRITE(*,*) 'detected real to complex'
            CALL dfftw_plan_many_dft_c2r(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                              A(1), NZ, 1, (IDZ+2)/2 , &
                              A(1), NZ, 1, IDZ+2 , &
                              EFFORT)
         ENDIF
      ELSE
         IF (GRID%RL_FFT%NCOL > 0) THEN
!            WRITE(*,*) 'detected inverse real to complex'
            CALL dfftw_plan_many_dft_r2c(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                              A(1), NZ, 1, IDZ+2 , &
                              A(1), NZ, 1, (IDZ+2)/2 , &
                              EFFORT)
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_FORWARD, EFFORT)
         ENDIF

         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_FORWARD, EFFORT)
         ENDIF
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
      IF (ISIGN==1) THEN
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            CALL dfftw_execute_dft_c2r(planz,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
      ELSE
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            CALL dfftw_execute_dft_r2c(planz,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL  dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
      ENDIF

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      IF (GRID%RC%NCOL > 0) THEN
         call dfftw_destroy_plan(planx)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         call dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         call dfftw_destroy_plan(planz)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbrc_plan_mpi')

      RETURN
    END SUBROUTINE FFTBRC_PLAN_MPI
#endif // MPI
