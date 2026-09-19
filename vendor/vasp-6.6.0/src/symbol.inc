#define GBBSE
!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part --------------------------------
#if defined(USELIBXC5) && ! defined(USELIBXC)
#define USELIBXC
#endif
#ifdef single_BLAS
!
!   use single precission BLAS and LAPACK calls (usually on Cray)
!
#define ZGEMM  CGEMM
#define ZGEMV  CGEMV
#define ZDSCAL CSSCAL
#define ZDOTC  CDOTC
#define ZAXPY  CAXPY
#define ZHEEV  CHEEV
#define ZHEGV  CHEGV
#define ZTRMM  CTRMM
#define ZTRTRI CTRTRI
#define ZPOTRF CPOTRF
#define PZPOTRF PCPOTRF
#define PZTRTRI PCTRTRI
#define PZHEEVX PCHEEVX
#define ZCOPY   CCOPY
#define ZGETRF  CGETRF
#define ZGETRS  CGETRS
#define ZGEFA   CGEFA
#define ZGECO   CGECO
#define ZGEDI   CGEDI
#define ZHEEVX  CHEEVX

#define DGEMM  SGEMM
#define DGEMV  SGEMV
#define DGEGV  SGEGV
#define DDOT   SDOT
#define DSCAL  SSCAL
#define DAXPY  SAXPY
#define DSYEV  SSYEV
#define DSYGV  SSYGV
#define DTRMM  STRMM
#define DTRTRI STRTRI
#define DPOTRF SPOTRF
#define DGETRF SGETRF
#define DGETRS SGETRS
#define PDPOTRF PSPOTRF
#define PDTRTRI PSTRTRI
#define PDSYEVX PSSYEVX
#define DCOPY   SCOPY
#define DSYEVX  SSYEVX

#endif
#ifdef NGZhalf
!
!   charge density: half grid mode Z direction
!
#define realmode
#define  zm(n)             (n/2+1)
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P     N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N3 /= 1) FACTM=2
#define  SETSQRT          IF (N3 /= 1) FACTM=SQRT(2._q)
#define  SETISQRT         IF (N3 /= 1) FACTM=1./SQRT(2._q)
#define  MULFACT          FACTM*
#elif defined(NGXhalf)
!
!   charge density: half grid mode X direction
!
#define realmode
#define  xm(n)             (n/2+1)
#define  xmin(n)            0
#define  zm(n)              n
#define  SETN1fromN1P     N1=N1P
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N1 /= 1) FACTM=2
#define  SETSQRT          IF (N1 /= 1) FACTM=SQRT(2._q)
#define  SETISQRT         IF (N1 /= 1) FACTM=1./SQRT(2._q)
#define  MULFACT          FACTM*
#else
!
!   charge density: full grid mode
!
#define  zm(n)              n
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P      N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1          FACTM=1
#define  SETFACT
#define  SETSQRT
#define  SETISQRT
#define  MULFACT
#endif
#ifdef realmode
!
!   charge density real
!
#define RGRID             REAL(q)
#define OVERLAP           REAL(q)
#define DIMREAL(n)        n*2
#define CHGCONJG(n)       (n)
#else
!
!   charge density complex
!
#define RGRID            COMPLEX(q)
#define OVERLAP          COMPLEX(q)
#define DIMREAL(n)       n
#define CHGCONJG(n)      CONJG(n)
#endif

#ifdef wNGZhalf
!
!   wavefunctions: half grid mode for Z direction
!
#define gammareal
#elif defined(wNGXhalf)
!
!   wavefunctions: half grid mode for X direction
!
#define gammareal
#else
!
!   wavefunctions: full grid mode
!
#endif
#ifdef  gammareal
!
!   wavefunctions real (gamma only)
!
#define GREAL(n)  REAL(n,KIND=q)
#define GCMPLX(n,m,q)  (n)
#define GDEF      REAL(q)
#define GDEFS     REAL(qs)
#define GDEFN     REAL(qn)
#define M_sum_g   M_sum_d
#define M_sumf_g  M_sumf_d
#define M_sum_g8  M_sum_d8
#define M_bcast_g M_bcast_d
#define M_bcast_g_from M_bcast_d_from
#define M_allgatherv_g M_allgatherv_d
#define CPROJ     GPROJ
#define GCONJG(n) (n)
#define GAIMAG(n) 0._q
#define GGEMM     DGEMM
#define GLANGE    DLANGE
#define GGEMV     DGEMV
#define GGETRF    DGETRF
#define GGETRS    DGETRS
#define GDOTC     DDOT
#define PGLANXX   PDLANSY
#define PGGETRF   PDGETRF
#define PGGETRI   PDGETRI
#define GDSCAL    DSCAL
#define GAXPY     DAXPY
#define GDOT      DDOT
#define PGGEMM    PDGEMM
#define PGGEMR2D  PDGEMR2D 
#define PGPOTRF   PDPOTRF
#define PGGEADD   PDGEADD
#define PGLACPY   PDLACPY
#define PGPOTRF   PDPOTRF
#define PGTRSM    PDTRSM
#define zero      0._q
#define one       1._q
#define m_        2*
#define ndata     1
#define ndata8    1_qi8
#define trans     'T'
#define gammaarg   NPLWKC,NINDPC,WVSCAL,WVSCAI
#define gammaargNK NPLWKC(NK),NINDPC(1,NK),WVSCAL(1,NK),WVSCAI(1,NK)
#else
!
!   wavefunctions complex
!
#define GREAL(n)  (n)
#define GCMPLX(n,m,q)  CMPLX(n,m,KIND=q)
#define GDEF      COMPLEX(q)
#define GDEFS     COMPLEX(qs)
#define GDEFN     COMPLEX(qn)
#define M_sumf_g  M_sumf_z
#define M_sum_g   M_sum_z
#define M_sum_g8  M_sum_z8
#define M_bcast_g M_bcast_z
#define M_bcast_g_from M_bcast_z_from
#define M_allgatherv_g M_allgatherv_z
#define GCONJG(n) CONJG(n)
#define GAIMAG(n) AIMAG(n)
#define GGEMM     ZGEMM
#define GLANGE    ZLANGE
#define GGEMV     ZGEMV
#define GGETRF    ZGETRF
#define GGETRS    ZGETRS
#define GDOTC     ZDOTC
#define PGLANXX   PZLANHE
#define PGGETRF   PZGETRF
#define PGGETRI   PZGETRI
#define GDSCAL    ZDSCAL
#define GAXPY     ZAXPY
#define GDOT      ZDOTC
#define PGGEMM    PZGEMM
#define PGGEMR2D  PZGEMR2D 
#define PGPOTRF   PZPOTRF
#define PGGEADD   PZGEADD
#define PGLACPY   PZLACPY
#define PGPOTRF   PZPOTRF
#define PGTRSM    PZTRSM
#define zero      (0._q,0._q)
#define one       (1._q,0._q)
#define trans     'C'
#define m_
#define ndata     2
#define ndata8    2_qi8
#define gammaarg
#define gammaargNK
#endif
!
!   common definitions
!
#define ztrans     'C'
#define dtrans     'T'
#define zzero      (0._q,0._q)
#define zone       (1._q,0._q)
#define dzero      0._q
#define done       1._q

#if defined(MPI)
!
!   mpi parallel macros
!
#define CALLMPI(x) CALL x
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_stop('VASP aborting ...'); stop

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#elif defined(MPI_CHAIN)
!
!   mpi nudged chain version
!
#define CALLMPI(x)
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_stop('VASP aborting ...'); stop

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#else
#define CALLMPI(x)
#define CALLMPI_C(x)
#define nodes
#define do_io
#define io_begin
#define io_end
#undef use_shmem
#undef shmem_rproj
#undef shmem_bcast_buffer
#endif
#ifdef F90_T3D
!
!   T3D macros
!
!#define IVDEP
#endif
#if defined(T3D_SMA)
!
!   T3D SHMEM communication
!
#endif
#ifdef CRAY_MPP
!
!   CRAY_MPP macros
!
#define MPI_double_precision MPI_real
#define MPI_double_complex   MPI_complex
#endif
#ifdef USE_ERF
!
!  use instrinsic ERF
!
#define ERRF(x) ERF(x)
#define ERRFC(x) ERFC(x)
#endif

#define dwrite DWRITE
#ifdef debug
!
!  debugging primitives
!
#define DWRITE WRITE(0,*)
#define DWRITE0 do_io WRITE(0,*)
#else
#define DWRITE  !
#define DWRITE0 !
#endif

#define CW CPTWFP

#ifdef nonlr_single
!
! single precision non local projection operators
!
#define RPROJ   RPROJS
#define qn      qs
#define DGEMVn  SGEMV
#define DGEMMn  SGEMM
#else
#define qn      q
#define DGEMVn  DGEMV
#define DGEMMn  DGEMM
#endif


#if  ! (defined(NGXhalf)  || defined(NGYhalf) ||  defined(NGZhalf))
#define RHOLM_complex
#endif

#ifdef dotiming
!
!  some timing primitives
!
#define usetmg REAL(q) :: TV,TC,TV0,TC0
#define statmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV0,TC0)
#define stotmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') 'timing ',TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#define stotmgt(X) CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') X,TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#define addtmg(T1,T2) CALLMPI(MPI_barrier(WDES%COMM%MPI_COMM,ierror)); CALL VTIME(TV,TC); T1=T1+TV-TV0;T2=T2+TC-TC0; TV0=TV;TC0=TC
#define addtmga(T1,T2) CALL VTIME(TV,TC); T1=T1+TV-TV0;T2=T2+TC-TC0; TV0=TV;TC0=TC
#define wrttmg(X,T1,T2) IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') X,T1,T2
#else
#define usetmg
#define statmg
#define stotmg
#define stotmgt(X)
#define addtmg(T1,T2)
#define addtmga(T1,T2)
#define wrttmg(X,T1,T2)
#endif
!
! OpenMP macros
!
#ifdef _OPENMP45
#define _OPENMP40
#define _OPENMP
#endif
#ifdef _OPENMP40
#define _OPENMP
#endif
#ifdef _OPENMP
#define __omp_id               omp_id
#define __omp_get_max_threads  ,OMP_GET_MAX_THREADS()
#define __omp_get_num_threads  ,OMP_GET_NUM_THREADS()
#define __omp_get_thread_num   ,OMP_GET_THREAD_NUM()+1
#define __omp_arg(x)           ,x
#define __omp_arg_id           __omp_arg(__omp_id)
#define __omp_add(x)           x
#else
#define __omp_id               1
#define __omp_get_max_threads
#define __omp_get_num_threads
#define __omp_get_thread_num
#define __omp_arg(x)
#define __omp_arg_id
#define __omp_add(x)
#endif
!
! profiling macros
!
#ifdef PROFILING
#define PROFILING_START(x)     CALL START_PROFILING(x)
#define PROFILING_STOP(x)      CALL STOP_PROFILING(x)
#else
#define PROFILING_START(x)
#define PROFILING_STOP(x)
#endif
!
! shmem macros
!
#ifdef shmem_rproj
#define use_shmem
#endif
#ifdef shmem_bcast_buffer
#define use_shmem
#endif
#ifdef shmem_elphon
#define use_shmem
#endif
!
! quadruple precision
!
#ifdef qd_emulate
#define QDPREAL                TYPE(qd_real)
#define QDPCOMPLEX             TYPE(qd_complex)
#define QDPD(x)                qdreal(REAL(x,KIND=q))
#define TOPRINT(x)             TOREAL(x)
#define QDPC(x,y )             to_qdc_qd2( x,y )
#define QDPC2QDPREAL( x )      to_qd_qdc( x )
#ifdef __GNUC__
#define QDPP(x)                qdreal('x')
#define QDFP(x)                qdreal('x')
#else
#define QDPP(x)                qdreal(#x)
#define QDFP(x)                qdreal(#x)
#endif
#else
#define TOPRINT(x)             x
#define QDPC2QDPREAL( x )      REAL( x, qd )
#define QDPREAL                REAL(qd)
#define QDPCOMPLEX             COMPLEX(qd)
#define QDPD(x)                REAL(x,KIND=qd)
#define QDPP(x)                REAL(x,KIND=qd)
#define QDPC(x,y)              CMPLX(x,y,KIND=qd)
#define IDENTITYMACRO(x)       x
#define QDFP(x)                IDENTITYMACRO(x)_qd
#endif
!
! for the VASP2WANNIER90 interface
!
#ifdef VASP2WANNIER90v2
#define VASP2WANNIER90
#endif
!
! for the DFTD4 interface
!
#ifdef DFTD4_API_V3
#define DFTD4
#endif
!
!
! CUDA includes
!
#ifdef CUDA_GPU
#include "cuda_fft.inc"
#include "cuda_macros.inc"
#include "cuda_profiling.inc"
#endif
!
! SIMD related definitions
!
#include "simd.inc"
!
! memalign macros
!
#if   defined(memalign64)
#define memalign
#define __align_BTS            64
#define __align_DBL            8
#define __align_DIM(x)         ((x+__align_DBL-1)/__align_DBL)*__align_DBL
#elif defined(memalign32)
#define memalign
#define __align_BTS            32
#define __align_DBL            4
#define __align_DIM(x)         ((x+__align_DBL-1)/__align_DBL)*__align_DBL
#else
#define __align_DIM(x)         x
#endif

!
! Macros for PGI/NV HPC compilers version specific code
!
#ifdef __PGIC__
#define NV_MAJOR_VER __PGIC__
#define NV_MINOR_VER __PGIC_MINOR__
#else
#ifdef __PGIF90__
#define NV_MAJOR_VER __PGIF90__
#define NV_MINOR_VER __PGIF90_MINOR__
#elif __NVCOMPILER
#define NV_MAJOR_VER __NVCOMPILER_MAJOR__
#define NV_MINOR_VER __NVCOMPILER_MINOR__
#endif
#endif

#ifdef NV_MAJOR_VER
#define PGI_BEFORE_XX_X 1
#define PGI_BEFORE_24_3 (NV_MAJOR_VER < 24 || (NV_MAJOR_VER == 24 && NV_MINOR_VER < 3))
#define PGI_BEFORE_24_5 (NV_MAJOR_VER < 24 || (NV_MAJOR_VER == 24 && NV_MINOR_VER < 5))
#define PGI_BEFORE_24_7 (NV_MAJOR_VER < 24 || (NV_MAJOR_VER == 24 && NV_MINOR_VER < 7))
#define PGI_BEFORE_24_11 (NV_MAJOR_VER < 24 || (NV_MAJOR_VER == 24 && NV_MINOR_VER < 11))
#define PGI_BEFORE_25_1 (NV_MAJOR_VER < 25 || (NV_MAJOR_VER == 25 && NV_MINOR_VER < 1))
#define PGI_BEFORE_25_9 (NV_MAJOR_VER < 25 || (NV_MAJOR_VER == 25 && NV_MINOR_VER < 9))
#else
#define PGI_BEFORE_XX_X 0
#define PGI_BEFORE_24_3 0
#define PGI_BEFORE_24_5 0
#define PGI_BEFORE_24_7 0
#define PGI_BEFORE_24_11 0
#define PGI_BEFORE_25_1 0
#define PGI_BEFORE_25_9 0
#endif

#if PGI_BEFORE_24_5
#undef CUBLASMP
#endif

#if PGI_BEFORE_24_3
#undef CUSOLVERMP
#endif

!
! OpenACC macros
!
#ifdef ACC_OFFLOAD
#define SET_ACC_ASYNC_Q(x)     CALL ACC_SET_ASYNC_Q(x)
#define __IF_ASYNC__           IF(OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)

#ifdef gammareal
#define WGCOPY                 WDCOPY
#define WGGEMV                 WDGEMV
#define WGGEMM                 WDGEMM
#define WGDOTC                 WDDOT
#define WGDSCAL                WDSCAL
#else
#define WGCOPY                 WZCOPY
#define WGGEMV                 WZGEMV
#define WGGEMM                 WZGEMM
#define WGDOTC                 WZDOTC
#define WGDSCAL                WZDSCAL
#endif

#define __TUTOR__(x, y, u, v)  CALL VTUTOR_OFFLOAD_LOG(x, y)

#else
#define SET_ACC_ASYNC_Q(x)
#define __IF_ASYNC__
#define WPZPOTRF               PZPOTRF
#define WPZHEEV                PZHEEV
#define WPSSYEVD               PSSYEVD
#define WPDSYEVD               PDSYEVD
#define WPDSYEVX               PDSYEVX
#define WPZHEEVD               PZHEEVD
#define WPZHEEVX               PZHEEVX
#define WPCHEEVX               PCHEEVX
#define WPSSYEVX               PSSYEVX
#define WPZGEMR2D              PZGEMR2D
#define WPDGEMR2D              PDGEMR2D
#define WPZTRMR2D              PZTRMR2D
#define WPDTRMR2D              PDTRMR2D
#define WPSGEMR2D              PSGEMR2D
#define WPCGEMR2D              PCGEMR2D
#define WPSTRMR2D              PSTRMR2D
#define WPCTRMR2D              PCTRMR2D
#define WPZTRMR2D              PZTRMR2D
#define WPDTRMR2D              PDTRMR2D
#define WPZGEMM                PZGEMM
#define WPDGEMM                PDGEMM
#define WPCGEMM                PDGEMM
#define WPSGEMM                PSGEMM
#define WPSTRADD               PSTRADD
#define WPCTRADD               PCTRADD
#define WPZTRADD               PZTRADD
#define WPDTRADD               PDTRADD
#define WPCGEADD               PCGEADD
#define WPSGEADD               PSGEADD
#define WPZGEADD               PZGEADD
#define WPDGEADD               PDGEADD
#ifdef gammareal
#define WGCOPY                 DCOPY
#define WGGEMV                 DGEMV
#define WGGEMM                 DGEMM
#define WGDOTC                 DDOT
#define WGDSCAL                DSCAL
#else
#define WGCOPY                 ZCOPY
#define WGGEMV                 ZGEMV
#define WGGEMM                 ZGEMM
#define WGDOTC                 ZDOTC
#define WGDSCAL                ZDSCAL
#endif
#define __TUTOR__(x, y, u, v)  CALL vtutor%write(x, y, filename = u, linenumber = v)
#endif

!
! routines replaced in LAPACK >=3.6
!
#ifdef LAPACK36
#define DGEGV DGGEV
#endif

!
! Macros for HDF5 error check
!
#ifdef VASP_HDF5
#define VH5_CHECK(e) call vh5_error(e,__FILE__,__LINE__)
#else
#define VH5_CHECK(e)
#endif

!
! Macros for GNU version specific code
!
#ifdef __GNUC__
#define GNU_BEFORE_7_5 (__GNUC__ < 7 || (__GNUC__ == 7 && __GNUC_MINOR__ < 5))
#endif

!
! For machine learning
!
#if defined(MPI) || defined(MPI_CHAIN)
#define ML_AVAILABLE
#endif

!
! Extra safe initializations (+ overflow protections)
!
#ifdef extra_safe_init
#define DOESI
#else
#define DOESI !
#endif

!
! Macros for memory estimation
!
#define VMEM_ARR_SIZE(x)   VMEM_DIM_SIZE(x,SIZE(x,KIND=qi8))
#define VMEM_DIM_SIZE(x,y) VMEM_SIZE(x)*y
#define VMEM_SIZE(x)       0.125_q*STORAGE_SIZE(x)

!
! Offloading related macros
!
#if defined(ACC_OFFLOAD) || defined(OMP_OFFLOAD)
#define OFFLOADING
#endif

!
! combined OpenMP and OpenACC macros
!
#if defined(_OPENMP) && ! defined(OFFLOADING)
#define __NOOFF_omp_arg(x)     ,x
#else
#define __NOOFF_omp_arg(x)
#endif

#ifdef OFFLOADING
! Line is included when   OFFLOADING
#define __GPU
! Line is a comment when  OFFLOADING
#define __CPU          !!
#define PUSH_OFFLOAD_ON(x)     CALL OFFLOAD_ON_PUSH(x)
#define POP_OFFLOAD_ON         CALL OFFLOAD_ON_POP
#define __IF_OFFLOAD_OFF__(x)  (x.AND.(.NOT.OFFLOAD_ON))
#define ONCPU(x)               CALL OFFLOAD_ON_PUSH(.FALSE.,.FALSE.) ; x ; CALL OFFLOAD_ON_POP(.FALSE.)
!make sure wrappers are correctly set for
!gamma and non-gamma versions:
#ifdef gammareal
#define WGCOPY                 WDCOPY
#define WGGEMV                 WDGEMV
#define WGGEMM                 WDGEMM
#define WGDOTC                 WDDOT
#define WGDSCAL                WDSCAL
#else
#define WGCOPY                 WZCOPY
#define WGGEMV                 WZGEMV
#define WGGEMM                 WZGEMM
#define WGDOTC                 WZDOTC
#define WGDSCAL                WZDSCAL
#endif // gammareal
#else
! Line is included when  !OFFLOADING
#define __CPU
! Line is a comment when !OFFLOADING
#define __GPU          !!
#define PUSH_OFFLOAD_ON(x)
#define POP_OFFLOAD_ON
#define __IF_OFFLOAD_OFF__(x)  x
#define ONCPU(x)               x
!replace blas, lapack, and fft wrapper calls with
!normal calls if no offloading is used:
#define WDCOPY                 DCOPY
#define WZCOPY                 ZCOPY
#define WCCOPY                 CCOPY
#define WSCOPY                 SCOPY
#define WDAXPY                 DAXPY
#define WZAXPY                 ZAXPY
#define WDSCAL                 DSCAL
#define WZSCAL                 ZSCAL
#define WZDSCAL                ZDSCAL
#define WDDOT                  DDOT
#define WZDOTC                 ZDOTC
#define WCDOTC                 CDOTC
#define WDGEMV                 DGEMV
#define WZGEMV                 ZGEMV
#define WDGEMM                 DGEMM
#define WSGEMM                 SGEMM
#define WZGEMM                 ZGEMM
#define WCGEMM                 CGEMM
#define WDTRMM                 DTRMM
#define WZTRMM                 ZTRMM
#define WDLACPY                DLACPY
#define WDPOTRF                DPOTRF
#define WZPOTRF                ZPOTRF
#define WDGETRF                DGETRF
#define WZGETRF                ZGETRF
#define WDGETRS                DGETRS
#define WZGETRS                ZGETRS
#define WDTRTRI                DTRTRI
#define WZTRTRI                ZTRTRI
#define WDSYGV                 DSYGV
#define WDSYEV                 DSYEV
#define WDSYEVX                DSYEVX
#define WSSYEVX                SSYEVX
#define WZHEGV                 ZHEGV
#define WZHEEV                 ZHEEV
#define WZHEEVX                ZHEEVX
#define WCHEEVX                CHEEVX
#ifdef gammareal
#define WGCOPY                 DCOPY
#define WGGEMV                 DGEMV
#define WGGEMM                 DGEMM
#define WGDOTC                 DDOT
#define WGDSCAL                DSCAL
#else
#define WGCOPY                 ZCOPY
#define WGGEMV                 ZGEMV
#define WGGEMM                 ZGEMM
#define WGDOTC                 ZDOTC
#define WGDSCAL                ZDSCAL
#endif
#define WFFTBAS                FFTBAS_PLAN
#define WFFTBRC                FFTBRC_PLAN
#define WFFTBAS_MU             FFTBAS_PLAN_MU
#define WFFTBRC_MU             FFTBRC_PLAN_MU
#define WFFTBAS_MPI            FFTBAS_PLAN_MPI
#define WFFTBRC_MPI            FFTBRC_PLAN_MPI
#endif

#if defined(_OPENACC) && !defined(ACC_OFFLOAD)
! Line is included when   _OPENACC and !ACC_OFFLOAD
#define CPACC
#else
! Line is a comment when !_OPENACC or   ACC_OFFLOAD
#define CPACC          !!
#endif

#ifdef ACC_OFFLOAD
! As of this point _OPENACC will be set when ACC_OFFLOAD
#define _OPENACC
! Line is included when   ACC_OFFLOAD
#define GPACC
#else
! Line is a comment when !ACC_OFFLOAD
#define GPACC          !!
#endif

#ifdef _OPENACC
! Line is included when   _OPENACC
#define DOACC
! Line is a comment when  _OPENACC
#define NOACC          !!
#else
! Line is included when  !_OPENACC
#define NOACC
! Line is a comment when !_OPENACC
#define DOACC          !!
#endif

#if defined(_OPENMP) && !defined(OMP_OFFLOAD)
! Line is included when   _OPENMP and !OMP_OFFLOAD
#define CPOMP
#else
! Line is a comment when !_OPENMP or   OMP_OFFLOAD
#define CPOMP          !!
#endif

#ifdef OMP_OFFLOAD
! As of this point _OPENMP will be set when OMP_OFFLOAD
#define _OPENMP
! Line is included when   OMP_OFFLOAD
#define GPOMP
! Line is a comment when  OMP_OFFLOAD
#define NOOFF          !!
#else
! Line is a comment when !OMP_OFFLOAD
#define GPOMP          !!
! Line is included when   OMP_OFFLOAD
#define NOOFF
#endif

#ifdef _OPENMP
! Line is included when   _OPENMP
#define DOOMP
! Line is a comment when  _OPENMP
#define NOOMP          !!
#else
! Line is included when  !_OPENMP
#define NOOMP
! Line is a comment when !_OPENMP
#define DOOMP          !!
#endif

!
! Macros for Intel compilers version specific code
!
#ifdef __INTEL_COMPILER_BUILD_DATE
#define IFX_BEFORE_YYYYMMDD 1
#if __INTEL_COMPILER_BUILD_DATE < 20241008
#define IFX_BEFORE_20241008 1
#endif
#if __INTEL_COMPILER_BUILD_DATE < 20250605
#define IFX_BEFORE_20250605 1
#endif
#if __INTEL_COMPILER_BUILD_DATE < 20251001
#define IFX_BEFORE_20251001 1
#endif
#else
#define IFX_BEFORE_YYYYMMDD 0
#define IFX_BEFORE_20241008 0
#define IFX_BEFORE_20250605 0
#define IFX_BEFORE_20251001 0
#endif

#if IFX_BEFORE_YYYYMMDD && defined(OFFLOADING)
!
! Workaround1: prevent nesting of OMP TARGET regions
#define IFXWAR1
!              lines preceded by IFXW1 will be inactive
#define IFXW1          !!
!
! Workaround2: update dope vector on GPU
#define IFXWAR2
!              lines preceded by IFXW2 will be active
#define IFXW2
!
! Workaround3: remove present checks
#define IFXWAR3
!              lines preceded by IFXW3 will be inactive
#define IFXW3          !!
! Workaround4: PRIVATE inside of OMP PARALLEL DO does not work use FIRSTPRIVATE possibly only a problem in ifx 2025! but also it is a problem in ifx <2025.2 to use loop bounds that are variables defined in the ! active region 
#define IFXWAR4
!              lines preceded by IFXW4 will be inactive
#define IFXW4          !!
! Workaround5: Problems with ERF(X) in active region!
#define IFXWAR5
!              lines preceded by IFXW5 will be inactive
#define IFXW5          !!
!              lines preceded by IFXW5 will be inactive
! Workaround6: limit threads in nested parallel regions with reduction clauses, as intel spawns too many threads otherwise
#define IFXWAR6
!
! Workaround7: Array assignments  A = B in active regions only work if the data is contingous in memory.
#define IFXWAR7
!              lines preceded by IFXW7 will be inactive
#define IFXW7          !!
!
! Workaround9: optional array dummy arguments may not be declared with explicit sizes.
#define IFXWAR9
#else
#define IFXW1
#define IFXW2          !!
#define IFXW3
#define IFXW4
#define IFXW5
#define IFXW7
#endif

#if IFX_BEFORE_YYYYMMDD && defined(OFFLOADING)
! activate workdistribute constructs for ifx version 2025.3.0 and newer, used in array assignments
#define IFXW8
#else
! deactivate workdistrbute constructs if version is not known
#define IFXW8          !!
#endif
#if IFX_BEFORE_20251001 && defined(OFFLOADING)
! deactivate workdistrbute constructs for ifx version before 2025.3
#define IFXW8          !!
#endif

!
! Macros for cray compilers version specific code
!

#ifdef __DCRAYFTN_VERSION
#ifdef OFFLOADING
! to correctly attach derived type pointers (check ACC_ATTACH_ASYNC calls in openacc.F for comparison)
#define FTNW3
#define CRAYWAR3
! to circumvent a possible compilerbug when using MPICH communication.
#define CRAYWAR4
! to circumvent a compiler error for OMP single
#define FTNW5
#define CRAYWAR5

#if __DCRAYFTN_VERSION <= 20
! Workaround1: prevent inlining of functions
#define CRAYWAR1
! Workaround2: GG_base.F needs to be compiled without omp support. check GG_base.F for example
#define CRAYWAR2
#endif

#if __DCRAYFTN_VERSION >= 20
! Workaround1: prevent inlining of functions
#define CRAYWAR4
#endif

#endif
#else
! lines preceded by FTNW3 will be inactives
#define FTNW3          !!
! lines preceded by FTNW5 will be inactives
#define FTNW5          !!
#endif


!
! General OpenMP offloading workarounds for all compilers
!
#ifdef OFFLOADING
#ifdef OMP_OFFLOAD
! the firstprivate default mapping of scalar values fails for optional parameters of subroutines
! if they are not present both in cray and intel. This is currently part of the spec, maybe fixed in OMP6.1
#define OFFWAR_OPT_PARAM
! Fortran array syntax copies fail on Intel (see IFXW7) and sometimes AMD if the data in the array
! are non-contingous.
#define OFFWAR_ARRAY_SYNTAX
#endif
#endif
