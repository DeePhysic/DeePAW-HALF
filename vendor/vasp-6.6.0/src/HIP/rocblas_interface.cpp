#include <rocblas.h>
#include <rocsolver.h>
#include <hip_runtime.h>
#include <rocm-core/rocm_version.h>
#include <iostream>
// USE_ROCSOLVER_XTRTRI activates the rocSOLVER xtrtri function instead of
// the rocBLAS xtrtri implementation. Performance needs testing.
// Toggle must match src/crayhip.F
// #define USE_ROCSOLVER_XTRTRI

#define HIP_ERR(stat)                                             \
    if(stat != hipSuccess)                                        \
    {                                                             \
        std::cerr << "HIP error: " << hipGetErrorString(stat) <<  \
        " in file " << __FILE__ << ":" << __LINE__ << std::endl;  \
        exit(-1);                                                 \
    }                                                             \

rocblas_operation findop(int mode) {
    switch (mode) {
        case 111 :
            return rocblas_operation_none;
        case 112 :
            return rocblas_operation_transpose;
        case 113 :
            return rocblas_operation_conjugate_transpose;
        default :
            printf("ERROR invalid rocblas_operation %d\n",mode);
            exit(1);
    }
    return rocblas_operation_none;
}
rocblas_fill findfill(int mode) {
    switch (mode) {
        case 121 :
            return rocblas_fill_upper;
        case 122 :
            return rocblas_fill_lower;
        case 123 :
            return rocblas_fill_full;
        printf("ERROR invalid rocblas_fill %d\n",mode);
            exit(1);
    }
    return rocblas_fill_upper;
}
rocblas_diagonal finddiag(int mode) {
    switch (mode) {
        case 131 :
            return rocblas_diagonal_non_unit;
        case 132 :
            return rocblas_diagonal_unit;
        printf("ERROR invalid rocblas_diagonal %d\n",mode);
            exit(1);
    }
    return rocblas_diagonal_non_unit;
}
rocblas_side findside(int mode) {
    switch (mode) {
        case 141 :
            return rocblas_side_left;
        case 142 :
            return rocblas_side_right;
        printf("ERROR invalid rocblas_side %d\n",mode);
            exit(1);
    }
    return rocblas_side_left;
}
rocblas_eform finditype(int mode) {
    switch (mode) {
        case 1 :
            return rocblas_eform_ax;
        case 2 :
            return rocblas_eform_abx;
        case 3 :
            return rocblas_eform_bax;
        default :
            printf("ERROR invalid rocblas_eform %d\n",mode);
            exit(1);
    }
    return rocblas_eform_ax;
}
rocblas_evect findjobz(int mode) {
    switch (mode) {
        case 211 :
            return rocblas_evect_original;
        case 213 :
            return rocblas_evect_none;
        default :
            printf("ERROR invalid rocblas_evect %d\n",mode);
            exit(1);
    }
    return rocblas_evect_none;
}
rocblas_erange finderange(int mode) {
    switch (mode) {
        case 231 :
            return rocblas_erange_all;
        case 232 :
            return rocblas_erange_value;
         case 233 :
            return rocblas_erange_index;
        default :
            printf("ERROR invalid rocblas_erange %d\n",mode);
            exit(1);
    }
}
extern "C" {
    void hip_init_rocblas(void *ptr) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_create_handle(handle);
    }
    void hip_finish_rocblas(void *ptr) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_destroy_handle(*handle);
    }
    void hip_rocblas_set_stream(void *ptr, void *stream) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        hipStream_t *hipstream = (hipStream_t *) stream;
        rocblas_set_stream(*handle,*hipstream);
    }
    void hip_rocblas_synchronize(void *ptr){
        rocblas_handle *handle = (rocblas_handle *) ptr;
        hipStream_t hipstream;
        rocblas_get_stream(*handle, &hipstream);
        HIP_ERR(hipStreamSynchronize(hipstream));
    }
    // ROCBLAS
    void hip_dcopy(void *ptr, int N, double *dx, int incx, double *dy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dcopy(*handle, N, dx, incx, dy, incy);
    }
    void hip_dcopy_strided_batched(void *ptr, int N, double *dx, int incx, int stridex, double *dy,
            int incy, int stridey, int batch_count) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dcopy_strided_batched(*handle, N, dx, incx, stridex, dy, incy, stridey, batch_count);
    }
    void hip_zcopy(void *ptr, int N, double _Complex *zx, int incx, double _Complex *zy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_double_complex *zy2 = reinterpret_cast<rocblas_double_complex*>(zy);
        rocblas_zcopy(*handle, N, zx2, incx, zy2, incy);
    }
    void hip_zcopy_strided_batched(void *ptr, int N, double _Complex *zx, int incx, int stridex,
            double _Complex *zy, int incy, int stridey, int batch_count) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_double_complex *zy2 = reinterpret_cast<rocblas_double_complex*>(zy);
        rocblas_zcopy_strided_batched(*handle, N, zx2, incx, stridex, zy2, incy, stridey, batch_count);
    }
    void hip_scopy(void *ptr, int N, float *dx, int incx, float *dy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_scopy(*handle, N, dx, incx, dy, incy);
    }
    void hip_ccopy(void *ptr, int N, float _Complex *zx, int incx, float _Complex *zy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *zx2 = reinterpret_cast<rocblas_float_complex*>(zx);
        rocblas_float_complex *zy2 = reinterpret_cast<rocblas_float_complex*>(zy);
        rocblas_ccopy(*handle, N, zx2, incx, zy2, incy);
    }
    void hip_daxpy(void *ptr, int N, double da, double *dx, int incx, double *dy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_daxpy(*handle, N, &da, dx, incx, dy, incy);
    }
    void hip_zaxpy(void *ptr, int N, double _Complex za, double _Complex *zx, int incx, double _Complex *zy, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_double_complex *zy2 = reinterpret_cast<rocblas_double_complex*>(zy);
        rocblas_double_complex *za2 = reinterpret_cast<rocblas_double_complex*>(&za);
        rocblas_zaxpy(*handle, N, za2, zx2, incx, zy2, incy);
    }
    void hip_dscal(void *ptr, int N, double da, double *dx, int incx) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dscal(*handle, N, &da, dx, incx);
    }
    void hip_zscal(void *ptr, int N, double _Complex za, double _Complex *zx, int incx) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_double_complex *za2 = reinterpret_cast<rocblas_double_complex*>(&za);
        rocblas_zscal(*handle, N, za2, zx2, incx);
    }
    void hip_zdscal(void *ptr, int N, double da, double _Complex *zx, int incx) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_zdscal(*handle, N, &da, zx2, incx);
    }
    void hip_ddot(void *ptr, int N, double *dx, int incx, double *dy, int incy, double *ddot) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_ddot(*handle, N, dx, incx, dy, incy, ddot);
    }
    void hip_zdotc(void *ptr, int N, double _Complex *zx, int incx, double _Complex *zy, int incy, double _Complex *zdotc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
        rocblas_double_complex *zy2 = reinterpret_cast<rocblas_double_complex*>(zy);
        rocblas_double_complex *zdotc2 = reinterpret_cast<rocblas_double_complex*>(zdotc);
        rocblas_zdotc(*handle, N, zx2, incx, zy2, incy, zdotc2);
    }
    void hip_cdotc(void *ptr, int N, float _Complex *zx, int incx, float _Complex *zy, int incy, float _Complex *cdotc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *zx2 = reinterpret_cast<rocblas_float_complex*>(zx);
        rocblas_float_complex *zy2 = reinterpret_cast<rocblas_float_complex*>(zy);
        rocblas_float_complex *cdotc2 = reinterpret_cast<rocblas_float_complex*>(cdotc);
        rocblas_cdotc(*handle, N, zx2, incx, zy2, incy, cdotc2);
    }
    void hip_dgemv(void *ptr, int mode, int m, int n, double alpha, double *A, int lda, double *x, int incx,
            double beta, double *y, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dgemv(*handle,findop(mode),m,n,&alpha,A,lda,x,incx,&beta,y,incy);
    }
    void hip_zgemv(void *ptr, int mode, int m, int n, double _Complex alpha, double _Complex *A, int lda,
            double _Complex *x, int incx, double _Complex beta, double _Complex *y, int incy) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *x2 = reinterpret_cast<rocblas_double_complex*>(x);
        rocblas_double_complex *y2 = reinterpret_cast<rocblas_double_complex*>(y);
        rocblas_double_complex *alpha2 = reinterpret_cast<rocblas_double_complex*>(&alpha);
        rocblas_double_complex *beta2 = reinterpret_cast<rocblas_double_complex*>(&beta);
        rocblas_zgemv(*handle,findop(mode),m,n,alpha2,A2,lda,x2,incx,beta2,y2,incy);
    }
    void hip_dgemm(void *ptr, int modeA, int modeB, int m, int n, int k, double alpha, double *A, int lda,
            double *B, int ldb, double beta, double *C, int ldc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dgemm(*handle,findop(modeA),findop(modeB),m,n,k,&alpha,A,lda,B,ldb,&beta,C,ldc);
    }
    void hip_zgemm(void *ptr, int modeA, int modeB, int m, int n, int k, double _Complex alpha, double _Complex *A, int lda,
            double _Complex *B, int ldb, double _Complex beta, double _Complex *C, int ldc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocblas_double_complex *C2 = reinterpret_cast<rocblas_double_complex*>(C);
        rocblas_double_complex *alpha2 = reinterpret_cast<rocblas_double_complex*>(&alpha);
        rocblas_double_complex *beta2 = reinterpret_cast<rocblas_double_complex*>(&beta);
        rocblas_zgemm(*handle,findop(modeA),findop(modeB),m,n,k,alpha2,A2,lda,B2,ldb,beta2,C2,ldc);
    }
    void hip_sgemm(void *ptr, int modeA, int modeB, int m, int n, int k, float alpha, float *A, int lda,
            float *B, int ldb, float beta, float *C, int ldc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_sgemm(*handle,findop(modeA),findop(modeB),m,n,k,&alpha,A,lda,B,ldb,&beta,C,ldc);
    }
    void hip_cgemm(void *ptr, int modeA, int modeB, int m, int n, int k, float _Complex alpha, float _Complex *A, int lda,
            float _Complex *B, int ldb, float _Complex beta, float _Complex *C, int ldc) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *A2 = reinterpret_cast<rocblas_float_complex*>(A);
        rocblas_float_complex *B2 = reinterpret_cast<rocblas_float_complex*>(B);
        rocblas_float_complex *C2 = reinterpret_cast<rocblas_float_complex*>(C);
        rocblas_float_complex *alpha2 = reinterpret_cast<rocblas_float_complex*>(&alpha);
        rocblas_float_complex *beta2 = reinterpret_cast<rocblas_float_complex*>(&beta);
        rocblas_cgemm(*handle,findop(modeA),findop(modeB),m,n,k,alpha2,A2,lda,B2,ldb,beta2,C2,ldc);
    }
    void hip_dtrmm(void *ptr, int side, int uplo, int transa, int diag, int m, int n, double alpha, double *A, int lda,
            double *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
#if ROCM_VERSION_MAJOR >= 6
        rocblas_dtrmm(*handle,findside(side),findfill(uplo),findop(transa),finddiag(diag),m,n,&alpha,A,lda,B,ldb,B,ldb);
#else
        rocblas_dtrmm(*handle,findside(side),findfill(uplo),findop(transa),finddiag(diag),m,n,&alpha,A,lda,B,ldb);
#endif
    }
    void hip_ztrmm(void *ptr, int side, int uplo, int transa, int diag, int m, int n, double _Complex alpha, double _Complex *A,
            int lda, double _Complex *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocblas_double_complex *alpha2 = reinterpret_cast<rocblas_double_complex*>(&alpha);
#if ROCM_VERSION_MAJOR >= 6
        rocblas_ztrmm(*handle,findside(side),findfill(uplo),findop(transa),finddiag(diag),m,n,alpha2,A2,lda,B2,ldb,B2,ldb);
#else
        rocblas_ztrmm(*handle,findside(side),findfill(uplo),findop(transa),finddiag(diag),m,n,alpha2,A2,lda,B2,ldb);
#endif
    }
    // ROCSOLVER
    void hip_dpotrf(void *ptr, int uplo, int n, double *A, int lda, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dpotrf(*handle,findfill(uplo),n,A,lda,info);
    }
    void hip_zpotrf(void *ptr, int uplo, int n, double _Complex *A, int lda, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocsolver_zpotrf(*handle,findfill(uplo),n,A2,lda,info);
    }
    void hip_dgetrf(void *ptr, int m, int n, double *A, int lda, int *ipiv, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dgetrf(*handle,m,n,A,lda,ipiv,info);
    }
    void hip_zgetrf(void *ptr, int m, int n, double _Complex *A, int lda, int *ipiv, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocsolver_zgetrf(*handle,m,n,A2,lda,ipiv,info);
    }
    void hip_dgetrs(void *ptr, int trans, int n, int nrhs, double *A, int lda, int *ipiv,
            double *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dgetrs(*handle,findop(trans),n,nrhs,A,lda,ipiv,B,ldb);
    }
    void hip_zgetrs(void *ptr, int trans, int n, int nrhs, double _Complex *A, int lda, int *ipiv,
            double _Complex *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocsolver_zgetrs(*handle,findop(trans),n,nrhs,A2,lda,ipiv,B2,ldb);
    }
#ifdef USE_ROCSOLVER_XTRTRI
    // The following ?TRTRI implementations use the newer rocsolver subroutines which matches the LAPACK API.
    void hip_dtrtri(void *ptr, int fill, int diag, int n, double *A, int lda, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dtrtri(*handle,findfill(fill),finddiag(diag),n,A,lda,info);
    }
    void hip_ztrtri(void *ptr, int fill, int diag, int n, double _Complex *A, int lda, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocsolver_ztrtri(*handle,findfill(fill),finddiag(diag),n,A2,lda,info);
    }
#else
    void hip_dtrtri(void *ptr, int fill, int diag, int n, double *A, int lda, double *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dtrtri(*handle,findfill(fill),finddiag(diag),n,A,lda,B,ldb);
    }
    void hip_ztrtri(void *ptr, int fill, int diag, int n, double _Complex *A, int lda,  double _Complex *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocblas_ztrtri(*handle,findfill(fill),finddiag(diag),n,A2,lda,B2,ldb);
    }
#endif
    void hip_dtrsm(void *ptr, int side, int fill, int trans, int diag, int m, int n, double alpha, double *A, int lda,
            double *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_dtrsm(*handle,findside(side),findfill(fill),findop(trans),finddiag(diag),m,n,&alpha,A,lda,B,ldb);
    }
    void hip_ztrsm(void *ptr, int side, int fill, int trans, int diag, int m, int n, double _Complex alpha,
            double _Complex *A, int lda, double _Complex *B, int ldb) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocblas_double_complex *alpha2 = reinterpret_cast<rocblas_double_complex*>(&alpha);
        rocblas_ztrsm(*handle,findside(side),findfill(fill),findop(trans),finddiag(diag),m,n,alpha2,A2,lda,B2,ldb);
    }
    void hip_dsygv(void *ptr, int itype, int jobz, int fill, int n, double *A, int lda, double *B, int ldb,
            double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dsygv(*handle,finditype(itype),findjobz(jobz),findfill(fill),n,A,lda,B,ldb,w,work,info);
    }
    void hip_dsygvd(void *ptr, int itype, int jobz, int fill, int n, double *A, int lda, double *B, int ldb,
            double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dsygvd(*handle,finditype(itype),findjobz(jobz),findfill(fill),n,A,lda,B,ldb,w,work,info);
    }
    // eigenvalue solver for general hermitian matrices using robust QR decomposition
    void hip_zhegv(void *ptr, int itype, int jobz, int fill, int n, double _Complex *A, int lda, double _Complex *B,
            int ldb, double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocsolver_zhegv(*handle,finditype(itype),findjobz(jobz),findfill(fill),n,A2,lda,B2,ldb,w,work,info);
    }
    // eigenvalue solver for general hermitian matrices using faster divide-and-conquer algorithm STEVD
    void hip_zhegvd(void *ptr, int itype, int jobz, int fill, int n, double _Complex *A, int lda, double _Complex *B,
            int ldb, double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
        rocsolver_zhegvd(*handle,finditype(itype),findjobz(jobz),findfill(fill),n,A2,lda,B2,ldb,w,work,info);
    }
    void hip_zhegvj(void *ptr, int itype, int jobz, int fill, int n, double _Complex *A, int lda, double _Complex *B,
            int ldb, double *res, int *n_sw, double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *B2 = reinterpret_cast<rocblas_double_complex*>(B);
	//MYCOMMENT: comment for now to make it compile with rocm5.0.2
        //rocsolver_zhegvj(*handle,finditype(itype),findjobz(jobz),findfill(fill),n,A2,lda,B2,ldb,0,res,10,n_sw,w,info);
    }
    void hip_dsyev(void *ptr, int jobz, int fill, int n, double *A, int lda, double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dsyev(*handle,findjobz(jobz),findfill(fill),n,A,lda,w,work,info);
    }
    // eigenvalue solver using faster divide-and-conquer algorithm STEVD
    void hip_dsyevd(void *ptr, int jobz, int fill, int n, double *A, int lda, double *w, double *work, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dsyevd(*handle,findjobz(jobz),findfill(fill),n,A,lda,w,work,info);
    }
    void hip_zheev(void *ptr, int jobz, int fill, int n, double _Complex *A, int lda, double *w, double *work,
            int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocsolver_zheev(*handle,findjobz(jobz),findfill(fill),n,A2,lda,w,work,info);
    }
    // eigenvalue solver using faster divide-and-conquer algorithm STEVD
    void hip_zheevd(void *ptr, int jobz, int fill, int n, double _Complex *A, int lda, double *w, double *work,
            int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocsolver_zheevd(*handle,findjobz(jobz),findfill(fill),n,A2,lda,w,work,info);
    }
    void hip_cheev(void *ptr, int jobz, int fill, int n, float _Complex *A, int lda, float *w, float *work,
            int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *A2 = reinterpret_cast<rocblas_float_complex*>(A);
        rocsolver_cheev(*handle,findjobz(jobz),findfill(fill),n,A2,lda,w,work,info);
    }
    // eigenvalue solver using faster divide-and-conquer algorithm STEVD
    void hip_cheevd(void *ptr, int jobz, int fill, int n, float _Complex *A, int lda, float *w, float *work,
            int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *A2 = reinterpret_cast<rocblas_float_complex*>(A);
        rocsolver_cheevd(*handle,findjobz(jobz),findfill(fill),n,A2,lda,w,work,info);
    }
    void hip_ssyevx(void *ptr, int jobz, int erange, int fill, int n, float *A, int lda, float vl, float vu,
            int il, int iu, float abstol, int *nev, float *W, float *Z, int ldz, int *ifail, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_ssyevx(*handle,findjobz(jobz),finderange(erange),findfill(fill),n,A,lda,vl,vu,il,iu,abstol,nev,W,Z,ldz,ifail,info);
    }
    void hip_dsyevx(void *ptr, int jobz, int erange, int fill, int n, double *A, int lda, double vl, double vu,
            int il, int iu, double abstol, int *nev, double *W, double *Z, int ldz, int *ifail, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocsolver_dsyevx(*handle,findjobz(jobz),finderange(erange),findfill(fill),n,A,lda,vl,vu,il,iu,abstol,nev,W,Z,ldz,ifail,info);
    }
    void hip_zheevx(void *ptr, int jobz, int erange, int fill, int n, double _Complex *A, int lda, double vl, double vu,
            int il, int iu, double abstol, int *nev, double *W, double _Complex *Z, int ldz, int *ifail, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_double_complex *A2 = reinterpret_cast<rocblas_double_complex*>(A);
        rocblas_double_complex *Z2 = reinterpret_cast<rocblas_double_complex*>(Z);
        rocsolver_zheevx(*handle,findjobz(jobz),finderange(erange),findfill(fill),n,A2,lda,vl,vu,il,iu,abstol,nev,W,Z2,ldz,ifail,info);
    }
    void hip_cheevx(void *ptr, int jobz, int erange, int fill, int n, float _Complex *A, int lda, double vl, double vu,
            int il, int iu, double abstol, int *nev, float *W, float _Complex *Z, int ldz, int *ifail, int *info) {
        rocblas_handle *handle = (rocblas_handle *) ptr;
        rocblas_float_complex *A2 = reinterpret_cast<rocblas_float_complex*>(A);
        rocblas_float_complex *Z2 = reinterpret_cast<rocblas_float_complex*>(Z);
        rocsolver_cheevx(*handle,findjobz(jobz),finderange(erange),findfill(fill),n,A2,lda,vl,vu,il,iu,abstol,nev,W,Z2,ldz,ifail,info);
    }
}
