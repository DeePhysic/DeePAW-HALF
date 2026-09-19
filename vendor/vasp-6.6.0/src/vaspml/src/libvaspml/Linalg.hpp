#ifndef LINALG_HPP
#define LINALG_HPP

#include "ParallelEnvironment.hpp" // IWYU pragma: keep
#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm> // min
#include <type_traits>

/*================================================================================================+
 | Compile without CBLAS and LAPACKE C wrappers for BLAS and LAPACK.
 +================================================================================================*/
#ifdef VASPML_NO_CBLAS
#warning "Compiling VASPml without CBLAS/LAPACKE, this is not recommended, expect poor performance!"

/*================================================================================================+
 | Standard build, assume existence of CBLAS and LAPACKE, or alternatively, cuBLAS or similar.
 |
 | The actual BLAS and LAPACK libraries which are linked to the application are determined at the
 | linking stage. However, for some choices there are some peculiarities regarding the correct
 | headers to include. Currently, three main options here:
 |
 | * If Intel oneAPI compiler is used we make use of MKL automatically, except if explicitly asked
 |   not to use it with VASPML_NO_MKL.
 | * Explicitly ask for using Intel's MKL via setting VASPML_USE_MKL in the makefile, requires
 |   special header names.
 | * Use standard headers for CBLAS and LAPACKE, some specialties:
 |   * NVIDIA brings their own OpenBLAS but in "lp64" subdirectory.
 |   * NEC does not seem to offer LAPACKE.
 +================================================================================================*/
#else // VASPML_NO_CBLAS

// Sanity-check MKL-specific flags.
#if defined(VASPML_USE_MKL) && defined(VASPML_NO_MKL)
#error "Conflicting preprocessor macros encountered: VASPML_USE_MKL + VASPML_NO_MKL!"
#endif

// Automatically use MKL if an Intel compiler is used, except if explicitly disabled.
#ifndef VASPML_NO_MKL
#if defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_COMPILER)
#define VASPML_USE_MKL
#endif
#endif

/*================================================================================================+
 | If MKL is requested we need to use their non-standard header names.
 +================================================================================================*/
#if defined(VASPML_USE_MKL)

#include <mkl_cblas.h>
#include <mkl_lapacke.h>

/*================================================================================================+
 | Standard CBLAS and LAPACKE header names.
 +================================================================================================*/
#else

/*================================================================================================+
 | NVIDIA-specific headers.
 +================================================================================================*/
#if defined(__NVCOMPILER)

#include <lp64/cblas.h>
#include <lp64/lapacke.h>

/*================================================================================================+
 | NEC-specific headers.
 +================================================================================================*/
#elif defined(__NEC__)

#include <cblas.h>
// NEC does not provide a C LAPACK interface, we have to use Fortran functions instead.

/*================================================================================================+
 | Standard headers.
 +================================================================================================*/
#else // Compiler-specific headers.

#include <cblas.h>
#include <lapacke.h>

#endif // Compiler-specific headers.

#endif // Intel MKL.

#endif // VASPML_NO_CBLAS

/*================================================================================================+
 | Headers for GPU offloading of CBLAS/LAPACKE calls.
 |
 | In addition to regular CBLAS/LAPACKE there could be libraries for offloading to GPU. Whether the
 | code is compiled with GPU support has already been determined in ParallelEnvironment.hpp. Hence,
 | we reuse macros from there to determine the GPU status.
 +================================================================================================*/

/*================================================================================================+
 | NVIDIA-specific headers.
 +================================================================================================*/
#ifdef VASPML_PALGO_GPU_NVIDIA

#include <cublas_v2.h>

#endif // Variants of VASPML_PALGO_GPU

namespace vaspml
{
template<class T>
class ShmemArray2D;
template<class T>
class ShmemArray;
}

namespace vaspml::linalg
{

/**************************************************************************************************
 * VASPml internal representation of transpose state.
 **************************************************************************************************/
enum class Transpose : Int
{
    NOTRANS = 111,
    TRANS = 112
};

/*================================================================================================+
 | Compile without CBLAS and LAPACKE C wrappers for BLAS and LAPACK.
 +================================================================================================*/
#ifdef VASPML_NO_CBLAS

using IntHost = Int;

/*================================================================================================+
 | Standard build, further specify integer types and transposing mechanism.
 +================================================================================================*/
#else // VASPML_NO_CBLAS

// Intel MKL provides its own integer type.
#if defined(VASPML_USE_MKL)
using IntHost = MKL_INT;
using TransposeHost = enum CBLAS_TRANSPOSE;
// NEC uses its own integer type.
#elif defined(__NEC__)
using IntHost = cblas_int_t;
using TransposeHost = CBLAS_TRANSPOSE;
// OpenBLAS commes with its own integer type.
#elif defined(OPENBLAS_VERSION)
using IntHost = blasint;
using TransposeHost = enum CBLAS_TRANSPOSE;
// Others may use the reference implementation type.
#else
using IntHost = int;
using TransposeHost = CBLAS_TRANSPOSE;
#endif // Library-specific using statements.

// Now with the transpose type defined we can make shortcuts.
inline TransposeHost noTransHost = CblasNoTrans;
inline TransposeHost transHost = CblasTrans;
constexpr IntHost    one = 1;

/// Function to convert VASPml transpose state to host library transpose specifier.
TransposeHost toHost(Transpose transpose);

#endif // VASPML_NO_CBLAS

/*================================================================================================+
 | If GPU libraries are used, add integer types and transposing mechanism separately.
 +================================================================================================*/
#ifdef VASPML_PALGO_GPU_NVIDIA

using IntDevice = Int;
using TransposeDevice = cublasOperation_t;

inline TransposeDevice noTransDevice = CUBLAS_OP_N;
inline TransposeDevice transDevice = CUBLAS_OP_T;
constexpr IntDevice    oneDevice = 1;

/// Function to convert VASPml transpose state to device library transpose specifier.
TransposeDevice toDevice(Transpose transpose);

#endif // Variants of VASPML_PALGO_GPU

template<typename T>
const Real* get_data_pointer_const(const T& data)
{
    const Real* ptr;
    if constexpr (std::is_same<Vec1Real, T>::value)
    {
        const Real& first_element = data[0];
        ptr = &first_element;
    }
    else if constexpr (std::is_same<Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<const Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<ShmemArray<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else if constexpr (std::is_same<ShmemArray2D<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else
    {
        //const String type = typeid(T).name();
        //global::tutor.bug( "const Real* ERROR:get_data_pointer_const( const T& data )
        //Supplied type not supported\n"
        //                                            "Supplied type " + type );
        return 0;
    }

    return ptr;
}

template<typename T>
Real* get_data_pointer(T& data)
{
    Real* ptr;
    if constexpr (std::is_same<Vec1Real, T>::value)
    {
        Real& first_element = data[0];
        ptr = &first_element;
    }
    else if constexpr (std::is_same<Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<ShmemArray<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else if constexpr (std::is_same<ShmemArray2D<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else
    {
        //const String type = typeid(T).name();
        //global::tutor.bug( "Real* get_data_pointer( T& data ) Supplied type not supported\n"
        //                         "Supplied type " + type );
        return 0;
    }

    return ptr;
}
/*******************************************************************************************
 * compute L2 norm of a Real vector
 *
 * @param vector input vector of which norm is computed
 * @param size   number of elements of type T which are stored in the supplied vector
 *
 * @note depending on how the code is compiled the function
 * will use openblas or simple self implemented routines
 *******************************************************************************************/
template<typename T>
Real l2Norm(const T& vector, const Int size)
{
    const Real* ptr = get_data_pointer_const(vector);
#ifdef VASPML_NO_CBLAS

    Real norm = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : norm)
#endif
    for (Size i = 0; i < (Size)size; i++) { norm += ptr[i] * ptr[i]; }
    return std::sqrt(norm);

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        IntDevice n = size;
        Real      result = (Real)0.0;
#ifdef VASPML_PALGO_GPU_NVIDIA
        cublasDnrm2(global::parallel.get_handle(), n, ptr, oneDevice, &result);
#endif // Variants of VASPML_PALGO_GPU
        return result;
    }
#endif // VASPML_PALGO_GPU

    IntHost n = size;

    return cblas_dnrm2(n, ptr, one);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * @brief Computes the outer product of two vectors and stores the result in a matrix.
 *
 * This function calculates the outer product of two vectors `vecA` and `vecB`, scales it by a
 *scalar `alpha`, and stores the result in the provided matrix `result`. The computation can
 *leverage optimized BLAS or cuBLAS libraries if enabled via preprocessor directives.
 *
 * @tparam T Type of the first input vector `vecA`.
 * @tparam U Type of the second input vector `vecB`.
 * @tparam V Type of the output matrix `result`.
 *
 * @param m Number of rows in the resulting matrix.
 * @param n Number of columns in the resulting matrix.
 * @param alpha Scalar multiplier for the outer product.
 * @param vecA Input vector of size `m`.
 * @param vecB Input vector of size `n`.
 * @param result Output matrix of size `m x n` to store the computed outer product.
 *
 * @note The function uses `cblas_dger` for BLAS-based computation if `VASPML_USE_CBLAS` is defined,
 *       or `cublasDscal` for cuBLAS-based computation if `VASPML_USE_CUBLAS` is defined.
 *       If neither is defined, the function does not perform any computation.
 *
 * @warning Ensure that the input vectors and output matrix have valid memory allocations
 *          and sufficient sizes to avoid undefined behavior.
 *******************************************************************************************/
template<typename T, typename U, typename V>
void outerProduct(const Int  m,
                  const Int  n,
                  const Real alpha,
                  const T&   vecA,
                  const U&   vecB,
                  V&         result)
{
#ifdef VASPML_NO_CBLAS

    global::tutor::error("This code path of VASPml requires CBLAS/LAPACKE (DGER).");

#else // VASPML_NO_CBLAS

    const Real* ptrA = get_data_pointer_const(vecA);
    const Real* ptrB = get_data_pointer_const(vecB);
    Real*       ptrResult = get_data_pointer(result);

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        //IntDevice m1 = m;
        //IntDevice n1 = n;
#ifdef VASPML_PALGO_GPU_NVIDIA
        global::tutor.error("This code path of VASPml requests outer product in CUBLAS on the "
                            "GPU which is not yet available.");
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    IntHost m1 = m;
    IntHost n1 = n;
    cblas_dger(CblasRowMajor, m1, n1, alpha, ptrA, one, ptrB, one, ptrResult, n1);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * rescale vector by a scalar
 *
 * @param vector  input/output vector which will be rescaled
 * @param factor  Real parameter by which the vector will be rescaled
 *
 * @note depending on how the code is compiled the function
 * will use openblas or simple self implemented routines
 *******************************************************************************************/
template<typename T>
void scaleVector(T& vector, const Real factor, const Int size)
{
    Real* ptr = get_data_pointer(vector);
#ifdef VASPML_NO_CBLAS

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size i = 0; i < (Size)size; i++) { ptr[i] *= factor; }

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        IntDevice n = size;
#ifdef VASPML_PALGO_GPU_NVIDIA
        cublasDscal(global::parallel.get_handle(), n, &factor, ptr, oneDevice);
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    IntHost n = size;
    cblas_dscal(n, factor, ptr, one);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * matrix matrix multiplication
 *
 * @param m   leading dimension of A and C
 * @param k   second dimension of A leading dimension of B
 * @param n   second dimension of B and C
 * @param A   @f[ \mathbf{A}\in \mathbb{R}^{ m\times k} @f]
 * @param B   @f[ \mathbf{B}\in \mathbb{R}^{ k\times n} @f]
 * @param C   @f[ \mathbf{C}\in \mathbb{R}^{ m\times n} @f]
 *
 * in detail the function is doing the following computation
 * @f[
 C _{ij} = \sum_{k\_indx=1}^{k}A_{i,k\_indx}B_{k\_indx,j} @f]
 *
 * @note depending on how the code is compiled the function
 * will use opencblas, intel mkl cblas or simple self implemented routines
 *
 * @warning when using the opencblas or intel mkl clbas library the matrix C has to be set to zero
 *  before calling matMul
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matMul(const Int m, const Int n, const Int k, const T& A, const U& B, V& C)
{
    const Real* A_ptr = get_data_pointer_const(A);
    const Real* B_ptr = get_data_pointer_const(B);
    Real*       C_ptr = get_data_pointer(C);
#ifdef VASPML_NO_CBLAS

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size m_indx = 0; m_indx < (Size)m; m_indx++)
    {
        for (Size n_indx = 0; n_indx < (Size)n; n_indx++)
        {
            Real tmp = 0;
            for (Size k_indx = 0; k_indx < (Size)k; k_indx++)
            {
                tmp += A_ptr[m_indx * k + k_indx] * B_ptr[k_indx * n + n_indx];
            }
            C_ptr[m_indx * n + n_indx] += tmp;
        }
    }

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        IntDevice      nn = n;
        IntDevice      kk = k;
        IntDevice      mm = m;
        constexpr Real one = (Real)1;
#ifdef VASPML_PALGO_GPU_NVIDIA
        cublasDgemm(global::parallel.get_handle(),
                    noTransDevice,
                    noTransDevice,
                    nn,
                    mm,
                    kk,
                    &one,
                    B_ptr,
                    nn,
                    A_ptr,
                    kk,
                    &one,
                    C_ptr,
                    nn);
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    IntHost        nn = n;
    IntHost        kk = k;
    IntHost        mm = m;
    constexpr Real one = (Real)1;
    cblas_dgemm(CblasRowMajor,
                noTransHost,
                noTransHost,
                mm,
                nn,
                kk,
                one,
                A_ptr,
                kk,
                B_ptr,
                nn,
                one,
                C_ptr,
                nn);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * matrix-matrix multiplication C := alpha*(op(A) + A_offset)*(op(B) + B_offset) + beta*C + C_offset
 *
 * op(X) is either op(X) = X or op(X) = XT,
 * @param transA if transA=CblasNoTrans, then op(A) = A; if transA=CblasTrans, then op(A) = AT.
 * @param transB if transb=CblasNoTrans, then op(B) = B; if transb=CblasTrans, then op(B) = BT.
 * @param m Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m
 *must be at least zero.
 * @param n Specifies the number of columns of the matrix op(B) and the number of columns of the
 * matrix C.
 *          The value of n must be at least zero.
 * @param k Specifies the number of columns of the matrix op(A) and the number of rows of the matrix
 * op(B).
 *          The value of k must be at least zero.
 * @param alpha Specifies the scalar alpha.
 * @param A
 * <table>
 * <caption id="description matrix A">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">Array, size lda*m; Before entry, the leading k-by-m part of the array a must
 * contain the matrix A.
 *     <td>Array, size lda*k; Before entry, the leading m-by-k part of the array a must contain the
 * matrix A.
 * </table>
 * @param lda
 * <table>
 * <caption id="description lda">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">lda must be at least max(1, k). <td>lda must be at least max(1, m).
 * </table>
 * @param B
 * <table>
 * <caption id="description matrix B">Complex table</caption>
 * <tr><th>transB=CblasNoTrans <th>transB=CblasTrans
 * <tr><td rowspan="2">Array, size ldb by k; Before entry, the leading n-by-k part of the array b
 * must contain the matrix B
 *     <td>Array, size ldb by n; Before entry, the leading k-by-n part of the array b must contain
 * the matrix B
 * </table>
 * @param ldb
 * <table>
 * <caption id="description ldb">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">ldb must be at least max(1, n).<td>ldb must be at least max(1, k).
 * </table>
 * @param beta Specifies the scalar beta. When beta is equal to zero, then c need not be set on
 * input.
 * @param C Array, size ldc by m. Before entry, the leading n-by-m part of the array C must contain
 * the matrix C, *except when beta is equal to zero, in which case C need not be set on entry.
 * @param ldc must be at least max(1, n).
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matMul(Transpose  transA,
            Transpose  transB,
            const Int  m,
            const Int  n,
            const Int  k,
            const Real alpha,
            const T&   A,
            Int        lda,
            const U&   B,
            Int        ldb,
            Real       beta,
            V&         C,
            Int        ldc)
{
#ifdef VASPML_NO_CBLAS

    global::tutor::error("This code path of VASPml requires CBLAS/LAPACKE (DGEMM).");

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        const Real*     A_ptr = get_data_pointer_const(A);
        const Real*     B_ptr = get_data_pointer_const(B);
        Real*           C_ptr = get_data_pointer(C);
        IntDevice       nn = n;
        IntDevice       kk = k;
        IntDevice       mm = m;
        IntDevice       lda_blas = lda;
        IntDevice       ldb_blas = ldb;
        IntDevice       ldc_blas = ldc;
        TransposeDevice transADevice = toDevice(transA);
        TransposeDevice transBDevice = toDevice(transB);
#ifdef VASPML_PALGO_GPU_NVIDIA
        // here the parameters have to be jumbled slightly to get the column major
        // format to which is needed for the cublas calls
        cublasDgemm(global::parallel.get_handle(),
                    transBDevice,
                    transADevice,
                    nn,
                    mm,
                    kk,
                    &alpha,
                    B_ptr,
                    ldb_blas,
                    A_ptr,
                    lda_blas,
                    &beta,
                    C_ptr,
                    ldc_blas);
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    const Real*   A_ptr = get_data_pointer_const(A);
    const Real*   B_ptr = get_data_pointer_const(B);
    Real*         C_ptr = get_data_pointer(C);
    IntHost       nn = n;
    IntHost       kk = k;
    IntHost       mm = m;
    IntHost       lda_blas = lda;
    IntHost       ldb_blas = ldb;
    IntHost       ldc_blas = ldc;
    TransposeHost transAHost = toHost(transA);
    TransposeHost transBHost = toHost(transB);
    cblas_dgemm(CblasRowMajor,
                transAHost,
                transBHost,
                mm,
                nn,
                kk,
                alpha,
                A_ptr,
                lda_blas,
                B_ptr,
                ldb_blas,
                beta,
                C_ptr,
                ldc_blas);
#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * @brief Computes the product of a matrix transpose and itself: A^T * A.
 *
 * This function performs the matrix multiplication `ATA = alpha * (A^T * A) + beta * ATA`,
 * where `A` is an `n x m` matrix, `ATA` is an `m x m` matrix, and `alpha` and `beta` are scalars.
 *
 * @tparam T Type of the input matrix `A`.
 * @tparam U Type of the output matrix `ATA`.
 * @param n Number of rows in the input matrix `A`.
 * @param m Number of columns in the input matrix `A`.
 * @param alpha Scalar multiplier for the product `A^T * A`.
 * @param A Input matrix of size `n x m`.
 * @param beta Scalar multiplier for the existing values in `ATA`.
 * @param ATA Output matrix of size `m x m` to store the result.
 *******************************************************************************************/
template<typename T, typename U>
void computeATransA(const Int n, const Int m, const Real alpha, const T& A, const Real beta, U& ATA)
{
    linalg::matMul(Transpose::TRANS, Transpose::NOTRANS, m, m, n, alpha, A, m, A, m, beta, ATA, m);
}
/*******************************************************************************************
 *
 * rescale vector by factor and add this to another vector
 *
 * @param factor rescaling factor of first vector
 * @param scaleVector vector to be rescaled by factor
 * @param addVector vector to which the rescaled vector is added. Contains result
 *
 * @f[addVector[i]= factor * scaleVector[i] + addVector[i] @f]
 *
 *******************************************************************************************/
template<typename T, typename U>
void scaleVectorPlusVector(const Real factor, const U& scaleVector, T& addVector, const Int size)
{
    const Real* scaleVector_ptr = get_data_pointer_const(scaleVector);
    Real*       addVector_ptr = get_data_pointer(addVector);
#ifdef VASPML_NO_CBLAS

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size i = 0; i < (Size)size; i++) { addVector_ptr[i] += factor * scaleVector_ptr[i]; }

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        IntDevice n = size;
#ifdef VASPML_PALGO_GPU_NVIDIA
        cublasDaxpy(global::parallel.get_handle(),
                    n,
                    &factor,
                    scaleVector_ptr,
                    oneDevice,
                    addVector_ptr,
                    oneDevice);
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    IntHost n = size;
    cblas_daxpy(n, factor, scaleVector_ptr, one, addVector_ptr, one);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * compute dot product of two input vectors
 *
 * @param vectorA first input vector
 * @param vectorB second input vector has to match size of vectorA
 * @param size the number of elements of type T stored in vectorA and VectorB
 *
 * @f[dotProduct=\sum_{i=1}^{N} vectorA[i] * vectorB[i] @f]
 *
 * @note vectorA and vectorB have to be of the same size
 *******************************************************************************************/
template<typename T, typename U>
Real dotProduct(const T& vectorA, const U& vectorB, const Int size)
{

    const Real* vectorA_ptr = get_data_pointer_const(vectorA);
    const Real* vectorB_ptr = get_data_pointer_const(vectorB);
#ifdef VASPML_NO_CBLAS

    Real norm = (Real)0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : norm)
#endif
    for (Size i = 0; i < (Size)size; i++) { norm += vectorA_ptr[i] * vectorB_ptr[i]; }
    return norm;

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        IntDevice n = size;
        Real      norm = 0;
#ifdef VASPML_PALGO_GPU_NVIDIA
        cublasDdot(global::parallel.get_handle(),
                   n,
                   vectorA_ptr,
                   oneDevice,
                   vectorB_ptr,
                   oneDevice,
                   &norm);
#endif // Variants of VASPML_PALGO_GPU
        return norm;
    }
#endif // VASPML_PALGO_GPU

    IntHost n = size;
    return cblas_ddot(n, vectorA_ptr, one, vectorB_ptr, one);

#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * matrix-Vector multiplication y := alpha*op(A)*x + beta*y
 *
 * op(X) is either op(X) = X or op(X) = XT,
 * @param trans if trans=CblasNoTrans, then op(A) = A; if transA=CblasTrans, then op(A) = AT.
 * @param m Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m
 *must be at least zero.
 * @param n Specifies the number of columns of the matrix op(A) and the number of rows of the vector
 *x
 * @param alpha Specifies the scalar alpha.
 * @param A
 * <table>
 * <caption id="description matrix A">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">Array, size lda*m; Before entry, the leading n-by-m part of the array a must
 * contain the matrix A.
 *     <td>Array, size lda*n; Before entry, the leading m-by-n part of the array a must contain the
 * matrix A.
 * </table>
 * @param lda
 * <table>
 * <caption id="description lda">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">lda must be at least max(1, n). <td>lda must be at least max(1, n).
 * </table>
 * @param x vector which is multiplied by A
 * @param incx On entry, incx specifies the increment for the elements of x. incx must not be zero.
 * @param beta Specifies the scalar beta. When beta is equal to zero, then y need not be set on
 *input.
 * @param incy On entry, incy specifies the increment for the elements of y. incx must not be zero.
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matVec(Transpose  trans,
            const Int  m,
            const Int  n,
            const Real alpha,
            const T&   A,
            const Int  lda,
            const U&   x,
            const Int  incx,
            const Real beta,
            V&         y,
            const Int  incy)
{
    const Real* A_ptr = get_data_pointer_const(A);
    const Real* x_ptr = get_data_pointer_const(x);
    Real*       y_ptr = get_data_pointer(y);
#ifdef VASPML_NO_CBLAS

    global::tutor::error("This code path of VASPml requires CBLAS/LAPACKE (DGEMV).");

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU_NVIDIA
    if (global::parallel.gpu())
    {
        IntDevice nn = n;
        IntDevice mm = m;
        IntDevice ldaInt = lda;
        IntDevice incxInt = incx;
        IntDevice incyInt = incy;
#ifdef VASPML_PALGO_GPU_NVIDIA
        // Since NVIDIA only allows for column major format which is not used in this C++ code
        // the transpose statement has to be inverted.
        if (trans == Transpose::NOTRANS) trans = Transpose::TRANS;
        else trans = Transpose::NOTRANS;
        TransposeDevice trans2 = toDevice(trans);
        cublasDgemv(global::parallel.get_handle(),
                    trans2,
                    nn,
                    mm,
                    &alpha,
                    A_ptr,
                    ldaInt,
                    x_ptr,
                    incxInt,
                    &beta,
                    y_ptr,
                    incyInt);
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

    IntHost       nn = n;
    IntHost       mm = m;
    IntHost       ldaInt = lda;
    IntHost       incxInt = incx;
    IntHost       incyInt = incy;
    TransposeHost transHost = toHost(trans);

    cblas_dgemv(CblasRowMajor,
                transHost,
                mm,
                nn,
                alpha,
                A_ptr,
                ldaInt,
                x_ptr,
                incxInt,
                beta,
                y_ptr,
                incyInt);
#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * Computes y := A^T * b
 *
 * @param A matrix which is supplied in non-transformed form
 * @param n first dimension of A. Number of rows
 * @param m second dimension of A. Number of columns
 * @param b input vector which is multiplied by A^T. Has to be of size n
 * @param y is the output vector which contains the result of A^T*b. Has to be of size m
 *
 *******************************************************************************************/
template<typename T, typename U, typename V>
void compute_ATransTimesVector(const Int  n,
                               const Int  m,
                               const Real alpha,
                               const T&   A,
                               const U&   b,
                               const Real beta,
                               V&         y)
{
    linalg::matVec(Transpose::TRANS, n, m, alpha, A, m, b, 1, beta, y, 1);
}
/*******************************************************************************************
 * compute inverse of a given square matrix by the use of LU factorization
 *
 * @param matrix matrix that will be inverted. Matrix is overwritten on output
 * @param n0 number of rows in supplied matrix
 * @param n1 number of columns in supplied matrix.
 *
 * @note the supplied matrix has to be in row major format. The only useful format for
 * matrices
 *
 * Step 1: \n
 * The LAPACKE_dgetrf function in the LAPACKE library performs an LU decomposition of a given
 *matrix. LU decomposition factors a matrix AA into two matrices: LL (lower triangular matrix) and
 *UU (upper triangular matrix), such that A=L×UA=L×U. This decomposition is useful for solving
 *linear systems, inverting matrices, and calculating determinants. Step 2: \n The LAPACKE_dgetri
 *function computes the inverse of a matrix using the LU decomposition that was previously obtained
 *from LAPACKE_dgetrf. Here's a step-by-step explanation of what the algorithm does: The algorithm
 *first computes the inverse of the upper triangular matrix U. It then solves for the inverse of the
 *original matrix by considering the lower triangular matrix L. Specifically, it performs forward
 *and backward substitution to find the solution.
 *******************************************************************************************/
template<typename T>
void computeInverseLU(T& matrix, const Int n0, const Int n1)
{
#ifdef VASPML_NO_CBLAS

    global::tutor.error("This code path of VASPml requires CBLAS/LAPACKE (DGETRI/DGETRF).");

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        global::tutor.error("This code path of VASPml requests LAPACKE (DGETRI/DGETRF) on the "
                            "GPU which is not yet available.");
#ifdef VASPML_PALGO_GPU_NVIDIA
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

#ifdef __NEC__
    global::tutor.error("This code path of VASPml requires LAPACKE (DGETRI/DGETRF) which is not "
                        "provided by the NEC compiler.");
#else
    if (n0 != n1)
    {
        global::tutor.bug("ERROR: computeInverseLU(T& matrix, const Int n0, const Int n1):\n"
                          "Matrix must be square (n0 == n1).\n");
    }
    Real*    matrix_ptr = get_data_pointer(matrix);
    IntHost* ipiv = new int[n0];
    IntHost  lda = n0;
    IntHost  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda, ipiv);
    if (info > 0)
    {
        global::tutor.error(
            "ERROR: computeInverseLU(T& matrix, const Int n0, const Int n1):\n"
            "LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda, ipiv) returned error code "
            + std::to_string(info) + "\n");
    }
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n0, matrix_ptr, lda, ipiv);
    if (info > 0)
    {
        global::tutor.error(
            "ERROR: computeInverseLU(T& matrix, const Int n0, const Int n1):\n"
            "LAPACKE_dgetri(LAPACK_ROW_MAJOR, n0, matrix_ptr, lda, ipiv) returned error code "
            + std::to_string(info) + "\n");
    }
    delete[] ipiv;
#endif // __NEC__
#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * @brief Computes the LU factorization of a given matrix.
 *
 * This function performs LU factorization on the input matrix using the LAPACK routine
 *`LAPACKE_dgetrf`. The factorization decomposes the matrix into a product of a lower triangular
 *matrix (L) and an upper triangular matrix (U). The pivot indices are stored in the provided
 *`pivot` vector.
 *
 * @tparam T The type of the matrix (e.g., a container supporting data access via
 *`get_data_pointer`).
 * @param matrix The matrix to be factorized. The factorization is performed in-place.
 * @param n0 The number of rows in the matrix.
 * @param n1 The number of columns in the matrix.
 * @param pivot A vector to store the pivot indices resulting from the factorization.
 *
 * @throws std::runtime_error If the function is compiled with the NEC compiler, as LU factorization
 *is not implemented for NEC.
 * @throws std::runtime_error If the LAPACK routine `LAPACKE_dgetrf` returns a positive error code,
 *indicating a singular matrix.
 *******************************************************************************************/
template<typename T>
void computeLUFactorization(T& matrix, const Int n0, const Int n1, Vec1Int& pivot)
{
#ifdef VASPML_NO_CBLAS

    global::tutor.error("This code path of VASPml requires CBLAS/LAPACKE (DGETRI/DGETRF).");

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        global::tutor.error("This code path of VASPml requests LAPACKE (DGETRF) on the "
                            "GPU which is not yet available.");
#ifdef VASPML_PALGO_GPU_NVIDIA
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

#ifdef __NEC__
    global::tutor.error("This code path of VASPml requires LAPACKE (DGETRF) which is not "
                        "provided by the NEC compiler.");
#else
    Real*    matrix_ptr = get_data_pointer(matrix);
    IntHost* ipiv = new int[n0];
    IntHost  lda = n0;
    IntHost  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda, ipiv);
    if (info > 0)
    {
        global::tutor.error(
            "ERROR: computeLUFactorization(T& matrix, const Int n0, const Int n1, U& pivot):\n"
            "LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda, ipiv) returned error code "
            + std::to_string(info) + "\n");
    }
    pivot.assign(ipiv, ipiv + n0);
#endif // __NEC__
#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * @brief Computes the determinant of a square matrix using LU factorization.
 *
 * This function calculates the determinant of a square matrix by performing LU
 * factorization and then computing the product of the diagonal elements of the
 * upper triangular matrix (U). The determinant's sign is adjusted based on the
 * number of row swaps during the factorization process.
 *
 * @tparam T The type of the matrix (e.g., dense matrix type).
 * @param matrix The input matrix whose determinant is to be computed. Must be square.
 * @param n0 The number of rows in the matrix.
 * @param n1 The number of columns in the matrix.
 * @return The determinant of the matrix.
 *
 * @throws std::runtime_error If the input matrix is not square (n0 != n1).
 *
 * @note This implementation requires the `VASPML_USE_CBLAS` macro to be defined.
 * @note For LAPACK-based LU factorization, the pivot indices are 1-based.
 *******************************************************************************************/
template<typename T>
Real computeDeterminantSquareMat(T& matrix, const Int n0, const Int n1)
{
    if (n0 != n1)
    {
        global::tutor.bug("ERROR: " + flf(VASPML_FLF)
                          + " square matrix is needed. Otherwise determinant can't be computed.");
    }
    Vec1Int pivot;
    computeLUFactorization(matrix, n0, n1, pivot);
    Real* matrix_ptr = get_data_pointer(matrix);
    Real  det = 1.0;
    Int   numRowSwaps = 0;
    // Count the number of row swaps based on pivot array
    for (Int i = 0; i < n0; i++)
    {
        if (pivot[i] != i + 1)
        { // pivot indices in LAPACK are 1-based
            numRowSwaps++;
        }
    }
    // Product of diagonal elements of U
    for (Int i = 0; i < n0; i++) { det *= matrix_ptr[i * n0 + i]; }
    // Adjust sign based on the number of row swaps
    if (numRowSwaps % 2 != 0) { det = -det; }

    return det;
}

/*******************************************************************************************
 * @brief Computes the determinant of the product of a matrix and its transpose (A * A^T).
 *
 * This function calculates the determinant of the product of a given matrix `A` and its transpose
 *`A^T`.
 * 1. Computes the product A * A^T using matrix multiplication.
 * 2. Performs LU factorization on the resulting matrix.
 * 3. Computes the determinant from the LU factorization, adjusting the sign based on the number of
 *row swaps.
 *
 * @tparam T The type of the input matrix.
 * @param matrix The input matrix `A` of dimensions `n1 x n0`.
 * @param n0 The number of rows in the input matrix `A`.
 * @param n1 The number of columns in the input matrix `A`.
 * @return The square root of the determinant of A * A^T.
 *
 * @note The determinant is computed as the product of the diagonal elements of the LU factorized
 *matrix, with the sign adjusted based on the number of row swaps.
 *******************************************************************************************/
template<typename T>
Real computeDeterminantATimesATrans(T& matrix, const Int n0, const Int n1)
{
    Vec1Real ata(n0 * n0);
    linalg::matMul(Transpose::NOTRANS,
                   Transpose::TRANS,
                   n0,
                   n0,
                   n1,
                   1.0,
                   matrix,
                   n1,
                   matrix,
                   n1,
                   0.0,
                   ata,
                   n0);
    Vec1Int pivot;
    computeLUFactorization(ata, n0, n0, pivot);
    Real* ataPtr = get_data_pointer(ata);
    Real  det = 1.0;
    Int   numRowSwaps = 0;
    // Count the number of row swaps based on pivot array
    for (Int i = 0; i < n0; i++)
    {
        if (pivot[i] != i + 1)
        { // pivot indices in LAPACK are 1-based
            numRowSwaps++;
        }
    }
    // Product of diagonal elements of U
    for (Int i = 0; i < n0; i++) { det *= ataPtr[i * n0 + i]; }
    // Adjust sign based on the number of row swaps
    if (numRowSwaps % 2 != 0) { det = -det; }

    return det;
}
/*******************************************************************************************
 * @brief Computes the Singular Value Decomposition (SVD) of a given matrix.
 *
 * This function performs SVD on the input matrix using LAPACK's `dgesvd` routine.
 * It decomposes the matrix into three components: U (left singular vectors),
 * S (singular values), and VT (right singular vectors, transposed). The function
 * also computes the super-diagonal elements of the bidiagonal matrix.
 *
 * @tparam T The type of the input matrix. Typically, this should be a matrix-like
 *           container that supports data access via `get_data_pointer`.
 * @param matrix The input matrix to be decomposed. It should be stored in row-major order.
 * @param n0 The number of rows in the input matrix.
 * @param n1 The number of columns in the input matrix.
 * @param U A vector to store the left singular vectors. It will be resized to `m * m`,
 *          where `m` is the number of rows in the input matrix.
 * @param S A vector to store the singular values. It will be resized to `min(m, n)`,
 *          where `m` is the number of rows and `n` is the number of columns in the input matrix.
 * @param VT A vector to store the transposed right singular vectors. It will be resized to `n * n`.
 * @param superb A vector to store the super-diagonal elements of the bidiagonal matrix.
 *               It will be resized to `min(m, n) - 1`.
 *
 * @throws std::runtime_error If the function is compiled with the NEC compiler,
 *                            as SVD is not implemented for this platform.
 * @note This function requires LAPACK and CBLAS to be available and linked.
 *       It uses the LAPACK_ROW_MAJOR storage format for matrix data.
 * @warning If the SVD computation does not converge, an error message is logged,
 *          and the function returns without modifying the output vectors.
 *
 * @pre The input matrix must be compatible with LAPACK's `dgesvd` routine.
 *
 * @post The vectors `U`, `S`, `VT`, and `superb` are resized and populated with
 *       the results of the SVD computation.
 *
 * @example
 * @code
 * Vec1Real matrix =
 * Vec1Real U, S, VT, superb;
 * Int n0 =  number of rows
 * Int n1 =  number of columns
 * computeSVD(matrix, n0, n1, U, S, VT, superb);
 * // U, S, VT, and superb now contain the SVD results.
 * @endcode
 *******************************************************************************************/
template<typename T>
void computeSVD(T&        matrix,
                const Int n0,
                const Int n1,
                Vec1Real& U,
                Vec1Real& S,
                Vec1Real& VT,
                Vec1Real& superb)
{
#ifdef VASPML_NO_CBLAS

    global::tutor.error("This code path of VASPml requires CBLAS/LAPACKE (DGESVD).");

#else // VASPML_NO_CBLAS

#ifdef VASPML_PALGO_GPU
    if (global::parallel.gpu())
    {
        global::tutor.error("This code path of VASPml requests LAPACKE (DGESVD) on the "
                            "GPU which is not yet available.");
#ifdef VASPML_PALGO_GPU_NVIDIA
#endif // Variants of VASPML_PALGO_GPU
        return;
    }
#endif // VASPML_PALGO_GPU

#ifdef __NEC__
    global::tutor.error("This code path of VASPml requires LAPACKE (DGESVD) which is not "
                        "provided by the NEC compiler.");
#else
    Real*   matrix_ptr = get_data_pointer(matrix);
    IntHost m = n0;
    IntHost n = n1;
    // Allocate space for SVD results
    S.resize(std::min(m, n)); // Singular values
    Size sU = (Size)m * (Size)m;
    U.resize(sU); // Left singular vectors
    Size sVT = (Size)n * (Size)n;
    VT.resize(sVT);                    // Right singular vectors (transposed)
    superb.resize(std::min(m, n) - 1); // Super-diagonal elements

    // Compute SVD using LAPACK's dgesvd
    Int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,
                              'A',
                              'A',
                              m,
                              n,
                              matrix_ptr,
                              n,
                              S.data(),
                              U.data(),
                              m,
                              VT.data(),
                              n,
                              superb.data());
    if (info > 0)
    {
        global::tutor.error("ERROR: computeSVD(computeSVD(T& matrix,const Int n0, const Int n1, "
                            "Vec1Real& U, Vec1Real& S, Vec1Real& VT, Vec1Real& superb):\n."
                            "SVD did not converge! Error code "
                            + std::to_string(info) + "\n");
        return;
    }
#endif // __NEC__
#endif // VASPML_NO_CBLAS
}
/*******************************************************************************************
 * @brief Computes the pseudo-inverse of a matrix using Singular Value Decomposition (SVD).
 *
 * This function computes the pseudo-inverse of a given matrix using SVD. It is implemented
 * for environments where the CBLAS library is available. The pseudo-inverse is calculated
 * by decomposing the matrix into U, S, and VT components, and then modifying the singular
 * values (S) to compute the pseudo-inverse.
 *
 * @tparam T The type of the matrix. Typically, this is a container that supports data access
 *           and manipulation (e.g., a custom matrix class or a vector of vectors).
 * @param matrix The input matrix to compute the pseudo-inverse for. This matrix will be modified
 *               in-place to store the pseudo-inverse.
 * @param m The number of rows in the input matrix.
 * @param n The number of columns in the input matrix.
 * @param regParam Regularization parameter (default is 0.0).
 *
 * @note Singular values smaller than 1e-10 are treated as zero to avoid numerical instability.
 *******************************************************************************************/
template<typename T>
void computePseudoInverseSVD(T& matrix, const Int m, const Int n, Real regParam = 0.0)
{
    Vec1Real U;
    Vec1Real S;
    Vec1Real VT;
    Vec1Real superb;

    computeSVD(matrix, m, n, U, S, VT, superb);

    Real* matrix_ptr = get_data_pointer(matrix);

    Real regSquare = regParam * regParam;
    for (Int i = 0; i < std::min(m, n); i++)
    {
        Real factor = 0;
        if (S[i] > 1e-10)
        {
            //factor = 1.0 / S[i];
            factor = S[i] / (S[i] * S[i] + regSquare);
        }
        for (Int k = 0; k < n; k++) { VT[i * n + k] = VT[i * n + k] * factor; }
    }
    matMul(Transpose::NOTRANS,
           Transpose::NOTRANS,
           m,
           n,
           std::min(m, n),
           1.0,
           U.data(),
           m,
           VT.data(),
           n,
           0.0,
           matrix_ptr,
           n);
}
/*******************************************************************************************
 * @brief Solves a least squares problem using Singular Value Decomposition (SVD).
 *
 * This function computes the solution to a least squares problem of the form Ax = b,
 * where A is a matrix, b is a vector of observations, and x is the vector of unknowns.
 * The solution is computed using the pseudo-inverse of the matrix A obtained via SVD.
 *
 * @tparam T Type of the matrix. Must support operations required by `computePseudoInverseSVD` and
 *`get_data_pointer`.
 * @tparam U Type of the vector b. Must support operations required by `get_data_pointer_const`.
 * @tparam V Type of the vector x. Must support operations required by `get_data_pointer`.
 *
 * @param[in,out] matrix The input matrix A (n0 x n1). This matrix is modified in-place to store its
 *pseudo-inverse.
 * @param[in] n0 Number of rows in the matrix A.
 * @param[in] n1 Number of columns in the matrix A.
 * @param[in] b The input vector b (size n0). Represents the observations.
 * @param[in] regParam Regularization parameter (default is 0.0). This parameter can be used to
 *stabilize the solution in case of ill-conditioned matrices.
 * @param[out] x The output vector x (size n1). Represents the solution to the least squares
 *problem.
 *
 * @note This function relies on the following helper functions:
 *       - `computePseudoInverseSVD`: Computes the pseudo-inverse of the matrix using SVD.
 *       - `get_data_pointer`: Retrieves a raw pointer to the data of a container.
 *       - `get_data_pointer_const`: Retrieves a raw pointer to the constant data of a container.
 *       - `matVec`: Performs matrix-vector multiplication.
 *
 * @example
 * // Example usage:
 * Matrix<double> A; // Define and initialize matrix A
 * Vector<double> b; // Define and initialize vector b
 * Vector<double> x; // Define vector x to store the solution
 * Int n0 = A.rows();
 * Int n1 = A.cols();
 * solveLeastSquaresSVD(A, n0, n1, b, x);
 *******************************************************************************************/
template<typename T, typename U, typename V>
void solveLeastSquaresSVD(T&          matrix,
                          const Int   n0,
                          const Int   n1,
                          const U&    b,
                          V&          x,
                          const Real& regParam = 0.0)
{
    computePseudoInverseSVD(matrix, n0, n1, regParam);
    Real*       x_ptr = get_data_pointer(x);
    Real*       matrix_ptr = get_data_pointer(matrix);
    const Real* b_ptr = get_data_pointer_const(b);
    matVec(Transpose::TRANS, n0, n1, 1.0, matrix_ptr, n1, b_ptr, 1, 0.0, x_ptr, 1);
}

/*================================================================================================+
 |
 | Functions available only for GPU.
 |
 +================================================================================================*/

#ifdef VASPML_PALGO_GPU
/*******************************************************************************************
 * This routine can be used to copy data from the GPU back to the CPU
 *
 * This function copies n elements from a vector x in GPU memory space to a vector
 * y in host memory space. Elements in both vectors are assumed to have a size of elemSize
 * bytes. The storage spacing between consecutive elements is given by incx for the source
 * vector and incy for the destination vector y.
 *
 * @param[in] data1 data which is currently located on the GPU
 * @param[inout] data2 data to which the GPU data is copied to be on CPU
 * @param[in] n is the number of elements in data1 and data2
 *
 * @note data1 and data2 can be the same variable
 *******************************************************************************************/
template<typename T>
void getVectorGPU(const T& data1, T& data2, const Int n)
{
    const Real* data1Ptr = get_data_pointer_const(data1);
    Real*       data2Ptr = get_data_pointer(data2);
#ifdef VASPML_PALGO_GPU_NVIDIA
    cublasGetVector(n, sizeof(T), data1Ptr, one, data2Ptr, one);
#endif
}
/*******************************************************************************************
 * This routine can be used to copy data from the CPU to the GPU
 *
 * This function copies n elements from a vector x in host memory space to a vector y in GPU
 * memory space. Elements in both vectors are assumed to have a size of elemSize bytes.
 * The storage spacing between consecutive elements is given by incx for the source vector x
 * and by incy for the destination vector y.
 *
 * @param[in] data1 data which is currently located on the CPU
 * @param[inout] data2 data to which the CPU data is copied to be on GPU
 * @param[in] n is the number of elements in data1 and data2
 *
 * @note data1 and data2 can be the same variable
 *******************************************************************************************/
template<typename T>
void setVectorGPU(const T& data1, T& data2, const Int n)
{
    const Real* data1Ptr = get_data_pointer_const(data1);
    Real*       data2Ptr = get_data_pointer(data2);
#ifdef VASPML_PALGO_GPU_NVIDIA
    cublasGetVector(n, sizeof(T), data1Ptr, one, data2Ptr, one);
#endif
}
#endif // VASPML_PALGO_GPU

} // namespace vaspml::linalg

#endif
