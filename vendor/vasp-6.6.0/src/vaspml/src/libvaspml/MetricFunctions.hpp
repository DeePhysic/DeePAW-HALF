#ifndef METRICFUNCTIONS_HPP
#define METRICFUNCTIONS_HPP

#include "Linalg.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "types.hpp"

// Avoid IWYU suggesting <stdlib.h> because of integer argument variant of std::abs (not used here).
// IWYU pragma: no_include <stdlib.h>
#include <algorithm>    // for max
#include <cmath>        // for sqrt, abs, pow
#include <type_traits>  // for is_arithmetic_v

namespace vaspml
{
/*******************************************************************************************
 * @class Metrics
 * @brief A utility class for computing various mathematical norms and metrics between vectors.
 *
 * @tparam T The type of elements in the vectors. Must be an arithmetic type.
 *
 * This class provides static and non-static methods to compute norms and metrics such as L1 norm,
 * L2 norm, Euclidean distance, polynomial kernel norm, Pearson correlation, and their variations.
 * It ensures that the input vectors are of the same size and supports only arithmetic types.
 *******************************************************************************************/
template<typename T>
class MetricFunctions
{
    static_assert(std::is_floating_point<T>::value,
                  "Cannot specialize MetricFunctions to given template type, a floating-point type "
                  "is required.");

  public:
    /*******************************************************************************************
     * @brief Computes the L2 norm (Euclidean distance) between two arrays.
     *
     * This function calculates the L2 norm between two arrays of type `T` with a given size `n`.
     *
     * @tparam T The data type of the elements in the arrays.
     * @param x Pointer to the first array.
     * @param y Pointer to the second array.
     * @param n The number of elements in the arrays.
     * @return The L2 norm (Euclidean distance) between the two arrays.
     *******************************************************************************************/
    static constexpr T l2NormPtr(const T* x, const T* y, const Size n)
    {
        T norm = std::max(linalg::dotProduct(x, x, n) + linalg::dotProduct(y, y, n)
                              - 2.0 * linalg::dotProduct(x, y, n),
                          0.0);
        return std::sqrt(norm);
    }
    /*******************************************************************************************
     * @brief Computes the L2 norm (Euclidean distance) between two vectors.
     *
     * This function calculates the L2 norm between two vectors of type `T`. It ensures that the
     *sizes of the input vectors match before performing the computation.
     *
     * @tparam T The data type of the elements in the vectors.
     * @param x The first vector.
     * @param y The second vector.
     * @return The L2 norm (Euclidean distance) between the two vectors.
     * @note If the sizes of the input vectors do not match, a debug message is logged.
     *******************************************************************************************/
    static constexpr T l2Norm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return l2NormPtr(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the L2 norm (Euclidean distance) between two arrays which are normalized to
     *unity.
     *
     * This function calculates the L2 norm between two arrays of type `T` with a given size `n`.
     *
     * @tparam T The data type of the elements in the arrays.
     * @param x Pointer to the first array.
     * @param y Pointer to the second array.
     * @param n The number of elements in the arrays.
     * @return The L2 norm (Euclidean distance) between the two arrays.
     *******************************************************************************************/
    static constexpr T l2NormPtrUnity(const T* x, const T* y, const Size n)
    {
        T norm = std::max(2.0 - 2.0 * linalg::dotProduct(x, y, n), 0.0);
        return std::sqrt(norm);
    }
    /*******************************************************************************************
     * @brief Computes the L2 norm (Euclidean distance) between two vectors.
     *
     * This function calculates the L2 norm between two vectors of type `T`. It ensures that the
     *sizes of the input vectors match before performing the computation.
     *
     * @tparam T The data type of the elements in the vectors.
     * @param x The first vector.
     * @param y The second vector.
     * @return The L2 norm (Euclidean distance) between the two vectors.
     * @note If the sizes of the input vectors do not match, a debug message is logged.
     *******************************************************************************************/
    static constexpr T l2NormUnity(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return l2NormPtrUnity(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the Euclidean distance between two points represented as raw pointers.
     *
     * This function calculates the Euclidean distance between two points `p1` and `p2` in an
     *`n`-dimensional space. The points are represented as raw pointers to arrays of type `T`.
     *
     * @tparam T The data type of the elements (e.g., `float`, `double`).
     * @param p1 Pointer to the first point (array of size `n`).
     * @param p2 Pointer to the second point (array of size `n`).
     * @param n The number of dimensions (size of the arrays).
     * @return The Euclidean distance between the two points.
     *******************************************************************************************/
    T euclideanDistancePtr(const T* p1, const T* p2, const Size n)
    {
        T norm = std::max(2.0 - 2.0 * linalg::dotProduct(p1, p2, n), 0.0);
        return std::sqrt(norm);
    }
    /*******************************************************************************************
     * @brief Computes the Euclidean distance between two vectors.
     *
     * This function calculates the Euclidean distance between two vectors `p1` and `p2`.
     * It ensures that the input vectors have the same size and delegates the computation
     * to the `euclideanDistancePtr` function.
     *
     * @tparam T The type of the elements in the vectors.
     * @param p1 The first input vector.
     * @param p2 The second input vector.
     * @return The Euclidean distance between the two vectors.
     *
     * @note If the sizes of `p1` and `p2` do not match, a debug-level error message is logged.
     *******************************************************************************************/
    T euclideanDistance(const std::vector<T>& p1, const std::vector<T>& p2)
    {
        VASPML_DEBUG_L1(
            if (p1.size() != p2.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return euclideanDistancePtr(p1.data(), p2.data(), p1.size());
    }
    /*******************************************************************************************
     * @brief Computes the L1 norm (Manhattan distance) between two arrays.
     *
     * This function calculates the L1 norm (sum of absolute differences) between
     * two arrays of type `T` with a given size `n`.
     *
     * @tparam T The data type of the elements in the arrays.
     * @param x Pointer to the first array.
     * @param y Pointer to the second array.
     * @param n The number of elements in the arrays.
     * @return The L1 norm (Manhattan distance) between the two arrays.
     *******************************************************************************************/
    static constexpr T l1NormPtr(const T* x, const T* y, const Size n)
    {
        T norm = 0;
        for (Size i = 0; i < n; i++) { norm += std::abs(x[i] - y[i]); }
        return norm;
    }
    /*******************************************************************************************
     * @brief Computes the L1 norm (Manhattan distance) between two vectors.
     *
     * This function calculates the L1 norm (sum of absolute differences) between
     * two vectors of type `T`. It validates that the sizes of the input vectors
     * match before performing the computation.
     *
     * @tparam T The data type of the elements in the vectors.
     * @param x The first input vector.
     * @param y The second input vector.
     * @return The L1 norm (Manhattan distance) between the two vectors.
     * @note If the sizes of the input vectors do not match, a debug message is
     *       logged, and the behavior is undefined.
     *******************************************************************************************/
    static constexpr T l1Norm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return l1NormPtr(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the polynomial kernel normalization between two vectors using raw pointers.
     *
     * This function calculates the polynomial kernel normalization between two vectors `x` and `y`
     * of size `n` using the formula:
     *
     * \f[
     * \text{polyKernelNorm} = \text{intPowMixAlgo}(\|x\|^2, p) + \text{intPowMixAlgo}(\|y\|^2, p)
     * - 2 \cdot \text{intPowMixAlgo}(\langle x, y \rangle, p)
     * \f]
     *
     * where \f$\|x\|^2\f$ and \f$\|y\|^2\f$ are the squared norms of `x` and `y`,
     * and \f$\langle x, y \rangle\f$ is the dot product of `x` and `y`.
     *
     * @tparam p The integer power used in the polynomial kernel computation.
     * @param x Pointer to the first input vector.
     * @param y Pointer to the second input vector.
     * @param n The size of the input vectors.
     * @return The computed polynomial kernel normalization value.
     *******************************************************************************************/
    template<Int p>
    static constexpr T polyKernelNormPtr(const T* x, const T* y, const Size n)
    {
        T dotPaa = linalg::dotProduct(x, x, n);
        T dotPbb = linalg::dotProduct(y, y, n);
        T dotPab = linalg::dotProduct(x, y, n);
        return math::intPowMixAlgo(dotPaa, p) + math::intPowMixAlgo(dotPbb, p)
             - 2.0 * math::intPowMixAlgo(dotPab, p);
    }
    /*******************************************************************************************
     * @brief Computes the polynomial kernel norm of two vectors.
     *
     * This function calculates the polynomial kernel norm of two input vectors `x` and `y`
     * using the specified polynomial degree `p`. It validates that the sizes of the input
     * vectors match before proceeding with the computation.
     *
     * @tparam p The polynomial degree used in the kernel computation.
     * @param x The first input vector.
     * @param y The second input vector.
     * @return The computed polynomial kernel norm.
     *
     * @note This function assumes that the input vectors are of the same size.
     *       If the sizes do not match, a debug-level error message is logged.
     *******************************************************************************************/
    template<Int p>
    static constexpr T polyKernelNorm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return polyKernelNormPtr<p>(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the Lp-norm of two arrays.
     *
     * This function calculates the Lp-norm of two input arrays `x` and `y` of size `n`
     * using the specified power `p`. The Lp-norm is computed as the p-th root of the
     * sum of the p-th powers of the absolute differences between corresponding elements.
     *
     * @tparam p The power used in the Lp-norm computation.
     * @param x Pointer to the first input array.
     * @param y Pointer to the second input array.
     * @param n The size of the input arrays.
     * @return The computed Lp-norm.
     *
     * @note This function assumes that the input arrays have the same size `n`.
     *******************************************************************************************/
    template<int p>
    static constexpr T lpNormPtr(const T* x, const T* y, const Size n)
    {
        T norm = 0;
        for (Size i = 0; i < n; i++) { norm += math::intPowMixAlgo(x[i] - y[i], p); }
        return std::pow(norm, 1.0 / (T)p);
    }
    /*******************************************************************************************
     * @brief Computes the Lp norm of the difference between two vectors.
     *
     * This function calculates the Lp norm (parameterized by `p`) of the difference
     * between two input vectors `x` and `y`. The Lp norm is a generalization of
     * vector norms, where `p` determines the type of norm (e.g., p=2 for Euclidean norm).
     *
     * @tparam p The order of the norm (e.g., 1 for Manhattan norm, 2 for Euclidean norm).
     * @tparam T The type of the elements in the vectors.
     * @param x The first input vector.
     * @param y The second input vector.
     * @return The computed Lp norm of the difference between `x` and `y`.
     *
     * @note The function assumes that `lpNormPtr` is a valid function that computes
     *       the Lp norm given raw pointers to the data and the size of the vectors.
     * @note If the sizes of `x` and `y` do not match, a debug-level error message
     *       is logged (if debugging is enabled).
     *******************************************************************************************/
    template<int p>
    static constexpr T lpNorm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return lpNormPtr(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the Pearson correlation coefficient for two arrays of data.
     *
     * This function calculates the Pearson correlation coefficient between two
     * arrays of data using pointers. The Pearson correlation coefficient measures
     * the linear relationship between two datasets.
     *
     * @tparam T The data type of the input arrays.
     * @param x Pointer to the first array of data.
     * @param y Pointer to the second array of data.
     * @param n The number of elements in the arrays.
     * @return The Pearson correlation coefficient as a value of type T.
     *******************************************************************************************/
    static constexpr T pearsonNormPtr(const T* x, const T* y, const Size n)
    {
        T r = math::stats::pearsonCorrPtr(x, y, n);
        return r;
    }
    /*******************************************************************************************
     * @brief Computes the Pearson correlation coefficient for two vectors of data.
     *
     * This function calculates the Pearson correlation coefficient between two
     * vectors of data. It includes a debug check to ensure that the sizes of the
     * input vectors match. If the sizes do not match, an error message is logged.
     *
     * @tparam T The data type of the input vectors.
     * @param x The first vector of data.
     * @param y The second vector of data.
     * @return The Pearson correlation coefficient as a value of type T.
     *
     * @note This function assumes that the input vectors have the same size.
     *       If the sizes differ, a debug error message is generated.
     *******************************************************************************************/
    static constexpr T pearsonNorm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        T r = math::stats::pearsonCorr(x, y);
        return r;
    }
    /*******************************************************************************************
     * @brief Computes the squared Pearson correlation coefficient between two arrays.
     *
     * @tparam T The data type of the input arrays.
     * @param x Pointer to the first array.
     * @param y Pointer to the second array.
     * @param n The size of the arrays.
     * @return The squared Pearson correlation coefficient.
     *******************************************************************************************/
    static constexpr T pearsonNorm2Ptr(const T* x, const T* y, const Size n)
    {
        T r = math::stats::pearsonCorrPtr(x, y, n);
        return r * r;
    }
    /*******************************************************************************************
     * @brief Computes the squared Pearson correlation coefficient between two vectors.
     *
     * @tparam T The data type of the input vectors.
     * @param x The first vector.
     * @param y The second vector.
     * @return The squared Pearson correlation coefficient.
     * @note Throws an error if the sizes of the input vectors do not match.
     *******************************************************************************************/
    static constexpr T pearsonNorm2(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        T r = math::stats::pearsonCorr(x, y);
        return r * r;
    }
    /*******************************************************************************************
     * @brief Computes the anti-Pearson correlation coefficient (1 - Pearson correlation) between
     *two arrays.
     *
     * @tparam T The data type of the input arrays.
     * @param x Pointer to the first array.
     * @param y Pointer to the second array.
     * @param n The size of the arrays.
     * @return The anti-Pearson correlation coefficient.
     *******************************************************************************************/
    static constexpr T antiPearsonNormPtr(const T* x, const T* y, const Size n)
    {
        T r = math::stats::pearsonCorrPtr(x, y, n);
        return 1.0 - r;
    }
    /*******************************************************************************************
     * @brief Computes the anti-Pearson normalization between two vectors.
     *
     * This function calculates the anti-Pearson normalization for two input vectors `x` and `y`.
     * It ensures that the sizes of the input vectors match before proceeding with the computation.
     *
     * @tparam T The data type of the elements in the vectors.
     * @param x The first input vector.
     * @param y The second input vector.
     * @return The computed anti-Pearson normalization value.
     *
     * @note If the sizes of `x` and `y` do not match, a debug-level error message is logged.
     * @warning This function assumes that the input vectors are valid and non-empty.
     *******************************************************************************************/
    static constexpr T antiPearsonNorm(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        return antiPearsonNormPtr(x.data(), y.data(), x.size());
    }
    /*******************************************************************************************
     * @brief Computes the squared anti-Pearson correlation coefficient for two data arrays.
     *
     * This function calculates the squared anti-Pearson correlation coefficient (1 - r^2)
     * for two input arrays of size `n`. The Pearson correlation coefficient is computed
     * using a pointer-based implementation.
     *
     * @tparam T The data type of the input arrays (e.g., float, double).
     * @param x Pointer to the first input array.
     * @param y Pointer to the second input array.
     * @param n The size of the input arrays.
     * @return The squared anti-Pearson correlation coefficient (1 - r^2).
     *******************************************************************************************/
    static constexpr T antiPearsonNorm2Ptr(const T* x, const T* y, const Size n)
    {
        T r = math::stats::pearsonCorrPtr(x, y, n);
        return 1.0 - r * r;
    }
    /*******************************************************************************************
     * @brief Computes the squared anti-Pearson correlation coefficient for two vectors.
     *
     * This function calculates the squared anti-Pearson correlation coefficient (1 - r^2)
     * for two input vectors. It ensures that the sizes of the input vectors match before
     * performing the computation. The Pearson correlation coefficient is computed using
     * a vector-based implementation.
     *
     * @tparam T The data type of the input vectors (e.g., float, double).
     * @param x The first input vector.
     * @param y The second input vector.
     * @return The squared anti-Pearson correlation coefficient (1 - r^2).
     * @throws std::runtime_error If the sizes of the input vectors do not match.
     *******************************************************************************************/
    static constexpr T antiPearsonNorm2(const std::vector<T>& x, const std::vector<T>& y)
    {
        VASPML_DEBUG_L1(
            if (x.size() != y.size())
                global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                                  + "Sizes of input vectors x and y don't match\n");
        );
        T r = math::stats::pearsonCorr(x, y);
        return 1.0 - r * r;
    }
};

} // namespace vaspml
#endif // METRICFUNCTIONS_HPP
