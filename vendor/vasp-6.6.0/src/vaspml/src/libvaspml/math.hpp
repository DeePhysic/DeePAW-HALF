#ifndef MATH_HPP
#define MATH_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm>    // for transform, count, max_element, min_element
#include <cassert>      // for assert
#include <cmath>        // for sqrt, ceil, pow, floor, abs, isnan
#include <functional>   // for multiplies, _1, bind, plus
#include <limits>       // for numeric_limits
#include <numeric>      // for accumulate, inner_product, partial_sum
#include <type_traits>  // for is_integral, is_arithmetic

namespace vaspml::math
{

/** Fills the given output vector with factorial numbers.
 *
 * @param f Factorial numbers output vector (determines maximum integer, see note).
 *
 * @note The size of the output vector on entry determines the maximum integer
 * number to which factorials should be computed. For factorials up to
 * @f$N_\max@f$ reserve space for one more entry in the output vector, i.e.
 * f.size() = @f$N_\max + 1@f$. The first entry is 0! = 1.
 */
void factorial(Vec1Real& f);
/** Fills the given output vector with double factorials for odd numbers.
 *
 * @param f Double factorial numbers output vector (determines maximum integer, see note).
 *
 * @note The size of the output vector on entry determines the maximum integer number to
 * which double factorials should be computed. To map between the output
 * vector index i and N use the expressions @f$N = 2i + 1@f$ and @f$i = (N - 1)/2@f$.
 * Hence, the size of the output vector should be f.size() = @f$(N_\max + 1) / 2@f$.
 *
 * Example: get N!! up to N = 7 requires a vector with size 4:
 *
 * * i = 0 => N = 1 => N!! = 1
 * * i = 1 => N = 3 => N!! = 3
 * * i = 2 => N = 5 => N!! = 15
 * * i = 3 => N = 7 => N!! = 105
 */
void doubleFactorialOdd(Vec1Real& f);
/** Computes the Clebsch-Gordan coefficient @f$\left<j_1 m_1 j_2 m_2|J M\right>@f$ using Racah
 * recursion relations.
 *
 * @param j1 Index @f$j_1@f$.
 * @param j2 Index @f$j_2@f$.
 * @param J Index @f$J@f$.
 * @param m1 Index @f$m_1@f$.
 * @param m2 Index @f$m_2@f$.
 * @param M Index @f$M@f$.
 * @param factorials Precomputed factorial numbers.
 * @return Clebsch-Gordan coefficient @f$\left<j_1 m_1 j_2 m_2|J M\right>@f$.
 *
 * @note Precomputed factorial numbers need to be passed to this function, use factorial() for
 * this.
 *
 * @todo Find out how many factorials we actually need, default in Fortran code up to 40! makes
 * little sense. I think 4*jmax + 2 should be sufficient. Applies to clebschGordanM0() as well.
 * Alternatively or additionally add debugging checks.
 *
 * @note It seems the original Fortran implementation was not designed to be passed parameters which
 * do not fullfill the triangular condition @f$|j_1 - j_2| \leq J \leq j_1 + j_2@f$ because the
 * variables `k1`, `k2`, and `k3` will cause an out-of-bounds access in the factorials array. This
 * implementation provides additional safeguards which `return 0.0` before accessing the array in
 * these cases.
 */
Real clebschGordan(Int j1, Int j2, Int J, Int m1, Int m2, Int M, const Vec1Real& factorials);
/** Computes the Clebsch-Gordan coefficient @f$\left<j_1 0 j_2 0|J 0\right>@f$ using Racah
 * recursion relations.
 *
 * @param j1 Index @f$j_1@f$.
 * @param j2 Index @f$j_2@f$.
 * @param J Index @f$J@f$.
 * @param factorials Precomputed factorial numbers.
 * @return Clebsch-Gordan coefficient @f$\left<j_1 0 j_2 0|J 0\right>@f$.
 *
 * @note Precomputed factorial numbers need to be passed to this function, use factorial() for
 * this.
 */
Real clebschGordanM0(Int j1, Int j2, Int J, const Vec1Real& factorials);
/** Computes Gaussian function for a vector of numbers.
 *
 * @param x x-values for which Gaussian should be computed.
 * @param f Gaussian values for given x-values.
 * @param width Gaussian width parameter.
 *
 * Gaussian is assumed to be centered around zero.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void gaussian(const Vec1Real& x, Vec1Real& f, Real width = 1.0);
/** Computes Gaussian function and its derivative for a vector of numbers.
 *
 * @param x x-values for which Gaussian should be computed.
 * @param f Gaussian values for given x-values.
 * @param df Gaussian derivative for given x-values.
 * @param width Gaussian width parameter.
 *
 * Gaussian is assumed to be centered around zero.
 *
 * @warning Size of output vectors must be equal to input vector on entry.
 */
void gaussianAndDerivative(const Vec1Real& x, Vec1Real& f, Vec1Real& df, Real width = 1.0);
/** Iteratively Computes spherical Bessel function of the first kind for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 *
 * @note This function gets automatically called by #sphericalBessel if the argument is below 1.0.
 */
Real sphericalBesselIterative(const Real x, const UInt nu);
/** Computes spherical Bessel function of the first kind for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 */
Real sphericalBessel(const Real x, const UInt nu);
/** Computes spherical Bessel function of the first kind and its derivative for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param jnu Output Bessel function value.
 * @param djnu Output Bessel function derivative.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 */
void sphericalBesselAndDerivative(const Real x, Real& jnu, Real& djnu, const UInt nu);
/** Computes spherical Bessel function of the first kind for a vector of numbers.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param jnu Bessel function values for given x-values.
 * @param nu Order of the Bessel function.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void sphericalBessel(const Vec1Real& x, Vec1Real& jnu, const UInt nu);
/** Computes spherical Bessel function of the first kind and its derivative for a vector of numbers.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param jnu Bessel function values for given x-values.
 * @param djnu Bessel function derivatives for given x-values.
 * @param nu Order of the Bessel function.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void sphericalBesselAndDerivative(const Vec1Real& x, Vec1Real& jnu, Vec1Real& djnu, const UInt nu);
/** Computes multiple roots of a spherical Bessel function of the first kind.
 *
 * @param roots x-values for which the Bessel function is zero, in ascending order (see note).
 * @param nu Order of the Bessel function.
 *
 * @note The size of the output vector on entry to this function determines
 * how many roots are calculated.
 */
void sphericalBesselRoots(Vec1Real& roots, const UInt nu);
/** Computes modified spherical Bessel function of the first kind and its derivative.
 *
 * @warning The output vectors' first dimensions determine the maximum order of the Bessel
 * functions. Hence, on entry the condition `inu.size()` = @f$l_\max + 1@f$ must be satisfied. The
 * second dimension must match with the size of the input vector. In summary,
 * * index 1: Bessel function orders @f$l = 0,\ldots,l_\max@f$.
 * * index 2: Corresponds to input vector x-values.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param inu Bessel function values for given x-values.
 * @param dinu Bessel function derivatives for given x-values.
 */
void modifiedSphericalBesselAndDerivative(const Vec1Real& x, Vec2Real& inu, Vec2Real& dinu);
/** Computes product of exponential and modified spherical Bessel function of the first kind.
 *
 * @param a2 Input scalar @f$a^2@f$ value.
 * @param ab Input vector @f$ab@f$, i.e., should be a product of scalar @f$a@f$ and vector @f$b@f$.
 * @param b2 Input vector @f$b^2@f$.
 * @param p Output product of the exponential and Bessel function expression.
 *
 * This function is used to compute parts of @f$h_{nl}@f$ found in expression
 * (19) in Phys. Rev. B 100, 014105 (see also https://arxiv.org/pdf/1904.12961.pdf).
 * Given the entire expression
 * @f[
 * h_{nl}\left(r\right) = \frac{4\pi}{\left(\sqrt{2\sigma_\mathrm{atom}^{2}\pi}\right)^{3}}
 * f_\mathrm{cut}\left(r\right)
 * \int  _{0}^{\infty} \chi_{nl}\left(r'\right)
 * \mathrm{exp}\left(-\frac{r^{2}+r'^{2}}{2\sigma_\mathrm{atom}^{2}}\right)
 * \iota_{l}\left(\frac{rr'}{\sigma_{\mathrm  {atom}}^{2}}\right)
 * r'^{2}dr'
 * @f]
 * this function computes the product of the exponential and modified spherical
 * Bessel function (first kind), i.e.,
 * @f[
 * \mathrm{exp}\left(-a^2-b^2\right) \iota_{l}\left(ab\right),
 * @f]
 * where @f$a^2@f$ is passed as a scalar value and @f$ab@f$ and @f$b^2@f$ are passed as vectors of
 * values on a grid. The function computes the product above for all points in the vectors @f$ab@f$
 * and @f$b^2@f$ and for all Bessel function orders up to a given maximum integer number.
 *
 * The reason for computing this product instead of the individual terms is that the modified
 * spherical Bessel function alone quickly diverges for large arguments. However, it can be
 * expressed in terms of hyperbolic sine and cosine, i.e., exponential functions. Multiplying with
 * the exponential function @f$\mathrm{exp}\left(-a^2-b^2\right)@f$ results in terms of
 * @f$\mathrm{exp}\left(-a^2-b^2 \pm ab\right)@f$ which allows the arguments to cancel before the
 * exponentiation is carried out.
 *
 * Expressions for small @f$ab@f$ values are derived from the limit of the modified spherical Bessel
 * function for small arguments:
 * @f[
 * \iota_{l}(x) \approx \frac{x^l}{(2l + 1)!!}
 * @f]
 *
 * @warning The output vector's first dimension determines the maximum order of the Bessel
 * functions. Hence, on entry the condition `p.size()` = @f$l_\max + 1@f$ must be satisfied. The
 * second dimension must match with the size of the input vectors. In summary,
 * * index 1: Bessel function orders @f$l = 0,\ldots,l_\max@f$.
 * * index 2: Corresponds to @f$ab@f$ and @f$b^2@f$ dimension.
 */
void expModSphBessel(const Real a2, const Vec1Real& ab, const Vec1Real& b2, Vec2Real& p);

/** wrapper compute modified spherical Bessel function for max order = 0
 * @param Int  n   maximal order of the bessel functions to compute
 * @param Int  n   number of grid points
 * @param Real r   array which contains the grid and is of size n
 * @param Real ri  scalar which stores an offset
 * @param Real rj  arrray containing a second grid, (Jonathan thinks this is redundant)
 * @param Real ril contains Bessel function on return
 */
void modifiedSphericalBesselInterface(const Int       n,
                                      const Vec1Real& r,
                                      const Real      ri,
                                      const Vec1Real& rj,
                                      Vec1Real&       ril);

/** wrapper compute modified spherical Bessel function for max order > 0
 * @param Int lmax maximal order of the bessel functions to compute
 * @param Int  n   number of grid points
 * @param Real r   array which contains the grid and is of size n
 * @param Real ri  scalar which stores an offset
 * @param Real rj  arrray containing a second grid, (Jonathan thinks this is redundant)
 * @param Real ril contains Bessel function on return
 */
void modifiedSphericalBesselInterface(const Int       lmax,
                                      const Int       n,
                                      const Vec1Real& r,
                                      const Real      ri,
                                      const Vec1Real& rj,
                                      Vec2Real&       ril);

/** compute spherical Bessel function and derivative for max order = 0
 * should be used from modifiedSphericalBesselInterface in math_c
 * @param Real derivative contains derivative of the function on output
 * @param Int n           number of grid points
 * @param Real r          array which contains the grid and is of size n
 * @param Real function   contains Bessel function on return
 */
void modifiedSphericalBesselDerivative(Vec1Real&       derivative,
                                       const Int       n,
                                       const Vec1Real& r,
                                       Vec1Real&       function);

/** compute integral over std::vector with trapez, and uniform
 * @param func function to integrate
 * @param dx   grid spacing
 */
template<typename T>
T trapezIntegral(const Vec1Real& func, const Real dx)
{
    Real integral = (Real)0.5 * func[0];
    integral += (Real)0.5 * func[func.size() - 1];
    integral += std::accumulate(func.begin() + 1, func.end() - 1, 0);

    return integral * dx;
}

/** compute integral of squared function on the spherical radial grid x
 * @param func function to integrate
 * @param x    radial grid hast to be uniform
 * @param dx   grid spacing
 */
template<typename T>
T trapezSphericalL1Integral(const std::vector<T>& func, const std::vector<T>& x, const Real dx)
{
    assert((func.size() == x.size())
           && "trapezSphericalL1Integral: dimension of function and xaxis do not match");
    Real tmp = x[0] * x[0];
    Real integral = (Real)0.5 * func[0] * tmp;
    tmp = x[func.size() - 1] * x[func.size() - 1];
    integral += (Real)0.5 * func[func.size() - 1] * tmp;
    for (auto i = 1; i < func.size() - 1; i++)
    {
        tmp = x[i] * x[i];
        integral += func[i] * tmp;
    }

    return integral * dx;
}

/** compute integral of squared function on the spherical radial grid x
 * @param func function to integrate
 * @param x    radial grid hast to be uniform
 * @param dx   grid spacing
 */
template<typename T>
T trapezSphericalL2Integral(const std::vector<T>& func, const std::vector<T>& x, const Real dx)
{
    assert((func.size() == x.size())
           && "trapezSphericalL2Integral: dimension of function and xaxis do not match");
    Real tmp = func[0] * x[0];
    Real integral = (Real)0.5 * tmp * tmp;
    tmp = func[func.size() - 1] * x[func.size() - 1];
    integral += (Real)0.5 * tmp * tmp;
    for (Size i = 1; i < func.size() - 1; i++)
    {
        tmp = func[i] * x[i];
        integral += tmp * tmp;
    }

    return integral * dx;
}
/** multiply vector times scalar, note a copy of the vector is made
 *
 * @param Real vec    input vector to multiply
 * @param Real scalar scalar to multiply vec with
 * @param returned will be the multiplied vector
 *
 */
template<typename T>
std::vector<T> vectorTimesScalar(const std::vector<T>& vec, const T scalar)
{
    std::vector<T> result = vec;

    std::transform(result.begin(),
                   result.end(),
                   result.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, scalar));

    return result;
}
/** multiply vector times scalar, no copy of the input vector will be made
 *@param vec    -> input vector to multiply element wise by scalar, vector will be altered
 *@param scalar -> sclar to multiply vector with
 */
template<typename T>
void vectorTimesScalarNoCopy(std::vector<T>& vec, const T scalar)
{
    std::transform(vec.begin(),
                   vec.end(),
                   vec.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, scalar));
}

/*******************************************************************************************
 * multiply 2D vector by a scalar without making a copy of the data.
 *
 * @param vec vector which will be rescaled in the routine. Note the original vector will
 * be destroyed.
 * @param scalar by which the 2D vector will be rescaled
 *******************************************************************************************/
template<typename T>
void vectorTimesScalarNoCopy(std::vector<std::vector<T>>& vec, const T scalar)
{
    std::for_each(vec.begin(),
                  vec.end(),
                  [scalar](std::vector<T>& innerVector)
                  {
                      std::transform(innerVector.begin(),
                                     innerVector.end(),
                                     innerVector.begin(),
                                     [scalar](T value) { return value * scalar; });
                  });
}
/** square vector elementwise
 * @param Real vec input vector to square
 * @param output is a Real vector containing the input vector squared
 */
template<typename T>
std::vector<T> squareVector(const std::vector<T>& vec)
{
    std::vector<T> result = vec;

    std::transform(result.begin(), result.end(), result.begin(), [](T x) { return x * x; });

    return result;
}
/** compute dot product between two vectors
 *@param Real vecA input vector has to have same length as vecB
 *@param Real vecB input vector to compute dot product with
 *@param return is the dot product between the two vectors
 */
template<typename T>
T dotProduct(const std::vector<T>& vecA, const std::vector<T>& vecB)
{
    assert((vecA.size() == vecB.size())
           && "dotProduct: dimension of function and xaxis do not match");

    return std::inner_product(vecA.begin(), vecA.end(), vecB.begin(), (T)0);
}
/**************************************************************************************************
 * Compute dot product between two vectors (only first n entries).
 *
 * @param vecA Input vector A (must have at least n elements).
 * @param vecB Input vector B (must have at least n elements).
 * @param return Dot product between the two vectors.
 **************************************************************************************************/
template<typename T>
T dotProduct(const std::vector<T>& vecA, const std::vector<T>& vecB, Size n)
{
    VASPML_DEBUG_L1(
        if (vecA.size() < n) global::tutor.bug("Vector A in dot product is too short");
        if (vecB.size() < n) global::tutor.bug("Vector B in dot product is too short");
    );

    return std::inner_product(vecA.begin(), vecA.begin() + n, vecB.begin(), (T)0);
}

template<typename T>
T dotProduct(const T*& vecA, const T*& vecB, const Size& n)
{
    T result = (T)0;
    for (Size i = 0; i < n; i++) result += vecA[i] * vecB[i];

    return result;
}

template<typename T>
T l2Norm(const std::vector<T>& vecA, const Size& n)
{
    T result = (T)0;
    for (Size i = 0; i < n; i++) result += vecA[i] * vecA[i];

    return std::sqrt(result);
}

/** do simpson integration on a integer grid
 * @param data array containing the the function to integrate
 */
template<typename T>
T simpsonIntegration(const std::vector<T>& data)
{

    T result = 0;
    T temp1 = 0;
    T temp2 = 0;

    for (auto i = 2; i < data.size(); i += 2)
    {
        temp1 += data[i];
        temp2 += data[i - 1];
    }
    temp1 *= 4.0;
    temp2 *= 2.0;
    result = temp1 + temp2 + data[data.size() - 1] + data[0];
    result /= 3.0;

    return result;
}

/** do simpson integration on a integer grid, 2-dimensional
 * @param data  2d vector containing the function values
 * @param returns the 2dimensional integral on integer grid
 */
template<typename T>
T simpsonIntegration(const std::vector<std::vector<T>>& data)
{
    std::vector<T> integral;
    integral.resize(data.size());
    Real result;

    for (auto i = 0; i < data.size(); i++) { integral[i] = D1_simps(data[i]); }
    result = D1_simps(integral);

    return result;
}

/** compute square of number
 * @param x input number to square
 */
template<typename T>
static inline Real computeSquare(T x)
{
    return x * x;
}

/** compute square of std::vector
 * @param vec_in input vector to square elementqise
 *        on output it contains the squared entries
 */
template<typename T>
void squareVector(std::vector<T>& vec_in)
{
    std::transform(vec_in.begin(), vec_in.end(), vec_in.begin(), computeSquare<T>);
}

/**
 * sum up all elements of the supplied vector
 *@param input vector of numeric type
 */
template<typename T>
T sumVector(const std::vector<T>& input)
{
    return std::accumulate(input.begin(), input.end(), (T)0);
}

/*******************************************************************************************
 * @brief Computes the partial sums of the elements in the input vector.
 *
 * @tparam T The type of the elements in the input vector. Must be a numeric type.
 * @param input A constant reference to the input vector of type `T`.
 * @return A vector of type `T` containing the partial sums of the input vector.
 *
 * @note If the input vector is empty, the function returns an empty vector.
 * @note The function includes a debug check to ensure that the template parameter `T`
 *       is a numeric type. If `T` is non-numeric, a debug error is triggered.
 *
 * @warning The debug check relies on the `VASPML_DEBUG_L1` macro and the `global::tutor.bug`
 *          function. Ensure these are properly defined in your codebase.
 *******************************************************************************************/
template<typename T>
std::vector<T> partialSum(const std::vector<T>& input)
{
    VASPML_DEBUG_L1(
        if (!std::is_arithmetic<T>())
        {
            global::tutor.bug("Error:" + flf(VASPML_FLF)
                              + " supplied template parameter is NON numeric!");
        }
    );
    if (input.empty()) { return {}; }
    std::vector<T> partialSums(input.size());
    std::partial_sum(input.cbegin(), input.cend(), partialSums.begin());

    return partialSums;
}

/*******************************************************************************************
 * @brief Computes the partial sums of the elements in the input vector starting at zero.
 *
 * @tparam T The type of the elements in the input vector. Must be a numeric type.
 * @param input A constant reference to the input vector of type `T`.
 * @return A vector of type `T` containing the partial sums of the input vector.
 *
 * @note If the input vector is empty, the function returns an empty vector.
 * @note The function includes a debug check to ensure that the template parameter `T`
 *       is a numeric type. If `T` is non-numeric, a debug error is triggered.
 *
 * @warning The debug check relies on the `VASPML_DEBUG_L1` macro and the `global::tutor.bug`
 *          function. Ensure these are properly defined in your codebase.
 *******************************************************************************************/
template<typename T>
std::vector<T> partialSum0(const std::vector<T>& input)
{
    VASPML_DEBUG_L1(
        if (!std::is_arithmetic<T>())
        {
            global::tutor.bug("Error:" + flf(VASPML_FLF)
                              + " supplied template parameter is NON numeric!");
        }
    );
    if (input.empty()) { return {}; }
    std::vector<T> partialSums(input.size());
    partialSums[0] = (T)0;
    std::partial_sum(input.cbegin(), input.cend() - 1, partialSums.begin() + 1);

    return partialSums;
}

/*******************************************************************************************
 * Multiplies two vectors elementwise. Computes Hadarmad product
 *
 * @param vectorA first input vector which has to match the size of vectorB
 * @param vectorB second vectr which is multiplied with vectorA
 *
 * @f[
 * r_{i} = a_{i} * b_{i}
 * @f]
 *******************************************************************************************/
template<typename T>
std::vector<T> elementwiseProduct(const std::vector<T>& vectorA, const std::vector<T>& vectorB)
{
    VASPML_DEBUG_L1(
        if (vectorA.size() != vectorB.size())
        {
            global::tutor.bug("ERROR: math::elementwiseProduct( const std::vector<T>& "
                              "vectorA, const std::vector<T>& vectorB ) \n"
                              "length of input vectors does not match");
        }
    );
    std::vector<T> result(vectorA.size());
    std::transform(vectorA.begin(),
                   vectorA.end(),
                   vectorB.begin(),
                   result.begin(),
                   std::multiplies<T>());

    return result;
}

/*******************************************************************************************
 * computes an integer power of some supplied value by explicit for loop
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 *
 * For small integer exponents an explicit loop to compute the nth power of a function
 * is often faster than the c++ std::pow function. For large integer powers the function
 * intPowSquareAlgo should be prefered.
 *******************************************************************************************/
template<typename T>
inline T intPowLoop(const T& base, const Int& exp)
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (exp < 0)
        {
            global::tutor.bug("ERROR: inline T intPowLoop( const T& base, const Int& exp )\n"
                              "exp < 0. Function can only compute powers exp >= 0");
        }
    );
#endif

    T result = 1.0;
    for (Int i = 0; i < exp; i++) { result *= base; }

    return result;
}

/*******************************************************************************************
 * computes an integer power of some supplied value by use of explicit cases
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 *
 * For small integer exponents cases with hard coding the power can be faster than
 * the c++ std::pow function. For large integer powers the function
 * intPowSquareAlgo is called. This function was in tests a little slower than the intPowLoop
 * but still faster than std::pow
 *******************************************************************************************/
template<typename T>
inline T intPowCase(const T& base, const Int& exp)
{
    VASPML_DEBUG_L1(
        if (exp < 0)
        {
            global::tutor.bug("inline T intPowCase( const T& base, const Int& exp ) \n"
                              "exp < 0. Function can only compute powers exp >= 0");
        }
    );
    switch (exp)
    {
    case 0:
        return 1;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    case 6:
        return base * base * base * base * base * base;
    case 7:
        return base * base * base * base * base * base * base;
    case 8:
        return base * base * base * base * base * base * base * base;
    case 9:
        return base * base * base * base * base * base * base * base * base;
    case 10:
        return base * base * base * base * base * base * base * base * base * base;
    default:
        return intPowLoop(base, exp);
    }
}

/*******************************************************************************************
 * computes an integer power of some supplied value by use of a square algorithm
 *
 * @param baseIn value of which the power will be computed
 * @param expIn determines the order of the power which is computed of base
 *
 * For large integer exponents the square algorithm will be faster than intPowCase, intPowLoop
 * and the std::pow algorithm. But for small integer powers the intPowLoop or intPowCase
 * should be preferred. Testing with gnu compiler showed around exp=16 this function is
 * fastest.
 *******************************************************************************************/
template<typename T>
inline T intPowSquareAlgo(const T& baseIn, const Int& expIn)
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (expIn < 0)
        {
            global::tutor.bug("inline T intPowSquareAlgo( const T& baseIn, const Int& expIn )\n"
                              "exp < 0. Function can only compute powers exp >= 0");
        }
    );
#endif
    T   result = 1.0;
    Int exp = expIn;
    T   base = baseIn;
    while (exp)
    {
        if (exp % 2 == 1) { result *= base; }
        base *= base;
        exp /= 2;
    }

    return result;
}

/*******************************************************************************************
 * computes integer power of some supplied value by use of a square algorithm or explicit loop
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 * For small integer exponents the explicit loop algorithm is used and for larger exponents
 * the intPowSquareAlgo is used.
 *******************************************************************************************/
template<typename T>
inline T intPowMixAlgo(const T& baseIn, const Int& expIn)
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (expIn < 0)
        {
            global::tutor.bug("inline T intPowMixAlgo( const T& baseIn, const Int& expIn )\n"
                              "exp < 0. Function can only compute powers exp >= 0");
        }
    );
#endif
    if (expIn < 10) return intPowLoop(baseIn, expIn);
    else return intPowSquareAlgo(baseIn, expIn);
}

/*******************************************************************************************
 * compute average value of input vector
 *
 * @param dataIn vector of which elements the average is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T average(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;

    return sumVector(dataIn) / (T)dataIn.size();
}

/*******************************************************************************************
 * compute minimum value of input vector
 *
 * @param dataIn vector of which elements the minimum is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T minimum(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;

    return *std::min_element(dataIn.cbegin(), dataIn.cend());
}

/*******************************************************************************************
 * compute maximum value of input vector
 *
 * @param dataIn vector of which elements the maximum is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T maximum(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;

    return *std::max_element(dataIn.cbegin(), dataIn.cend());
}

/*******************************************************************************************
 * compute prime numbers up to a given number n using Eratosthenes sieve
 *
 * the sieve of Eratosthenes is an ancient algorithm for
 * finding all prime numbers up to any given limit. Routine was taken from vasp and was
 * originally written by Merzuk Kaltak
 *
 * @param[in] n upper bound to which prime numbers are computed
 * @param[out] primSqrt stores square roots of prime numbers which when squared are smaller n
 * @param[out] stores prime numbers up to n
 *******************************************************************************************/
template<typename T>
void computePrimeNumbers(const T n, std::vector<T>& primSqrt, std::vector<T>& primes)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("computePrimeFactor( const T n, std::vector<T>& "
                              "primSqrt, std::vector<T>& primes )\n"
                              "template data type is not an integral type");
        }
    );

    T                 maxPrime = std::floor(std::sqrt((Real)n));
    std::vector<bool> a(n, true);
    for (T i = 2; i <= maxPrime; i++)
    {
        T j = i * i;
        while (true)
        {
            a[j - 1] = false;
            j = j + i;
            if (j > n) break;
        }
    }

    T size = std::count(a.begin(), a.end(), true);
    primes.resize(size - 1);
    T j = 1;
    // skipping 1
    for (T i = 2; i <= n; i++)
    {
        if (!a[i - 1]) continue;
        primes[j - 1] = i;
        j++;
        if (j > (T)primes.size()) break;
    }

    j = 0;
    for (T i = 0; i < (Int)primes.size(); i++)
    {
        if (primes[i] * primes[i] > n)
        {
            j++;
            break;
        }
        j++;
    }
    primSqrt.resize(j);
    for (T i = 0; i < j; i++) primSqrt[i] = primes[i];
}

/*******************************************************************************************
 * compute prime factorization of given input number
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param param[in] number for which the prime factorization is obtained
 * @param param[in]  primes prime numbers smaller than n, can be obtained by
 * computePrimeNumbers( const T n, std::vector<T>& primSqrt, std::vector<T>& primes )
 * @param[ out ] nexp prime factors of n
 *******************************************************************************************/
template<typename T>
void primeFactorization(const T n, const std::vector<T>& primes, std::vector<T>& nexp)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("void primeFactorization( const T n, std::vector<T>& nexp, "
                              "const std::vector<T>& primes )\n"
                              "template data type is not an integral type");
        }
    );
    T nprimes = primes.size();
    nexp.resize(primes.size(), (T)0);
    // 1 is not considered as prime
    T m = n;

    for (T ip = 0; ip < nprimes; ip++)
    {
        while (true)
        {
            if (m % primes[ip] == 0)
            {
                m = m / primes[ip];
                nexp[ip]++;
            }
            else { break; }
        }
    }
}
/*******************************************************************************************
 * returns true if supplied value is perfect square ie if it can be expressed as n*n=value
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[ in ] value value to be checked if it is expressible as n * n
 *******************************************************************************************/
template<typename T>
bool isPerfectSquare(const T value)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("bool isPerfectSquare( const T value )\n"
                              "template data type is not an integral type");
        }
    );
    bool result;

    T ilimit = std::ceil(std::sqrt((Real)value));

    T i;
    for (i = 1; i <= ilimit; i++)
    {
        if (i * i == value) break;
    }

    if (i > ilimit) result = false;
    else result = true;

    return result;
}
/*******************************************************************************************
 * if the supplied value can be expressed as perfect square this function returns n
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] value which will be expressed as the return value n * n
 *******************************************************************************************/
template<typename T>
T perfectSquare(const T value)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("T perfectSquare( const T value )\n"
                              "template data type is not an integral type");
        }
    );

    T ilimit = std::ceil(std::sqrt((Real)value));
    T result;
    for (result = 1; result <= ilimit; result++)
    {
        if (result * result == value) break;
    }

    return result;
}
/*******************************************************************************************
 * check if a supplied number is a prime number; returns true if number is prime
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] value value which will be check if it is a prime number
 * @param[in] primes prime numbers wich are smaller equal the supplied value; can be obtained
 * with void computePrimeNumbers( const T n, std::vector<T>& primSqrt, std::vector<T>& primes )
 *******************************************************************************************/
template<typename T>
bool isPrime(const T value, const std::vector<T>& primes)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("bool isPrime( const T value, const std::vector<T>& primes )\n"
                              "template data type is not an integral type");
        }
    );
    bool result;
    T    ilimit = (T)primes.size();
    T    i = 0;
    for (T i = 0; i < ilimit; i++)
    {
        if (value == primes[i]) break;
    }
    if (i >= ilimit) result = false;
    else result = true;

    return result;
}

/*******************************************************************************************
 * Factorize an integer number into a product of two numbers.
 *
 * The factorization is done under the constraint that the two numbers are as close as
 * possible to each other. This routine is very helpful when decomposing a parallel
 * environment of MPI into a scalapack grid. The number of processors will be decomposed
 * in a 2 dimensional process grid where the number of rows is as close as possible to
 * the number of columns.
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] nprocs is the number which will be decomposed into a product of nprow * npcol
 * @param[out] nprow is the first number
 *
 *******************************************************************************************/
template<typename T>
void fermatRazor(const T nprocs, T& nprow, T& npcol)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global::tutor.bug("void fermatRazor( const T nprocs, T nprow, T npcol )\n"
                              "template data type is not an integral type");
        }
    );
    if (nprocs == 1)
    {
        nprow = 1;
        npcol = 1;
        return;
    }

    std::vector<T> primes;
    std::vector<T> allPrimes;
    // obtain all primes up to sqrt( nprocs  )
    computePrimeNumbers(nprocs, primes, allPrimes);
    std::vector<T> nexp;
    // obtain prime factorization of NPROCS
    primeFactorization(nprocs, primes, nexp);

    // first check if NPROCS is a perfect square
    if (isPerfectSquare(nprocs))
    {
        nprow = perfectSquare(nprocs);
        npcol = nprow;
    }
    else
    {
        // the idea is to factorize nprocs = a*b
        // with restriction |A - B| -> 0
        // fermat sieve works well if N is odd
        T nodd = nprocs / (std::pow(primes[0], nexp[0]));
        // in case nprocs is even
        if (nodd == 1)
        {
            nprow = std::pow(primes[0], nexp[0] / 2);
            npcol = std::pow(primes[0], nexp[0] / 2 + nexp[0] % 2);
        }
        // in case nprocs is not even
        // use fermats sieve on nodd
        else
        {
            // perfect square nothing to do
            if (isPerfectSquare(nodd))
            {
                nprow = perfectSquare(nodd);
                npcol = perfectSquare(nodd);
            }
            else
            {
                if (isPrime(nodd, allPrimes))
                {
                    nprow = 1;
                    npcol = nodd;
                }
                else
                {
                    // Fermat's way to factorize an odd integer:
                    // nodd = a**2- b**2 = ( a - b ) * ( a + b )
                    // and write a**2 - n = b**2
                    // and increase A by 1 until A**2 - N is a perfect square
                    T a = std::ceil(std::sqrt((Real)nodd));
                    T b = a * a - nodd;
                    while (true)
                    {
                        if (isPerfectSquare(b)) break;
                        a++;
                        b = a * a - nodd;
                    }
                    nprow = a + perfectSquare(b);
                    npcol = a - perfectSquare(b);
                } // nodd is prime
                for (T i = 0; i < nexp[0]; i++)
                {
                    if (nprow < npcol) nprow = nprow * 2;
                    else npcol = npcol * 2;
                }
            } // perfect square nodd
        } // nprocs not even
    } // perfec square

    if (nprow * npcol != nprocs)
    {
        global::tutor.error("ERROR: void fermatRazor( const T nprocs, T& nprow, T& npcol )\n"
                            "failed "
                            + std::to_string(nprocs) + "  " + std::to_string(nprow)
                            + std::to_string(npcol) + "\n");
    }
}

namespace stats
{

/*******************************************************************************************
 * @brief Calculates the mean (average) of a vector of numeric values.
 *
 * This function computes the arithmetic mean of the elements in the input vector.
 * It requires the vector to contain numeric types (e.g., integers, floating-point numbers).
 * If the vector is empty, a warning is logged, and the behavior is undefined.
 *
 * @tparam T The type of elements in the input vector. Must be a numeric type.
 * @param data A vector containing the data for which the mean is to be calculated.
 * @return The mean of the elements in the vector as a `Real` type.
 *
 * @note If the input vector is empty, a warning is logged, and the function may return an undefined
 *value.
 * @note The function uses `std::accumulate` to compute the sum of the elements.
 * @warning If the input vector contains non-numeric types, a debug message is logged, and the
 *behavior is undefined.
 *******************************************************************************************/
template<typename T>
Real calculateMean(const std::vector<T>& data)
{
    VASPML_DEBUG_L1(
        if (!std::is_arithmetic<T>::value)
        {
            global::tutor.bug("calculateMean: Input vector must contain numeric types.");
        }
        if (data.empty())
        {
            global::tutor.error("Warning: Attempted to calculate mean of an empty vector.");
        }
    );
    Real sum = (Real)std::accumulate(data.begin(), data.end(), 0.0);

    return sum / (Real)data.size();
}

/*******************************************************************************************
 * @brief Calculates the sample variance of a vector of numeric values.
 *
 * This function computes the sample variance of the elements in the input vector.
 * The sample variance is calculated using the formula:
 * \f[
 * \text{Variance} = \frac{\sum_{i=1}^{N} (x_i - \text{mean})^2}{N-1}
 * \f]
 * If the mean is not provided, it is calculated internally using `calculateMean`.
 * If the vector contains only one element, the variance is defined as 0.
 *
 * @tparam T The type of elements in the input vector. Must be a numeric type.
 * @param data A vector containing the data for which the variance is to be calculated.
 * @param meanVal (Optional) The mean of the data. If not provided, it is calculated internally.
 * @return The sample variance of the elements in the vector as a `Real` type.
 *
 * @note If the input vector is empty, an error is logged, and the behavior is undefined.
 * @note The function uses `std::accumulate` with a lambda function to compute the sum of squared
 *differences.
 * @warning If the input vector contains non-numeric types, the behavior is undefined.
 * @warning If the vector size is 1, the variance is defined as 0.
 *******************************************************************************************/
template<typename T>
Real calculateVariance(const std::vector<T>& data,
                       Real                  meanVal = std::numeric_limits<Real>::max())
{

    VASPML_DEBUG_L1(
        if (!std::is_arithmetic<T>::value)
        {
            global::tutor.bug("calculateVariance: Input vector must contain numeric types.");
        }
        if (data.empty()) { global::tutor.error("calculateVariance: Input vector is empty."); }
    );
    if (data.size() == 1)
    {
        // Sample variance for a single data point is 0
        return 0.0;
    }
    if (meanVal == std::numeric_limits<Real>::max())
    {
        // Mean was NOT supplied, so calculate it internally (will return double)
        meanVal = calculateMean(data);
    }
    Real sumOfSquaredDifferences =
        std::accumulate(data.cbegin(),
                        data.cend(),
                        0.0,
                        [&meanVal](Real& accumulator, const T& current_value)
                        {
                            // Calculate the difference, promoting current_value to double
                            Real diff = static_cast<Real>(current_value) - meanVal;
                            // Square the difference using multiplication (faster than std::pow for
                            // exponent 2)
                            return accumulator + (diff * diff);
                        });

    // Using (N-1) for sample variance
    return sumOfSquaredDifferences / (Real)(data.size());
}
/*******************************************************************************************
 * @brief Computes the Root Mean Square Error (RMSE) between two vectors.
 *
 * This function calculates the RMSE, a measure of the differences between
 * corresponding elements of two input vectors. The RMSE is commonly used
 * to evaluate the accuracy of predictions in regression tasks.
 *
 * @tparam T The type of elements in the input vectors. Must be a numeric type.
 * @param setA The first input vector.
 * @param setB The second input vector. Must have the same size as `setA`.
 * @return The RMSE value as a `Real` type.
 *
 * @note The function assumes that the input vectors contain numeric types.
 *       If the sizes of the vectors differ or the type is non-numeric,
 *       debug-level checks will trigger error messages.
 *
 * @throws std::runtime_error If debug checks are enabled and the input
 *         vectors do not meet the requirements.
 *
 * @warning Ensure that the input vectors are of the same size and contain
 *          numeric types to avoid undefined behavior.
 *******************************************************************************************/
template<typename T>
Real calculateRMSE(const std::vector<T>& setA, const std::vector<T>& setB)
{
    VASPML_DEBUG_L1(
        if (setA.size() != setB.size())
        {
            global::tutor.bug("calculateRMSE: Input vectors must have the same size.");
        }
        if (!isNumberType<T>())
        {
            global::tutor.bug("calculateRMSE: Input vectors must have be numeric types.");
        }
    );

    Real sum = 0.0;
    for (Size i = 0; i < setA.size(); ++i)
    {
        Real diff = static_cast<Real>(setA[i]) - static_cast<Real>(setB[i]);
        sum += diff * diff;
    }

    return std::sqrt(sum / static_cast<Real>(setA.size()));
}
/*******************************************************************************************
 * @brief Computes the Mean Absolute Error (MAE) between two vectors.
 *
 * This function calculates the Mean Absolute Error (MAE) between two input vectors `setA` and
 *`setB`. The MAE is defined as the average of the absolute differences between corresponding
 *elements in the two vectors.
 *
 * @tparam T The type of elements in the input vectors. Must be a numeric type.
 * @param setA The first input vector.
 * @param setB The second input vector. Must have the same size as `setA`.
 * @return The Mean Absolute Error (MAE) as a `Real` value.
 *
 * @note The function performs runtime checks to ensure that:
 *       - The input vectors have the same size.
 *       - The type `T` is numeric.
 *       If these conditions are not met, debug messages are triggered.
 *
 * @warning This function assumes that the type `Real` is defined and represents a floating-point
 *type. Ensure that `isNumberType<T>()` is properly implemented to validate numeric types.
 *******************************************************************************************/
template<typename T>
Real calculateMAE(const std::vector<T>& setA, const std::vector<T>& setB)
{
    VASPML_DEBUG_L1(
        if (setA.size() != setB.size())
        {
            global::tutor.bug("calculateMAE: Input vectors must have the same size.");
        }
        if (!isNumberType<T>())
        {
            global::tutor.bug("calculateMAE: Input vectors must have be numeric types.");
        }
    );

    Real sum = 0.0;
    for (Size i = 0; i < setA.size(); ++i)
    {
        Real diff = static_cast<Real>(setA[i]) - static_cast<Real>(setB[i]);
        sum += std::abs(diff);
    }

    return sum / static_cast<Real>(setA.size());
}

/*******************************************************************************************
 * @brief Computes the Pearson correlation coefficient between two vectors.
 *
 * This function calculates the Pearson correlation coefficient, which measures
 * the linear correlation between two datasets. The result is a value between
 * -1 and 1, where 1 indicates a perfect positive linear relationship, -1 indicates
 * a perfect negative linear relationship, and 0 indicates no linear relationship.
 *
 * @tparam T The type of the elements in the input vectors. Must be a numeric type.
 * @param x The first input vector.
 * @param y The second input vector. Must have the same size as `x`.
 * @return The Pearson correlation coefficient as a `Real` value.
 *
 * @throws std::runtime_error If the input vectors are not numeric, have different sizes,
 *                            are empty, or if a division by zero occurs during computation.
 *
 * @note This function uses debug checks (VASPML_DEBUG_L1) to validate input and
 *       detect potential issues during computation.
 *
 * @warning The input vectors must not be empty, must have the same size, and must
 *          contain numeric types. Division by zero may occur if the input data
 *          has no variance.
 *******************************************************************************************/
template<typename T>
Real pearsonCorrPtr(const T* x, const T* y, const Size n)
{
    VASPML_DEBUG_L1(
        if (!isNumberType<T>())
        {
            global::tutor.bug("ERROR pearsonr: Input vectors must have be numeric types.");
        }
    );
    // Calculate means
    Real mean_x = (Real)std::accumulate(x, x + n, 0.0) / (Real)n;
    Real mean_y = (Real)std::accumulate(y, y + n, 0.0) / (Real)n;

    // Calculate numerator and denominators using std algorithms
    Real numerator = std::inner_product(x,
                                        x + n,
                                        y,
                                        0.0,
                                        std::plus<>(),
                                        [&mean_x, &mean_y](const T& xi, const T& yi)
                                        { return ((Real)xi - mean_x) * ((Real)yi - mean_y); });

    Real sum_sq_x = std::accumulate(x,
                                    x + n,
                                    0.0,
                                    [&mean_x](Real& acc, const T& xi)
                                    { return acc + ((Real)xi - mean_x) * ((Real)xi - mean_x); });

    Real sum_sq_y = std::accumulate(y,
                                    y + n,
                                    0.0,
                                    [&mean_y](Real acc, const T& yi)
                                    { return acc + ((Real)yi - mean_y) * ((Real)yi - mean_y); });

    Real denominator = std::sqrt(sum_sq_x * sum_sq_y);

    VASPML_DEBUG_L1(
        if (denominator == 0.0)
            global::tutor.error("ERROR pearsonr: Division by zero in Pearson calculation");
    );

    return numerator / denominator;
}

template<typename T>
Real pearsonCorr(const std::vector<T>& x, const std::vector<T>& y)
{
    VASPML_DEBUG_L1(
        if (!isNumberType<T>())
        {
            global::tutor.bug("ERROR pearsonr: Input vectors must have be numeric types.");
        }
        if (x.size() != y.size())
            global::tutor.bug("ERROR pearsonr: Vectors must be the same length\n");
        if (x.empty()) global::tutor.bug("ERROR pearsonr: Vectors must not be empty");
    );
    const Size n = x.size();
    //// Calculate means
    //Real mean_x = (Real)std::accumulate( x.cbegin(), x.cend(), 0.0 ) / (Real)n;
    //Real mean_y = (Real)std::accumulate( y.cbegin(), y.cend(), 0.0 ) / (Real)n;

    //// Calculate numerator and denominators using std algorithms
    //Real numerator = std::inner_product(x.cbegin(), x.cend(), y.cbegin(), 0.0,
    //    std::plus<>(),
    //    [&mean_x, &mean_y]( const T& xi, const T& yi) {
    //        return ((Real)xi - mean_x) * ((Real)yi - mean_y);
    //    });

    //Real sum_sq_x = std::accumulate(x.cbegin(), x.cend(), 0.0, [&mean_x](Real& acc, const T& xi) {
    //    return acc + ((Real)xi - mean_x) * ((Real)xi - mean_x);
    //});

    //Real sum_sq_y = std::accumulate(y.cbegin(), y.cend(), 0.0, [&mean_y](Real acc, const T& yi) {
    //    return acc + ((Real)yi - mean_y) * ((Real)yi - mean_y);
    //});

    //Real denominator = std::sqrt(sum_sq_x * sum_sq_y);

    //VASPML_DEBUG_L1(
    //    if (denominator == 0.0)
    //        global::tutor.error("ERROR pearsonr: Division by zero in Pearson calculation");
    //);
    //return numerator / denominator;

    return pearsonCorrPtr(x.data(), y.data(), n);
}

/*******************************************************************************************
 * @brief Computes the correlation matrix for a set of vectors.
 *
 * This function calculates the Pearson correlation coefficient between
 * each pair of vectors in the input and stores the results in a symmetric
 * correlation matrix.
 *
 * @tparam T The data type of the elements in the vectors (e.g., float, double).
 * @param vectors A 2D vector where each inner vector represents a data series.
 * @return A 2D vector representing the symmetric correlation matrix, where
 *         the element at (i, j) is the Pearson correlation coefficient
 *         between vectors[i] and vectors[j].
 *
 * @note The function assumes that the input vectors are of equal length
 *       and that the Pearson correlation coefficient is computed using
 *       `math::stats::pearsonCorr`.
 *******************************************************************************************/
template<typename T>
std::vector<T> computeCorrelationMatrix(const std::vector<std::vector<T>>& vectors)
{
    const Int                   numVectors = vectors.size();
    std::vector<std::vector<T>> corrMatrix(numVectors, std::vector<T>(numVectors));
    for (Int i = 0; i < numVectors; i++)
    {
        corrMatrix[i][i] = 1.0;
        for (Int j = i + 1; j < numVectors; j++)
        {
            T corr = math::stats::pearsonCorr(vectors[i], vectors[j]);
            corrMatrix[i][j] = corr;
            corrMatrix[j][i] = corr;
        }
    }
}
/*******************************************************************************************
 * @brief Computes a squared correlation matrix for a set of input vectors.
 *
 * This function calculates the squared Pearson correlation coefficients
 * between all pairs of input vectors and stores the results in a symmetric
 * correlation matrix. The diagonal elements of the matrix are set to 1.0.
 *
 * @tparam T The data type of the elements in the input vectors.
 * @param vectors A 2D vector containing the input vectors for which the
 *                correlation matrix is to be computed.
 * @return A 2D vector representing the squared correlation matrix. The
 *         matrix is symmetric, and the diagonal elements are 1.0.
 *
 * @note This function assumes that the `math::stats::pearsonCorr` function
 *       is available and computes the Pearson correlation coefficient
 *       between two vectors.
 *******************************************************************************************/
template<typename T>
std::vector<T> computeCorrelationMatrixSquare(const std::vector<std::vector<T>>& vectors)
{
    const Int                   numVectors = vectors.size();
    std::vector<std::vector<T>> corrMatrix(numVectors, std::vector<T>(numVectors));
    for (Int i = 0; i < numVectors; i++)
    {
        corrMatrix[i][i] = 1.0;
        for (Int j = i + 1; j < numVectors; j++)
        {
            T corr = math::stats::pearsonCorr(vectors[i], vectors[j]);
            corr *= corr;
            corrMatrix[i][j] = corr;
            corrMatrix[j][i] = corr;
        }
    }
}

} // namespace stats

} // namespace vaspml::math

#endif
