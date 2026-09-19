#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "ParallelEnvironment.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm>   // for min
#include <type_traits> // for is_floating_point_v

namespace vaspml
{

/**************************************************************************************************
 * This class sets up and computes a cubic spline interpolation.
 *
 * The routine is adapted from the cubic spline interpolation algorithm given in
 * ~~~
 * Numerical Recipes in C++: The Art of Scientific Computing 3rd Edition
 * by  William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
 * page 120
 * ISBN-10  0521880688
 * ISBN-13  978-0521880688
 * ~~~
 *
 * Users need to provide an input vector of @f$x@f$-values (`x`) and corresponding function values
 * @f$y(x)@f$ (`y`). Optionally, the first derivatives at the first (@f$y'(x_1)@f$, `dy1`) and/or
 * last (@f$y'(x_n)@f$, `dyn`) point may be added (clamped cubic splines). After setup, the various
 * `interpolate...()` functions may be used to retrieve the interpolated function values @f$f(x)@f$
 * (`fx`) and the corresponding derivatives @f$f'(x)@f$ (`dfx`).
 *
 * Some technical details about modifications to the original code in the literature:
 *
 * The book suggests these steps to evaluate, i.e. interpolate, the original function: assume the
 * x-value lies between @f$[x_{i}, x_{i+1}]@f$, first compute these intermediate quantities:
 * @f[
 * h = x_{i+1} - x_{i}, \quad
 * a = \frac{x_{i+1} - x}{h} \quad \text{and} \quad
 * b = \frac{x - x_{i}}{h}.
 * @f]
 * Then, evaluate the cubic spline with this expression:
 * @f[
 * f(x) = a y_{i} + b y_{i+1} + \frac{h^2}{6} \left[ (a^3-a) y''_{i} + (b^3-b) y''_{i+1} \right]
 * @f]
 * However, this implementation uses a different method: the idea is to evaluate the cubic
 * polynomials according to Horner's scheme (this requires least amount of floating-point
 * operations):
 * @f[
 * f(x) = ((A z + B) z + C) z + D = A z^3 + B z^2 + C z + D,
 * @f]
 * where @f$z := x - x_{i}@f$ and the parameters @f$A, B, C, D@f$ can be precomputed for each node
 * interval @f$i@f$. If we replace @f$a = 1 - b@f$ and regroups terms along powers of @f$z@f$ we
 * get:
 * @f[
 * \begin{split}
 * A_{i} &= \frac{1}{6h}(y''_{i+1} - y''_{i}) \\
 * B_{i} &= \frac{1}{2} y''_{i} \\
 * C_{i} &= \left[ \frac{1}{h}(y_{i+1} - y_{i}) - \frac{h}{3} y''_{i} - \frac{h}{6} y''_{i+1}
 *          \right] \\
 * D_{i} &= y_{i}
 * \end{split}
 * @f]
 * These coefficients are stored together with @f$x_{i}@f$ sequentially for all nodes in the array
 * `c` in this order:
 * @f[
 * \begin{split}
 * \texttt{c[5 * i + 0]} &= x_{i} \\
 * \texttt{c[5 * i + 1]} &= A_{i} \\
 * \texttt{c[5 * i + 2]} &= B_{i} \\
 * \texttt{c[5 * i + 3]} &= C_{i} \\
 * \texttt{c[5 * i + 4]} &= D_{i}
 * \end{split}
 * @f]
 *
 * @warning This implementation assumes an evenly spaced grid of nodes with x-values at
 * @f[
 * x_{i} = (i + 1)\,\Delta, \quad i = 0,\ldots,n-1
 * @f]
 * Hence, the cubic spline interpolation can be evaluated in the interval @f$[\Delta, (n+1)\,
 * \Delta)@f$. The grid limitation is not present in the setup, i.e., any (x,y)-values may be
 * supplied and coefficients will be correctly computed. However, the `interpolate...()` functions
 * do not yet support an irregular grid not starting at @f$\Delta@f$ due to the way the interval
 * index @f$i@f$ is computed.
 *
 * If the provided x-value for interpolation is out-of-range, either a bug message is issued
 * (`VASPML_DEBUG_L2`), or the returned value is huge (no debugging or compilation without
 * exceptions).
 **************************************************************************************************/
template<class T>
class CubicSpline
{
    // Excluding non-floating-point types.
    static_assert(std::is_floating_point_v<T>,
                  "CubicSpline requires a floating-point template type.");

  public:
    /// Default constructor, initializing empty object not ready to use (call setup() first).
    CubicSpline() {};
    /**********************************************************************************************
     * Constructor with full setup.
     *
     * @param x Vector of input x-values of nodes (points of support).
     * @param y Vector of input function values y(x) of nodes.
     * @param dy1 First derivative of function at first node.
     * @param dyn First derivative of function at last node.
     *
     * @note If `dy1` and/or `dyn` are set to the default value or larger natural boundary
     * conditions (second derivative equals zero) are applied to the corresponding end-point.
     * Otherwise, if `dy1` and/or `dyn` are provided and below the default value, clamped cubic
     * splines are used at the respective end. The values are interpreted as the first derivatives
     * of the user-provided function at the end node.
     **********************************************************************************************/
    CubicSpline(const std::vector<T>& x,
                const std::vector<T>& y,
                const T               dy1 = 1.0E99,
                const T               dyn = 1.0E99)
    {
        setup(x, y, dy1, dyn);
    }
    /**********************************************************************************************
     * Set up cubic spline instance given user-function values.
     *
     * Use to refresh existing, already set up instance. Or, call before use if class was
     * instantiated with empty constructor.
     **********************************************************************************************/
    void setup(const std::vector<T>& x,
               const std::vector<T>& y,
               const T               dy1 = 1.0E99,
               const T               dyn = 1.0E99);
    /**********************************************************************************************
     * Compute interpolated function value at given position (single point version).
     *
     * @param x Input x-value for which to compute interpolated function value.
     * @return Function value at given input x-value.
     **********************************************************************************************/
    T interpolate(const T x) const;
    /**********************************************************************************************
     * Compute interpolated function values at given positions (vector version).
     *
     * @param x Vector of input x-values for which to compute interpolated function values.
     * @return Vector of function values at given input x-values.
     **********************************************************************************************/
    std::vector<T> interpolate(const std::vector<T>& x) const;
    /**********************************************************************************************
     * Compute first derivative of splined function (single point version).
     *
     * @param x Input x-value for which to compute derivative of interpolated function.
     * @return Derivative at given input x-value.
     **********************************************************************************************/
    T interpolateDerivative(const T x) const;
    /**********************************************************************************************
     * Compute first derivative of splined function (vector version).
     *
     * @param x Vector of input x-values for which to compute derivatives of interpolated function.
     * @return Vector of derivatives at given input x-values.
     **********************************************************************************************/
    std::vector<T> interpolateDerivative(const std::vector<T>& x) const;
    /**********************************************************************************************
     * Compute value and first derivative of splined function (single point version).
     *
     * @param x Input x-value for which to compute function and derivative of interpolated function.
     * @param fx Returned function value at given input x-value.
     * @param dfx Returned derivative value at given input x-value.
     **********************************************************************************************/
    VASPML_EXEC_SPACE_SPECIFIER
    void interpolateWithDerivative(const T x, T& fx, T& dfx) const;
    /**********************************************************************************************
     * Compute values and first derivatives of splined function (vector version).
     *
     * @param x Vector of input x-values for which to compute functions and derivatives of
     *          interpolated function.
     * @param fx Returned vector containing function values at given input x-values.
     * @param dfx Returned vector containing derivative values at given input x-values.
     **********************************************************************************************/
    VASPML_EXEC_SPACE_SPECIFIER
    void interpolateWithDerivative(const std::vector<T>& x,
                                   std::vector<T>&       fx,
                                   std::vector<T>&       dfx) const;

  private:
    /**********************************************************************************************
     * Compute spline index from given x-value.
     **********************************************************************************************/
    inline Size index(const T x) const;

    /**********************************************************************************************
     * Spline coefficients for evaluation (interpolation) of nested (Horner) polynomial form.
     *
     * @note This vector contains one extra entry at the end: this last bin deliberately contains
     * huge numbers (out-of-range bin) which signals some error even in case of disabled debugging.
     **********************************************************************************************/
    std::vector<T> c;
    /// Grid distance @f$\Delta@f$.
    T delta;
    /// Number of nodes.
    Size n;
};

template<class T>
void CubicSpline<T>::setup(const std::vector<T>& x,
                           const std::vector<T>& y,
                           const T               dy1,
                           const T               dyn)
{
    VASPML_DEBUG_L1(
        if (x.size() != y.size())
        {
            global::tutor.bug(str("Cannot construct cubic spline interpolation: Input "
                                  "vectors have different sizes (x: %zu, y: %zu).",
                                  x.size(),
                                  y.size()));
        }
        if (x.size() < 3)
        {
            global::tutor.bug(str("Cannot construct cubic spline interpolation: Not enough "
                                  "node points (%zu).",
                                  x.size()));
        }
    );

    // Number of nodes.
    n = x.size();

    // Compute grid spacing.
    delta = x[1] - x[0];

    // Resize spline coefficient vector.
    // Reserve one extra bin for out-of-range detection without debugging enabled.
    c.clear();
    c.resize(5 * (n + 1), static_cast<T>(0));
    // The last bin ("out-of-range"-bin) deliberately contains garbage value.
    c[5 * n] = 1.E99;
    c[5 * n + 1] = 1.E99;
    c[5 * n + 2] = 1.E99;
    c[5 * n + 3] = 1.E99;
    c[5 * n + 4] = 1.E99;

    // Temporary variables.
    std::vector<T> y2(n);
    std::vector<T> u(n);
    T              qn;
    T              un;
    T              sig;
    T              p;

    // Boundary condition at first point.
    if (dy1 > 0.99E99)
    {
        y2[0] = static_cast<T>(0);
        u[0] = static_cast<T>(0);
    }
    else
    {
        y2[0] = static_cast<T>(-0.5);
        u[0] = (static_cast<T>(3) / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - dy1);
    }

    // Process spline nodes.
    for (Size i = 1; i < n - 1; i++)
    {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + static_cast<T>(2);
        y2[i] = (sig - static_cast<T>(1)) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (static_cast<T>(6) * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }

    // Boundary condition at last point.
    if (dyn > 0.99E99)
    {
        qn = static_cast<T>(0);
        un = static_cast<T>(0);
    }
    else
    {
        qn = static_cast<T>(0.5);
        un = (static_cast<T>(3) / (x[n - 1] - x[n - 2]))
           * (dyn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }

    // Process spline nodes again.
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + static_cast<T>(1));
    for (Int i = n - 2; i >= 0; --i) { y2[i] = y2[i] * y2[i + 1] + u[i]; }

    // Prepare spline coefficient vector.
    for (Size i = 0; i < n - 1; ++i)
    {
        T* coeff = &(c[5 * i]);
        sig = x[i + 1] - x[i];
        p = (y2[i + 1] - y2[i]) / static_cast<T>(6);
        coeff[0] = x[i];
        coeff[1] = p / sig;
        coeff[2] = y2[i] / static_cast<T>(2);
        coeff[3] = (y[i + 1] - y[i]) / sig - (coeff[2] + p) * sig;
        coeff[4] = y[i];
    }
    c[5 * (n - 1)] = x[n - 1];
    c[5 * (n - 1) + 4] = y[n - 1];

    return;
}

template<class T>
T CubicSpline<T>::interpolate(const T x) const
{
    const Size i = index(x);
    const T    dx = x - c[i];

    return ((c[i + 1] * dx + c[i + 2]) * dx + c[i + 3]) * dx + c[i + 4];
}

template<class T>
T CubicSpline<T>::interpolateDerivative(const T x) const
{
    const Size i = index(x);
    const T    dx = x - c[i];

    return (3 * c[i + 1] * dx + 2 * c[i + 2]) * dx + c[i + 3];
}

template<class T>
std::vector<T> CubicSpline<T>::interpolate(const std::vector<T>& x) const
{
    std::vector<T> fx(x.size());
    for (auto i = 0; i < x.size(); ++i) { fx[i] = interpolateValue(x[i]); }

    return fx;
}

template<class T>
std::vector<T> CubicSpline<T>::interpolateDerivative(const std::vector<T>& x) const
{
    std::vector<T> dfx(x.size());
    for (auto i = 0; i < x.size(); ++i) { dfx[i] = interpolateValueDerivative(x[i]); }

    return dfx;
}

template<class T>
void CubicSpline<T>::interpolateWithDerivative(const T x, T& fx, T& dfx) const
{
    const Size i = index(x);
    const T    dx = x - c[i];
    fx = ((c[i + 1] * dx + c[i + 2]) * dx + c[i + 3]) * dx + c[i + 4];
    dfx = (3 * c[i + 1] * dx + 2 * c[i + 2]) * dx + c[i + 3];

    return;
}

template<class T>
void CubicSpline<T>::interpolateWithDerivative(const std::vector<T>& x,
                                               std::vector<T>&       fx,
                                               std::vector<T>&       dfx) const
{
    fx.resize(x.size());
    dfx.resize(x.size());
    for (auto i = 0; i < x.size(); ++i) { computeSplineAndDerivative(x[i], fx[i], dfx[i]); }

    return;
}

template<class T>
Size CubicSpline<T>::index(const T x) const
{
    // Casting to Size comes with the advantage that (not too huge) negative values are also mapped
    // to huge values.
    Size i = static_cast<Size>(x / delta) - 1;
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L2(
        if (i > n - 1)
        {
            global::tutor.bug(str("Spline index out of range (i = %zu, max = %zu), cannot "
                                  "evalute cubic spline interpolation.",
                                  i,
                                  n - 1));
        }
    );
#endif
    // Any bin index >= n will be mapped to n which is the out-of-range bin.
    i = std::min(i, n);

    return 5 * i;
}

} // namespace vaspml

#endif
