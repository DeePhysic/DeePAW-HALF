#include "math.hpp"

#include <cmath>      // for sqrt, pow, exp, signbit, cos, sin, cosh, sinh, abs
#include <cstdlib>    // for abs
#include <stdexcept>  // for runtime_error

using namespace vaspml;

void math::factorial(Vec1Real& f)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (f.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute factorials:"
                                     " output vector is not allocated.");
        }
    );

    // First entry is 0! = 1.
    f[0] = 1.0;
    for (Size i = 1; i < f.size(); ++i) f[i] = i * f[i - 1];

    return;
}

void math::doubleFactorialOdd(Vec1Real& f)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (f.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute double factorials:"
                                     " output vector is not allocated.");
        }
    );

    f[0] = 1.0;
    for (Size i = 1; i < f.size(); ++i) f[i] = (2 * i + 1) * f[i - 1];

    return;
}

Real math::clebschGordan(Int j1, Int j2, Int J, Int m1, Int m2, Int M, const Vec1Real& factorials)
{
    if (M != m1 + m2) return 0.0;
    Int k1 = j1 + j2 - J;
    if (k1 < 0) return 0.0;
    Int k2 = J + j1 - j2;
    if (k2 < 0) return 0.0;
    Int k3 = J + j2 - j1;
    if (k3 < 0) return 0.0;
    Int  k4 = j1 + j2 + J + 1;
    Real t = (2 * J + 1) * factorials[k1] * factorials[k2] * factorials[k3] / factorials[k4];
    k1 = j1 + m1;
    k2 = j1 - m1;
    k3 = j2 + m2;
    k4 = j2 - m2;
    Int       k5 = J + M;
    const Int k6 = J - M;
    t = std::sqrt(t * factorials[k1] * factorials[k2] * factorials[k3] * factorials[k4]
                  * factorials[k5] * factorials[k6]);
    const Int n1 = std::max({j2 - J - m1, j1 - J + m2, 0}) + 1;
    const Int n2 = std::min({j1 + j2 - J, j1 - m1, j2 + m2}) + 1;
    if (n1 > n2) return 0.0;
    Real t1 = 0.0;
    Int  m;
    for (m = n1; m <= n2; ++m)
    {
        const Int n = m - 1;
        k1 = j1 + j2 - J - n;
        k2 = j1 - m1 - n;
        k3 = j2 + m2 - n;
        k4 = J - j2 + m1 + n;
        k5 = J - j1 - m2 + n;
        // Note: intentional integer division!
        t1 += (1 + 4 * (Int)(n / 2) - 2 * n)
            / (factorials[n] * factorials[k1] * factorials[k2] * factorials[k3] * factorials[k4]
               * factorials[k5]);
    }

    return t * t1;
}

Real math::clebschGordanM0(Int j1, Int j2, Int J, const Vec1Real& factorials)
{
    if (std::abs(j1 - j2) > J) return 0.0;
    if (J > j1 + j2) return 0.0;
    Int jt = j1 + j2 + J;
    if (jt % 2) return 0.0;
    Int  p = jt / 2;
    Real result = std::sqrt((2 * J + 1.0) / (jt + 1.0));
    result *= factorials[p] / std::sqrt(factorials[2 * p]);
    Int x = p - j1;
    result *= std::sqrt(factorials[2 * x]) / factorials[x];
    x = p - j2;
    result *= std::sqrt(factorials[2 * x]) / factorials[x];
    x = p - J;
    result *= std::sqrt(factorials[2 * x]) / factorials[x];
    if (x % 2) return -result;
    else return result;
}

void math::gaussian(const Vec1Real& x, Vec1Real& f, Real width)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (x.size() != f.size())
        {
            throw std::runtime_error("ERROR: Cannot compute Gaussians:"
                                     " output vector has incorrect size.");
        }
    );

    std::transform(x.begin(),
                   x.end(),
                   f.begin(),
                   [&width](const Real& xi) { return std::exp(-width * xi * xi); });

    return;
}

void math::gaussianAndDerivative(const Vec1Real& x, Vec1Real& f, Vec1Real& df, Real width)
{
    VASPML_DEBUG_L1(
        // Check if output vectors have correct size.
        if (x.size() != f.size() || x.size() != df.size())
        {
            throw std::runtime_error("ERROR: Cannot compute Gaussians and derivatives:"
                                     " output vectors have incorrect sizes.");
        }
    );

    std::transform(x.begin(),
                   x.end(),
                   f.begin(),
                   [&width](const Real& xi) { return std::exp(-width * xi * xi); });
    width *= -2.0;
    std::transform(x.begin(),
                   x.end(),
                   f.begin(),
                   df.begin(),
                   [&width](const Real& xi, const Real& fi) { return width * xi * fi; });
    // Alternatively, though a bit slower:
    //for (Size i = 0; i < x.size(); ++i)
    //{
    //    f[i] = std::exp(-width * x[i] * x[i]);
    //    df[i] = -2.0 * width * x[i] * f[i];
    //}

    return;
}

Real math::sphericalBesselIterative(const Real x, const UInt nu)
{
    const Real x2 = 0.5 * x * x;
    const Real minDifference = 1.0E-16;
    const Int  maxIterations = 40;

    Real factor = 1.0;
    for (UInt inu = 1; inu <= nu; inu++) factor *= x / (2 * inu + 1);

    Real sum = 1.0;
    Real d = 1.0;
    Int  i;
    for (i = 1; i <= maxIterations; ++i)
    {
        d = -d * x2 / (i * (2 * (nu + i) + 1));
        if (std::abs(d) < minDifference) break;
        sum += d;
    }

    VASPML_DEBUG_L1(
        if (i >= maxIterations + 1)
        {
            throw std::runtime_error("ERROR: Iterative spherical Bessel function calculator does "
                                     "not converge.");
        }
    );

    return factor * sum;
}

Real math::sphericalBessel(const Real x, const UInt nu)
{
    VASPML_DEBUG_L1(
        // Sanity-check if argument is below zero.
        if (x < 0.0)
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(x)
                                     + ") to spherical Bessel function is below zero.");
        }
    );

    const Real iterativeSolutionLimit = 1.0;
    Real       jnu;
    if (x < iterativeSolutionLimit) { jnu = sphericalBesselIterative(x, nu); }
    else
    {
        // nth-order spherical Bessel function.
        jnu = std::sin(x) / x;
        if (nu == 0) return jnu;
        // (n+1)th-order spherical Bessel function.
        Real jnup1 = (-std::cos(x) + jnu) / x;
        if (nu == 1) return jnup1;
        for (UInt inu = 2; inu <= nu; ++inu)
        {
            jnu = -jnu + (2 * inu - 1) * jnup1 / x;
            if (inu == nu) break;
            std::swap(jnu, jnup1);
        }
    }

    return jnu;
}

void math::sphericalBesselAndDerivative(const Real x, Real& jnu, Real& djnu, const UInt nu)
{
    VASPML_DEBUG_L1(
        // Sanity-check if argument is below zero.
        if (x < 0.0)
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(x)
                                     + ") to spherical Bessel function is below zero.");
        }
    );

    const Real iterativeSolutionLimit = 1.0;
    if (x < iterativeSolutionLimit)
    {
        jnu = sphericalBesselIterative(x, nu);
        djnu = sphericalBesselIterative(x, nu + 1);
        djnu = -djnu + jnu * nu / x;
    }
    else
    {
        // nth-order spherical Bessel function.
        jnu = std::sin(x) / x;
        // (n+1)th-order spherical Bessel function.
        djnu = (-std::cos(x) + jnu) / x;
        for (UInt inu = 2; inu <= nu + 1; ++inu)
        {
            const Real tmp = -jnu + (2 * inu - 1) * djnu / x;
            jnu = djnu;
            djnu = tmp;
        }
        djnu = -djnu + jnu * nu / x;
    }

    return;
}

void math::sphericalBessel(const Vec1Real& x, Vec1Real& jnu, const UInt nu)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (x.size() != jnu.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical Bessel function:"
                                     " output vector has incorrect size.");
        }
        // Sanity-check if any radius is below zero.
        auto xtest = std::find_if(x.begin(), x.end(), [](Real xi) { return xi < 0.0; });
        if (xtest != x.end())
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(*xtest)
                                     + ") to spherical Bessel function is below zero.");
        }
    );

    std::transform(x.begin(),
                   x.end(),
                   jnu.begin(),
                   [&nu](const Real& xi) { return sphericalBessel(xi, nu); });
    // Alternatively, but way slower:
    //[&nu](const Real& xi) { return std::sph_bessel(nu, xi); });

    return;
}

void math::sphericalBesselAndDerivative(const Vec1Real& x,
                                        Vec1Real&       jnu,
                                        Vec1Real&       djnu,
                                        const UInt      nu)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (x.size() != jnu.size() || x.size() != djnu.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical Bessel function:"
                                     " output vector has incorrect size.");
        }
        // Sanity-check if any radius is below zero.
        auto xtest = std::find_if(x.begin(), x.end(), [](Real xi) { return xi < 0.0; });
        if (xtest != x.end())
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(*xtest)
                                     + ") to spherical Bessel function is below zero.");
        }
    );

    for (Size i = 0; i < x.size(); ++i) { sphericalBesselAndDerivative(x[i], jnu[i], djnu[i], nu); }

    return;
}

void math::sphericalBesselRoots(Vec1Real& roots, const UInt nu)
{
    VASPML_DEBUG_L1(
        // Check if output vector has non-zero size.
        if (roots.size() == 0)
        {
            throw std::runtime_error("ERROR: Cannot compute spherical Bessel roots:"
                                     " output vector has size zero on entry (should be number"
                                     " of desired roots).");
        }
    );

    const Real breakCondition = 1.0E-10;
    const Real step = 0.1;

    Real x = step;
    Size n = 0;
    // For interval bisectioning.
    Real delta;
    Real ix;
    bool foundRoot;

    do {
        Real j1 = sphericalBessel(x, nu);
        foundRoot = false;
        do {
            // Coarse search, increase second x-value by step each time.
            x += step;
            Real j2 = sphericalBessel(x, nu);
            // Check if sign of function value differs.
            if (std::signbit(j1) != std::signbit(j2))
            {
                // Interval bisectioning.
                delta = step;
                ix = x;
                do {
                    delta /= 2.0;
                    if (std::signbit(j1) != std::signbit(j2)) ix -= delta;
                    else ix += delta;
                    j2 = sphericalBessel(ix, nu);
                }
                while (delta > breakCondition);
                // Save result in output vector.
                roots[n] = ix;
                n++;
                foundRoot = true;
            }
        }
        while (!foundRoot);
    }
    while (n < roots.size());

    return;
}

void math::modifiedSphericalBesselAndDerivative(const Vec1Real& x, Vec2Real& inu, Vec2Real& dinu)
{
    VASPML_DEBUG_L1(
        // Check if output vector is allocated.
        if (inu.size() <= 0 || dinu.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute modified spherical Bessel function:"
                                     " output vectors are not allocated.");
        }
        // Check if input vectors are allocated.
        if (x.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute modified spherical Bessel function:"
                                     " input vector is empty.");
        }
        // Check if output vectors have equal size.
        if (inu.size() != dinu.size())
        {
            throw std::runtime_error("ERROR: Cannot compute modified spherical Bessel function:"
                                     " output vectors have different size.");
        }
        // Check if output vector has correct size.
        for (Size i = 0; i < inu.size(); ++i)
        {
            if (inu.at(i).size() != x.size() || dinu.at(i).size() != x.size())
            {
                throw std::runtime_error("ERROR: Cannot compute modified spherical Bessel function:"
                                         "output vectors have incorrect size (index 2).");
            }
        }
        // Sanity-check if any argument is below zero.
        auto xtest = std::find_if(x.begin(), x.end(), [](Real xi) { return xi < 0.0; });
        if (xtest != x.end())
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(*xtest)
                                     + ") to modified spherical Bessel function is below zero.");
        }
    );

    const Size lMax = inu.size() - 1;
    const Real eps = 1.0E-8;
    const Real eps2 = std::sqrt(eps);
    const Real eps3 = std::pow(eps, 1.0 / 3.0);
    Vec1Real   fac(lMax + 1);
    math::doubleFactorialOdd(fac);

    for (Size i = 0; i < x.size(); ++i)
    {
        const Real xi = x[i];
        const Real xi2 = xi * xi;
        const Real sinhx = std::sinh(x[i]);
        const Real coshx = std::cosh(x[i]);
        // lMax == 0.
        if (xi > eps2)
        {
            inu[0][i] = sinhx / xi;
            dinu[0][i] = (xi * coshx - sinhx) / xi2;
        }
        else
        {
            inu[0][i] = 1.0;
            dinu[0][i] = 0.0;
        }
        if (lMax == 0) continue;

        // lMax == 1.
        if (xi > eps3)
        {
            inu[1][i] = (xi * coshx - sinhx) / xi2;
            dinu[1][i] = inu[0][i] - 3.0 * inu[1][i] / xi;
            dinu[1][i] = (inu[0][i] + 2.0 * dinu[1][i]) / 3.0;
        }
        else
        {
            inu[1][i] = xi / fac[1];
            dinu[1][i] = 1.0 / fac[1];
        }
    }

    // lMax > 1.
    for (Size l = 2; l <= lMax; ++l)
    {
        const Real epsl = std::pow(eps, 1.0 / (l + 2));
        for (Size i = 0; i < x.size(); ++i)
        {
            if (x[i] > epsl) inu[l][i] = inu[l - 2][i] - (2 * (l - 1) + 1) / x[i] * inu[l - 1][i];
            else inu[l][i] = std::pow(x[i], l) / fac[l];
        }
    }
    for (Size l = 2; l <= lMax; ++l)
    {
        const Real epsl = std::pow(eps, 1.0 / (l + 2));
        for (Size i = 0; i < x.size(); ++i)
        {
            if (x[i] > epsl)
            {
                dinu[l][i] = inu[l - 1][i] - (2 * l + 1) / x[i] * inu[l][i];
                dinu[l][i] = (l * inu[l - 1][i] + (l + 1) * dinu[l][i]) / (2 * l + 1);
            }
            else { dinu[l][i] = l * std::pow(x[i], l - 1) / fac[l]; }
        }
    }

    return;
}

void math::expModSphBessel(const Real a2, const Vec1Real& ab, const Vec1Real& b2, Vec2Real& p)
{
    VASPML_DEBUG_L1(
        // Check if output vector is allocated.
        if (p.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute expModSphBessel function:"
                                     " output vector is not allocated.");
        }
        // Check if input vectors are allocated.
        if (ab.size() <= 0 || b2.size() <= 0)
        {
            throw std::runtime_error("ERROR: Cannot compute expModSphBessel function:"
                                     " input vectors are empty.");
        }
        // Check if input vectors have equal size.
        if (ab.size() != b2.size())
        {
            throw std::runtime_error("ERROR: Cannot compute expModSphBessel function:"
                                     " input vectors have different size.");
        }
        // Check if output vector has correct size.
        for (auto const& pi : p)
        {
            if (pi.size() != ab.size())
            {
                throw std::runtime_error("ERROR: Cannot compute expModSphBesselfunction: "
                                         "output vector has incorrect size (index 2).");
            }
        }
        // Sanity-check if any argument is below zero.
        if (a2 < 0.0)
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(a2)
                                     + ") to expModSphBessel function is below zero (a2).");
        }
        auto abtest = std::find_if(ab.begin(), ab.end(), [](Real xi) { return xi < 0.0; });
        if (abtest != ab.end())
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(*abtest)
                                     + ") to expModSphBessel function is below zero (ab).");
        }
        auto b2test = std::find_if(b2.begin(), b2.end(), [](Real xi) { return xi < 0.0; });
        if (b2test != b2.end())
        {
            throw std::runtime_error("ERROR: Argument (" + std::to_string(*b2test)
                                     + ") to expModSphBessel function is below zero (b2).");
        }
    );

    const Size lMax = p.size() - 1;
    const Real eps = 1.0E-8;
    const Real eps2 = std::sqrt(eps);
    const Real eps3 = std::pow(eps, 1.0 / 3.0);
    Vec1Real   fac(lMax + 1);
    math::doubleFactorialOdd(fac);
    using std::exp;

    for (Size i = 0; i < ab.size(); ++i)
    {
        const Real a2b2 = a2 + b2[i];
        const Real expPlus = exp(ab[i] - a2b2);
        const Real expMinus = exp(-ab[i] - a2b2);
        // lMax == 0.
        if (ab[i] > eps2) p[0][i] = 0.5 * (expPlus - expMinus) / ab[i];
        else p[0][i] = exp(-a2b2);
        if (lMax == 0) continue;

        // lMax == 1.
        if (ab[i] > eps3)
        {
            p[1][i] =
                (0.5 * ab[i] * (expPlus + expMinus) - 0.5 * (expPlus - expMinus)) / (ab[i] * ab[i]);
        }
        else p[1][i] = ab[i] / fac[1] * exp(-a2b2);
    }

    // lMax > 1.
    for (Size l = 2; l <= lMax; ++l)
    {
        const Real epsl = std::pow(eps, 1.0 / (l + 2));
        for (Size i = 0; i < ab.size(); ++i)
        {
            if (ab[i] > epsl) p[l][i] = p[l - 2][i] - (2 * (l - 1) + 1) / ab[i] * p[l - 1][i];
            else p[l][i] = std::pow(ab[i], l) * exp(-a2 - b2[i]) / fac[l];
        }
    }

    return;
}
