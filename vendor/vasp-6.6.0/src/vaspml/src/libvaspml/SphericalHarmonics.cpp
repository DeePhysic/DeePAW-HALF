#include "SphericalHarmonics.hpp"

#include "ParallelEnvironment.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace vaspml;
using Integrate3SphericalHarmonicsData = math::SphericalHarmonics::Integrate3SphericalHarmonicsData;

math::SphericalHarmonics::SphericalHarmonics(Int lmax) : lmax(lmax), ldim((lmax + 1) * (lmax + 1))
{
    if (lmax > lhc) integrate3SphericalHarmonics(lmax, this->i3sh);
}

VASPML_EXEC_SPACE_SPECIFIER
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradientsNormalized(const Vec1Real& rhat,
                                                                               Vec1Real&       ylm,
                                                                               Vec1Real& ylmd) const
{
    Size n = rhat.size() / 3;
    computeSphericalHarmonicsAndGradientsNormalized(rhat.data(), n, ylm.data(), ylmd.data());
}

VASPML_EXEC_SPACE_SPECIFIER
void math::SphericalHarmonics::computeSphericalHarmonicAndGradient(const Real& x,
                                                                   const Real& y,
                                                                   const Real& z,
                                                                   const Real& rnorm,
                                                                   const Size& l,
                                                                   const Int&  m,
                                                                   Real&       ylm,
                                                                   Real&       ylmdx,
                                                                   Real&       ylmdy,
                                                                   Real&       ylmdz) const
{
    if (l > (Size)lhc)
    {
#ifdef VASPML_PALGO_GPU
        return;
#else
        global::tutor.bug("ERROR: Single spherical harmonics function not implemented "
                          "for l > "
                          + std::to_string(lhc) + ".");
#endif
    }

    if (l == 0)
    {
        ylm = 1.0 / std::sqrt(constants::PI4);
        ylmdx = 0.0;
        ylmdy = 0.0;
        ylmdz = 0.0;
    }
    else if (l == 1)
    {
        const Real pres3 = std::sqrt(3.0 / constants::PI4);
        if (m == -1)
        {
            ylm = pres3 * y;
            ylmdx = 0.0;
            ylmdy = pres3;
            ylmdz = 0.0;
        }
        else if (m == 0)
        {
            ylm = pres3 * z;
            ylmdx = 0.0;
            ylmdy = 0.0;
            ylmdz = pres3;
        }
        else if (m == 1)
        {
            ylm = pres3 * x;
            ylmdx = pres3;
            ylmdy = 0.0;
            ylmdz = 0.0;
        }
    }
    else if (l == 2)
    {
        const Real pres15 = std::sqrt(15.0 / constants::PI4);
        const Real pres5h = std::sqrt(5.0 / constants::PI4) / 2.0;
        if (m == -2)
        {
            ylm = pres15 * x * y;
            ylmdx = pres15 * y;
            ylmdy = pres15 * x;
            ylmdz = 0.0;
        }
        else if (m == -1)
        {
            ylm = pres15 * y * z;
            ylmdx = 0.0;
            ylmdy = pres15 * z;
            ylmdz = pres15 * y;
        }
        else if (m == 0)
        {
            ylm = pres5h * (3.0 * z * z - 1.0);
            ylmdx = 0.0;
            ylmdy = 0.0;
            ylmdz = pres5h * 6.0 * z;
        }
        else if (m == 1)
        {
            ylm = pres15 * x * z;
            ylmdx = pres15 * z;
            ylmdy = 0.0;
            ylmdz = pres15 * x;
        }
        else if (m == 2)
        {
            ylm = pres15 * 0.5 * (x * x - y * y);
            ylmdx = pres15 * x;
            ylmdy = -pres15 * y;
            ylmdz = 0.0;
        }
    }
    else if (l == 3)
    {
        const Real pres7 = std::sqrt(7.0 / constants::PI) / 4.0;
        const Real pres21 = std::sqrt(21.0 / constants::PI2) / 4.0;
        const Real pres35 = std::sqrt(35.0 / constants::PI2) / 4.0;
        const Real pres105 = std::sqrt(105.0 / constants::PI) / 4.0;
        const Real x2 = x * x;
        const Real y2 = y * y;
        const Real z2 = z * z;
        const Real r2 = x2 + y2 + z2;
        if (m == -3)
        {
            ylm = pres35 * y * (3.0 * x2 - y2);
            ylmdx = 6.0 * pres35 * x * y;
            ylmdy = 3.0 * pres35 * (x2 - y2);
            ylmdz = 0.0;
        }
        else if (m == -2)
        {
            ylm = 2.0 * pres105 * x * y * z;
            ylmdx = 2.0 * pres105 * y * z;
            ylmdy = 2.0 * pres105 * x * z;
            ylmdz = 2.0 * pres105 * x * y;
        }
        else if (m == -1)
        {
            ylm = pres21 * y * (5.0 * z2 - r2);
            ylmdx = -2.0 * pres21 * x * y;
            ylmdy = pres21 * (4.0 * z2 - x2 - 3.0 * y2);
            ylmdz = 8.0 * pres21 * y * z;
        }
        else if (m == 0)
        {
            ylm = pres7 * (5.0 * z * z2 - 3.0 * z * r2);
            ylmdx = -6.0 * pres7 * x * z;
            ylmdy = -6.0 * pres7 * y * z;
            ylmdz = 3.0 * pres7 * (3.0 * z2 - r2);
        }
        else if (m == 1)
        {
            ylm = pres21 * x * (5.0 * z2 - r2);
            ylmdx = pres21 * (4.0 * z2 - y2 - 3.0 * x2);
            ylmdy = -2.0 * pres21 * x * y;
            ylmdz = 8.0 * pres21 * x * z;
        }
        else if (m == 2)
        {
            ylm = pres105 * (x2 - y2) * z;
            ylmdx = 2.0 * pres105 * x * z;
            ylmdy = -2.0 * pres105 * y * z;
            ylmdz = pres105 * (x2 - y2);
        }
        else if (m == 3)
        {
            ylm = pres35 * x * (x2 - 3.0 * y2);
            ylmdx = 3.0 * pres35 * (x2 - y2);
            ylmdy = -6.0 * pres35 * x * y;
            ylmdz = 0.0;
        }
    }
    else if (l == 4)
    {
        const Real pres1 = 0.1875 / std::sqrt(constants::PI);
        const Real pres5a = std::sqrt(5.0 / constants::PI) * 0.75;
        const Real pres5b = std::sqrt(5.0 / constants::PI2) * 0.75;
        const Real pres35a = std::sqrt(35.0 / constants::PI) * 0.75;
        const Real pres35b = std::sqrt(35.0 / constants::PI2) * 0.75;
        const Real x2 = x * x;
        const Real y2 = y * y;
        const Real z2 = z * z;
        const Real r2 = x2 + y2 + z2;
        if (m == -4)
        {
            ylm = pres35a * x * y * (x2 - y2);
            ylmdx = pres35a * y * (3 * x2 - y2);
            ylmdy = pres35a * x * (x2 - 3 * y2);
            ylmdz = 0;
        }
        else if (m == -3)
        {
            ylm = pres35b * y * z * (3 * x2 - y2);
            ylmdx = 6 * pres35b * x * y * z;
            ylmdy = pres35b * z * (3 * x2 - 3 * y2);
            ylmdz = pres35b * y * (3 * x2 - y2);
        }
        else if (m == -2)
        {
            ylm = pres5a * x * y * (-r2 + 7 * z2);
            ylmdx = pres5a * y * (-r2 + 7 * z2);
            ylmdy = pres5a * x * (-r2 + 7 * z2);
            ylmdz = 14 * pres5a * x * y * z;
        }
        else if (m == -1)
        {
            ylm = pres5b * y * z * (-3 * r2 + 7 * z2);
            ylmdx = 0;
            ylmdy = pres5b * z * (-3 * r2 + 7 * z2);
            ylmdz = 3 * pres5b * y * (-r2 + 7 * z2);
        }
        else if (m == 0)
        {
            ylm = pres1 * (-30 * r2 * z2 + 3 * (r2 * r2) + 35 * (z2 * z2));
            ylmdx = 0;
            ylmdy = 0;
            ylmdz = 20 * pres1 * z * (-3 * r2 + 7 * z2);
        }
        else if (m == 1)
        {
            ylm = pres5b * x * z * (-3 * r2 + 7 * z2);
            ylmdx = pres5b * z * (-3 * r2 + 7 * z2);
            ylmdy = 0;
            ylmdz = 3 * pres5b * x * (-r2 + 7 * z2);
        }
        else if (m == 2)
        {
            ylm = 0.5 * pres5a * (-r2 + 7 * z2) * (x2 - y2);
            ylmdx = 1.0 * pres5a * x * (-r2 + 7 * z2);
            ylmdy = -1.0 * pres5a * y * (-r2 + 7 * z2);
            ylmdz = 7.0 * pres5a * z * (x2 - y2);
        }
        else if (m == 3)
        {
            ylm = pres35b * x * z * (x2 - 3 * y2);
            ylmdx = pres35b * z * (3 * x2 - 3 * y2);
            ylmdy = -6 * pres35b * x * y * z;
            ylmdz = pres35b * x * (x2 - 3 * y2);
        }
        else if (m == 4)
        {
            ylm = 0.25 * pres35a * (x2 * (x2 - 3 * y2) - y2 * (3 * x2 - y2));
            ylmdx = 0.5 * pres35a * x * (2 * x2 - 6 * y2);
            ylmdy = 0.5 * pres35a * y * (-6 * x2 + 2 * y2);
            ylmdz = 0;
        }
    }

    // Finally, compute derivative for unnormalized input vector.
    const Real invrnorm = 1.0 / rnorm;
    const Real tmpx = ylmdx;
    const Real tmpy = ylmdy;
    const Real tmpz = ylmdz;
    const Real tmps = tmpx * x + tmpy * y + tmpz * z;
    ylmdx = (tmpx - tmps * x) * invrnorm;
    ylmdy = (tmpy - tmps * y) * invrnorm;
    ylmdz = (tmpz - tmps * z) * invrnorm;

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradientsNormalized(const Real* rhat,
                                                                               const Size& nPoints,
                                                                               Real*       ylm,
                                                                               Real* ylmd) const
{
    if (lmax < 0) return;
    const Real pre = 1.0 / std::sqrt(constants::PI4);
    const Size n = nPoints;

    // l = 0
    for (Size i = 0; i < n * ldim; i += ldim)
    {
        ylm[i] = pre;
        ylmd[3 * i] = 0.0;
        ylmd[3 * i + 1] = 0.0;
        ylmd[3 * i + 2] = 0.0;
    }
    if (lmax < 1) return;

    // l = 1
    const Real pres3 = std::sqrt(3.0 / constants::PI4);
    for (Size i = 0; i < n; ++i)
    {
        Size       s = i * ldim + 1;
        const Size a = 3 * i;
        // clang-format off
        ylm[s    ] = pres3 * rhat[a + 1]; // l = 1, m = -1
        ylm[s + 1] = pres3 * rhat[a + 2]; // l = 1, m = 0
        ylm[s + 2] = pres3 * rhat[a    ]; // l = 1, m = 1
        s = 3 * i * ldim + 3;
        // l = 1, m = -1
        ylmd[s    ] = 0.0;
        ylmd[s + 1] = pres3;
        ylmd[s + 2] = 0.0;
        // l = 1, m = 0
        ylmd[s + 3] = 0.0;
        ylmd[s + 4] = 0.0;
        ylmd[s + 5] = pres3;
        // l = 1, m = 1
        ylmd[s + 6] = pres3;
        ylmd[s + 7] = 0.0;
        ylmd[s + 8] = 0.0;
        // clang-format on
    }
    if (lmax < 2) return;

    // l = 2
    const Real pres15 = std::sqrt(15.0 / constants::PI4);
    const Real pres5h = std::sqrt(5.0 / constants::PI4) / 2.0;
    for (Size i = 0; i < n; ++i)
    {
        Size       s = i * ldim + 4;
        const Real x = rhat[3 * i];
        const Real y = rhat[3 * i + 1];
        const Real z = rhat[3 * i + 2];
        // clang-format off
        ylm[s    ] = pres15 * x * y;                 // l = 2, m = -2
        ylm[s + 1] = pres15 * y * z;                 // l = 2, m = -1
        ylm[s + 2] = pres5h * (3.0 * z * z - 1.0);   // l = 2, m = 0
        ylm[s + 3] = pres15 * x * z;                 // l = 2, m = 1
        ylm[s + 4] = pres15 * 0.5 * (x * x - y * y); // l = 2, m = 2
        s = 3 * i * ldim + 12;
        // l = 2, m = -2
        ylmd[s     ] = pres15 * y;
        ylmd[s +  1] = pres15 * x;
        ylmd[s +  2] = 0.0;
        // l = 2, m = -1
        ylmd[s +  3] = 0.0;
        ylmd[s +  4] = pres15 * z;
        ylmd[s +  5] = pres15 * y;
        // l = 2, m = 0
        ylmd[s +  6] = 0.0;
        ylmd[s +  7] = 0.0;
        ylmd[s +  8] = pres5h * 6.0 * z;
        // l = 2, m = 1
        ylmd[s +  9] = pres15 * z;
        ylmd[s + 10] = 0.0;
        ylmd[s + 11] = pres15 * x;
        // l = 2, m = 2
        ylmd[s + 12] = pres15 * x;
        ylmd[s + 13] = -pres15 * y;
        ylmd[s + 14] = 0.0;
        // clang-format on
    }
    if (lmax < 3) return;

    // l = 3
    const Real pres7 = std::sqrt(7.0 / constants::PI) / 4.0;
    const Real pres21 = std::sqrt(21.0 / constants::PI2) / 4.0;
    const Real pres35 = std::sqrt(35.0 / constants::PI2) / 4.0;
    const Real pres105 = std::sqrt(105.0 / constants::PI) / 4.0;
    for (Size i = 0; i < n; ++i)
    {
        Size       s = i * ldim + 9;
        const Real x = rhat[3 * i];
        const Real y = rhat[3 * i + 1];
        const Real z = rhat[3 * i + 2];
        const Real x2 = x * x;
        const Real y2 = y * y;
        const Real z2 = z * z;
        const Real r2 = x2 + y2 + z2;
        // clang-format off
        ylm[s    ] = pres35 * y * (3.0 * x2 - y2);          // l = 3, m = -3
        ylm[s + 1] = 2.0 * pres105 * x * y * z;             // l = 3, m = -2
        ylm[s + 2] = pres21 * y * (5.0 * z2 - r2);          // l = 3, m = -1
        ylm[s + 3] = pres7 * (5.0 * z * z2 - 3.0 * z * r2); // l = 3, m = 0
        ylm[s + 4] = pres21 * x * (5.0 * z2 - r2);          // l = 3, m = 1
        ylm[s + 5] = pres105 * (x2 - y2) * z;               // l = 3, m = 2
        ylm[s + 6] = pres35 * x * (x2 - 3.0 * y2);          // l = 3, m = 3
        s = 3 * i * ldim + 27;
        // l = 3, m = -3
        ylmd[s     ] = 6.0 * pres35 * x * y;
        ylmd[s +  1] = 3.0 * pres35 * (x2 - y2);
        ylmd[s +  2] = 0.0;
        // l = 3, m = -2
        ylmd[s +  3] = 2.0 * pres105 * y * z;
        ylmd[s +  4] = 2.0 * pres105 * x * z;
        ylmd[s +  5] = 2.0 * pres105 * x * y;
        // l = 3, m = -1
        ylmd[s +  6] = -2.0 * pres21 * x * y;
        ylmd[s +  7] = pres21 * (4.0 * z2 - x2 - 3.0 * y2);
        ylmd[s +  8] = 8.0 * pres21 * y * z;
        // l = 3, m = 0
        ylmd[s +  9] = -6.0 * pres7 * x * z;
        ylmd[s + 10] = -6.0 * pres7 * y * z;
        ylmd[s + 11] = 3.0 * pres7 * (3.0 * z2 - r2);
        // l = 3, m = 1
        ylmd[s + 12] = pres21 * (4.0 * z2 - y2 - 3.0 * x2);
        ylmd[s + 13] = -2.0 * pres21 * x * y;
        ylmd[s + 14] = 8.0 * pres21 * x * z;
        // l = 3, m = 2
        ylmd[s + 15] = 2.0 * pres105 * x * z;
        ylmd[s + 16] = -2.0 * pres105 * y * z;
        ylmd[s + 17] = pres105 * (x2 - y2);
        // l = 3, m = 3
        ylmd[s + 18] = 3.0 * pres35 * (x2 - y2);
        ylmd[s + 19] = -6.0 * pres35 * x * y;
        ylmd[s + 20] = 0.0;
        // clang-format on
    }
    if (lmax < 4) return;

    // l = 4 Note: derivatives generated via sympy.
    const Real pres1 = 0.1875 / std::sqrt(constants::PI);
    const Real pres5a = std::sqrt(5.0 / constants::PI) * 0.75;
    const Real pres5b = std::sqrt(5.0 / constants::PI2) * 0.75;
    const Real pres35a = std::sqrt(35.0 / constants::PI) * 0.75;
    const Real pres35b = std::sqrt(35.0 / constants::PI2) * 0.75;
    for (Size i = 0; i < n; ++i)
    {
        Size       s = i * ldim + 16;
        const Real x = rhat[3 * i];
        const Real y = rhat[3 * i + 1];
        const Real z = rhat[3 * i + 2];
        const Real x2 = x * x;
        const Real y2 = y * y;
        const Real z2 = z * z;
        const Real r2 = x2 + y2 + z2;
        // clang-format off
        ylm[s    ] = pres35a * x * y * (x2 - y2);                                // l = 4, m = -4
        ylm[s + 1] = pres35b * y * z * (3 * x2 - y2);                            // l = 4, m = -3
        ylm[s + 2] = pres5a * x * y * (-r2 + 7 * z2);                            // l = 4, m = -2
        ylm[s + 3] = pres5b * y * z * (-3 * r2 + 7 * z2);                        // l = 4, m = -1
        ylm[s + 4] = pres1 * (-30 * r2 * z2 + 3 * (r2 * r2) + 35 * (z2 * z2));   // l = 4, m = 0
        ylm[s + 5] = pres5b * x * z * (-3 * r2 + 7 * z2);                        // l = 4, m = 1
        ylm[s + 6] = 0.5 * pres5a * (-r2 + 7 * z2) * (x2 - y2);                  // l = 4, m = 2
        ylm[s + 7] = pres35b * x * z * (x2 - 3 * y2);                            // l = 4, m = 3
        ylm[s + 8] = 0.25 * pres35a * (x2 * (x2 - 3 * y2) - y2 * (3 * x2 - y2)); // l = 4, m = 4
        s = 3 * i * ldim + 48;
        // l = 4, m = -4
        ylmd[s +  0] = pres35a * y * (3 * x2 - y2);
        ylmd[s +  1] = pres35a * x * (x2 - 3 * y2);
        ylmd[s +  2] = 0;
        // l = 4, m = -3
        ylmd[s +  3] = 6 * pres35b * x * y * z;
        ylmd[s +  4] = pres35b * z * (3 * x2 - 3 * y2);
        ylmd[s +  5] = pres35b * y * (3 * x2 - y2);
        // l = 4, m = -2
        ylmd[s +  6] = pres5a * y * (-r2 + 7 * z2);
        ylmd[s +  7] = pres5a * x * (-r2 + 7 * z2);
        ylmd[s +  8] = 14 * pres5a * x * y * z;
        // l = 4, m = -1
        ylmd[s +  9] = 0;
        ylmd[s + 10] = pres5b * z * (-3 * r2 + 7 * z2);
        ylmd[s + 11] = 3 * pres5b * y * (-r2 + 7 * z2);
        // l = 4, m = 0
        ylmd[s + 12] = 0;
        ylmd[s + 13] = 0;
        ylmd[s + 14] = 20 * pres1 * z * (-3 * r2 + 7 * z2);
        // l = 4, m = 1
        ylmd[s + 15] = pres5b * z * (-3 * r2 + 7 * z2);
        ylmd[s + 16] = 0;
        ylmd[s + 17] = 3 * pres5b * x * (-r2 + 7 * z2);
        // l = 4, m = 2
        ylmd[s + 18] = 1.0 * pres5a * x * (-r2 + 7 * z2);
        ylmd[s + 19] = -1.0 * pres5a * y * (-r2 + 7 * z2);
        ylmd[s + 20] = 7.0 * pres5a * z * (x2 - y2);
        // l = 4, m = 3
        ylmd[s + 21] = pres35b * z * (3 * x2 - 3 * y2);
        ylmd[s + 22] = -6 * pres35b * x * y * z;
        ylmd[s + 23] = pres35b * x * (x2 - 3 * y2);
        // l = 4, m = 4
        ylmd[s + 24] = 0.5 * pres35a * x * (2 * x2 - 6 * y2);
        ylmd[s + 25] = 0.5 * pres35a * y * (-6 * x2 + 2 * y2);
        ylmd[s + 26] = 0;
        // clang-format on
    }
    if (lmax < 5) return;

    // Early exit if i3sh data set is not prepared for higher l's.
    if (this->i3sh.lmax < lhc + 1) return;

    // Zero all remaining entries for l > lhc.
    for (Size i = 0; i < n; ++i)
    {
        for (Size ilm = (lhc + 1) * (lhc + 1); ilm < (Size)ldim; ++ilm)
        {
            ylm[i * ldim + ilm] = 0.0;
        }
        for (Size ilm = (lhc + 1) * (lhc + 1) * 3; ilm < 3 * (Size)ldim; ++ilm)
        {
            ylmd[3 * i * ldim + ilm] = 0.0;
        }
    }

    // For l > lhc we use (a variant of) Clebsch-Gordan coefficients to recursively obtain the
    // remaining spherical harmonics.
    Size l2 = 1;
    for (Size l1 = lhc; l1 < (Size)i3sh.lmax; ++l1)
    {
        Size       lmindx = i3sh.lookup[l1][l2];
        const Size lsum = l1 + l2;
        for (Size m1 = 0; m1 <= 2 * l1; ++m1)
        {
            for (Size m2 = 0; m2 <= 2 * l2; ++m2)
            {
                for (Int ic = i3sh.indcg[lmindx]; ic < i3sh.indcg[lmindx + 1]; ++ic)
                {
                    const Size lm = i3sh.js[ic];
                    if (lm + 1 > lsum * lsum && lm + 1 <= (lsum + 1) * (lsum + 1))
                    {
                        const Real c = i3sh.ylm3i[ic];
                        for (Size i = 0; i < n; ++i)
                        {
                            const Size s = i * ldim;
                            Size       s1 = s + l1 * l1 + m1;
                            Size       s2 = s + l2 * l2 + m2;
                            const Size slm = 3 * (s + lm);
                            const Real ylm1 = ylm[s1];
                            const Real ylm2 = ylm[s2];
                            s1 *= 3;
                            s2 *= 3;
                            ylm[s + lm] += c * ylm1 * ylm2;
                            // clang-format off
                            ylmd[slm    ] += c * (ylmd[s1    ] * ylm2 + ylmd[s2    ] * ylm1);
                            ylmd[slm + 1] += c * (ylmd[s1 + 1] * ylm2 + ylmd[s2 + 1] * ylm1);
                            ylmd[slm + 2] += c * (ylmd[s1 + 2] * ylm2 + ylmd[s2 + 2] * ylm1);
                            // clang-format on
                        }
                    }
                }
                lmindx++;
            }
        }
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradients(const Vec1Real& rhat,
                                                                     const Vec1Real& rnorm,
                                                                     Vec1Real&       ylm,
                                                                     Vec1Real&       ylmd) const
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (rhat.size() != 3 * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: input vectors "
                                     "do not match in size.");
        }
        if (ylm.size() != ldim * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: incorrect size "
                                     "of output vector."
                                     + std::to_string(ylm.size()) + " "
                                     + std::to_string(ldim * rnorm.size()));
        }
        if (ylmd.size() != 3 * ldim * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: incorrect size "
                                     "of derivative output vector.");
        }
    );
#endif

    // First, compute spherical harmonics and gradients for normalized vectors.
    computeSphericalHarmonicsAndGradientsNormalized(rhat, ylm, ylmd);

    // Finally, compute gradients with respect to unnormalized coordinates.
    for (Size i = 0; i < rnorm.size(); ++i)
    {
        Size const a = 3 * i;
        const Real x = rhat[a];
        const Real y = rhat[a + 1];
        const Real z = rhat[a + 2];
        const Real invrnorm = 1.0 / rnorm[i];
        Size       s = a * ldim;
        for (Size ilm = 0; ilm < (Size)ldim; ++ilm)
        {
            const Real ylmdx = ylmd[s];
            const Real ylmdy = ylmd[s + 1];
            const Real ylmdz = ylmd[s + 2];
            const Real ylmdxa = ylmdx * x + ylmdy * y + ylmdz * z;
            // clang-format off
            ylmd[s    ] = (ylmdx - ylmdxa * x) * invrnorm;
            ylmd[s + 1] = (ylmdy - ylmdxa * y) * invrnorm;
            ylmd[s + 2] = (ylmdz - ylmdxa * z) * invrnorm;
            s += 3;
            // clang-format on
        }
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradients(const Real* rhat,
                                                                     const Real* rnorm,
                                                                     const Size& nPoints,
                                                                     Real*       ylm,
                                                                     Real*       ylmd) const
{
    // First, compute spherical harmonics and gradients for normalized vectors.
    computeSphericalHarmonicsAndGradientsNormalized(rhat, nPoints, ylm, ylmd);

    // Finally, compute gradients with respect to unnormalized coordinates.
    for (Size i = 0; i < nPoints; ++i)
    {
        Size const a = 3 * i;
        const Real x = rhat[a];
        const Real y = rhat[a + 1];
        const Real z = rhat[a + 2];
        const Real invrnorm = 1.0 / rnorm[i];
        Size       s = a * ldim;
        for (Size ilm = 0; ilm < (Size)ldim; ++ilm)
        {
            const Real ylmdx = ylmd[s];
            const Real ylmdy = ylmd[s + 1];
            const Real ylmdz = ylmd[s + 2];
            const Real ylmdxa = ylmdx * x + ylmdy * y + ylmdz * z;
            // clang-format off
            ylmd[s    ] = (ylmdx - ylmdxa * x) * invrnorm;
            ylmd[s + 1] = (ylmdy - ylmdxa * y) * invrnorm;
            ylmd[s + 2] = (ylmdz - ylmdxa * z) * invrnorm;
            s += 3;
            // clang-format on
        }
    }

    return;
}

void math::SphericalHarmonics::integrate3SphericalHarmonics(Int                               lmax,
                                                            Integrate3SphericalHarmonicsData& data)
{
    VASPML_DEBUG_L1(
        if (data.ylm3.size() != 0 || data.ylm3i.size() != 0 || data.jl.size() != 0
            || data.js.size() != 0 || data.indcg.size() != 0 || data.lookup.size() != 0)
        {
            throw std::runtime_error("ERROR: Cannot integrate three spherical harmonics:"
                                     " data structure is not empty on entry to function.");
        }
    );

    using namespace constants;
    data.lmax = lmax;
    Vec1Real fn(4 * lmax + 2);
    math::factorial(fn);
    // Small helper function, if n is even return 1, else -1.
    auto fs = [](const Int n) { return n % 2 == 0 ? 1 : -1; };

    for (Int l1 = 0; l1 <= lmax; ++l1)
    {
        for (Int l2 = 0; l2 <= lmax; ++l2)
        {
            Real k2 = sqrt((2 * l1 + 1) * (2 * l2 + 1));
            for (Int m1 = -l1; m1 <= l1; ++m1)
            {
                for (Int m2 = -l2; m2 <= l2; ++m2)
                {
                    data.indcg.push_back(data.ylm3.size());

                    const Int n1 = std::abs(m1);
                    const Int n2 = std::abs(m2);
                    const Int s1 = m1 < 0 ? 1 : 0;
                    const Int s2 = m2 < 0 ? 1 : 0;
                    const Int t1 = m1 == 0 ? 1 : 0;
                    const Int t2 = m2 == 0 ? 1 : 0;

                    // There are potentially two non-zero m3 values.
                    Int m3;
                    Int m3p;
                    Int nm3;
                    if (m1 * m2 < 0)
                    {
                        m3 = -n1 - n2;
                        m3p = -std::abs(n1 - n2);
                        if (m3p == 0) nm3 = 1;
                        else nm3 = 2;
                    }
                    else if (m1 * m2 == 0)
                    {
                        m3 = m1 + m2;
                        m3p = 0;
                        nm3 = 1;
                    }
                    else
                    {
                        m3 = n1 + n2;
                        m3p = std::abs(n1 - n2);
                        nm3 = 2;
                    }

                    do {
                        const Int n3 = std::abs(m3);
                        const Int s3 = m3 < 0 ? 1 : 0;
                        const Int t3 = m3 == 0 ? 1 : 0;

                        const Real q1 = 0.5 * k2 * fs(n3 + (s1 + s2 + s3) / 2);
                        const Real q2 = 1.0 / std::pow(std::sqrt(2.0), 1 + t1 + t2 + t3);

                        for (Int l3 = std::abs(l1 - l2); l3 <= l1 + l2; l3 += 2)
                        {
                            if (n3 > l3) continue;
                            Real t = 0.0;
                            // clang-format off
                            if (n1 + n2 == -n3) t += clebschGordanM0(l1, l2, l3, fn);
                            if (n1 + n2 ==  n3) t += clebschGordan(l1, l2, l3, n1, n2, n3, fn)
                                                   * fs(n3 + s3);
                            if (n1 - n2 == -n3) t += clebschGordan(l1, l2, l3, n1, -n2, -n3, fn)
                                                   * fs(n2 + s2);
                            if (n1 - n2 ==  n3) t += clebschGordan(l1, l2, l3, -n1, n2, -n3, fn)
                                                   * fs(n1 + s1);
                            // clang-format on
                            const Real t0 = clebschGordanM0(l1, l2, l3, fn);
                            const Real tmp = SQRT_PI * std::sqrt(2 * l3 + 1);
                            data.ylm3.push_back(q1 * q2 * t * t0 / tmp);
                            if (std::abs(t0) < 1.0E-15) data.ylm3i.push_back(0.0);
                            else data.ylm3i.push_back(t * q2 * tmp / (q1 * t0));
                            data.jl.push_back(l3);
                            data.js.push_back(l3 * (l3 + 1) + m3);
                        }
                        nm3--;
                        m3 = m3p;
                    }
                    while (nm3 > 0);
                }
            }
        }
    }

    data.indcg.push_back(data.ylm3.size());
    Int lmind = 0;
    for (Int l1 = 0; l1 <= lmax; ++l1)
    {
        data.lookup.push_back(Vec1Int());
        for (Int l2 = 0; l2 <= lmax; ++l2)
        {
            data.lookup.back().push_back(lmind);
            lmind += (2 * l1 + 1) * (2 * l2 + 1);
        }
    }

    return;
}

Int math::SphericalHarmonics::get_lmax() const
{
    return lmax;
}

Int math::SphericalHarmonics::get_ldim() const
{
    return ldim;
}

Int math::SphericalHarmonics::get_lhc() const
{
    return lhc;
}

const math::SphericalHarmonics::Integrate3SphericalHarmonicsData&
math::SphericalHarmonics::get_i3sh() const
{
    return i3sh;
}
