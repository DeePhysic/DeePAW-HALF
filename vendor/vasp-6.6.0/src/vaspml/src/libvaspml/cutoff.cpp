#include "cutoff.hpp"

#include "SmartEnum.hpp"
#include "constants.hpp"
#include "debug.hpp"

#include <algorithm> // for transform, find_if, fill
#include <cmath>     // for cos, exp, sin
#include <stdexcept> // for runtime_error

using namespace vaspml;

template<>
SmartEnum<math::CutoffType>::EnumMap SmartEnum<math::CutoffType>::mapEnumsNames()
{
    return SmartEnum<math::CutoffType>::EnumMap{
        {math::CutoffType::NONE,      "NONE"     },
        {math::CutoffType::HARD,      "HARD"     },
        {math::CutoffType::BP,        "BP"       },
        {math::CutoffType::MIWA,      "MIWA"     },
        {math::CutoffType::JINNOUCHI, "JINNOUCHI"},
        {math::CutoffType::WMC,       "WMC"      },
        {math::CutoffType::POLY2,     "POLY2"    },
        {math::CutoffType::BUMP,      "BUMP"     }
    };
}

void math::cutoff(const Vec1Real& r, Vec1Real& fc, Real rcut, CutoffType cutoffType)
{
    VASPML_DEBUG_L1(
        // Check if output vector has correct size.
        if (r.size() != fc.size())
        {
            throw std::runtime_error("ERROR: Cannot compute cutoff function:"
                                     " output vector has incorrect size.");
        }
        // Sanity-check if any radius is below zero.
        auto rtest = std::find_if(r.begin(), r.end(), [](Real ri) { return ri < 0.0; });
        if (rtest != r.end())
        {
            throw std::runtime_error("ERROR: Radius argument (" + std::to_string(*rtest)
                                     + ") to cutoff function is below zero.");
        }
    );

    if (cutoffType == CutoffType::NONE) { std::fill(fc.begin(), fc.end(), 1.0); }
    else if (cutoffType == CutoffType::HARD)
    {
        std::transform(r.begin(),
                       r.end(),
                       fc.begin(),
                       [&rcut](Real ri)
                       {
                           if (ri < rcut) return 1.0;
                           else return 0.0;
                       });
    }
    else if (cutoffType == CutoffType::BP)
    {
        const Real piOverRcut = constants::PI / rcut;
        std::transform(r.begin(),
                       r.end(),
                       fc.begin(),
                       [&rcut, &piOverRcut](Real ri)
                       {
                           if (ri > rcut) return 0.0;
                           return 0.5 * (std::cos(piOverRcut * ri) + 1.0);
                       });
    }
    else if (cutoffType == CutoffType::MIWA)
    {
        const Real fourOverRcut = 4.0 / rcut;
        std::transform(r.begin(),
                       r.end(),
                       fc.begin(),
                       [&rcut, &fourOverRcut](Real ri)
                       {
                           if (ri > rcut) return 0.0;
                           const Real x = fourOverRcut * ri - 3.0;
                           if (x < -1.0) return 1.0;
                           return x * (0.25 * x * x - 0.75) + 0.5;
                       });
    }
    else if (cutoffType == CutoffType::POLY2)
    {
        std::transform(r.begin(),
                       r.end(),
                       fc.begin(),
                       [&rcut](Real ri)
                       {
                           if (ri > rcut) return 0.0;
                           const Real x = ri / rcut;
                           return ((15.0 - 6.0 * x) * x - 10.0) * x * x * x + 1.0;
                       });
    }
    else if (cutoffType == CutoffType::BUMP)
    {
        std::transform(r.begin(),
                       r.end(),
                       fc.begin(),
                       [&rcut](Real ri)
                       {
                           if (ri > rcut) return 0.0;
                           const Real x = ri / rcut;
                           return constants::E * exp(1.0 / (x * x - 1.0));
                       });
    }
    else if (cutoffType == CutoffType::JINNOUCHI || cutoffType == CutoffType::WMC)
    {
        throw std::runtime_error("ERROR: Cutoff function type \"" + toString(cutoffType)
                                 + "\" not implemented yet.");
    }

    return;
}

void math::cutoffAndDerivative(const Vec1Real& r,
                               Vec1Real&       fc,
                               Vec1Real&       dfc,
                               Real            rcut,
                               CutoffType      cutoffType)
{
    VASPML_DEBUG_L1(
        // Check if output vectors have correct size.
        if (r.size() != fc.size() || r.size() != dfc.size())
        {
            throw std::runtime_error("ERROR: Cannot compute cutoff function and derivatives:"
                                     " output vectors have incorrect sizes.");
        }
        // Sanity-check if any radius is below zero.
        auto rtest = std::find_if(r.begin(), r.end(), [](Real ri) { return ri < 0.0; });
        if (rtest != r.end())
        {
            throw std::runtime_error("ERROR: Radius argument (" + std::to_string(*rtest)
                                     + ") to cutoff function is below zero.");
        }
    );

    if (cutoffType == CutoffType::NONE)
    {
        for (Size i = 0; i < r.size(); ++i)
        {
            fc[i] = 1.0;
            dfc[i] = 0.0;
        }
    }
    else if (cutoffType == CutoffType::HARD)
    {
        for (Size i = 0; i < r.size(); ++i)
        {
            if (r[i] < rcut) fc[i] = 1.0;
            else fc[i] = 0.0;
            dfc[i] = 0.0;
        }
    }
    else if (cutoffType == CutoffType::BP)
    {
        const Real w = constants::PI / rcut;
        for (Size i = 0; i < r.size(); ++i)
        {
            if (r[i] > rcut)
            {
                fc[i] = 0.0;
                dfc[i] = 0.0;
            }
            else
            {
                const Real temp = w * r[i];
                fc[i] = 0.5 * (std::cos(temp) + 1.0);
                dfc[i] = -0.5 * w * std::sin(temp);
            }
        }
    }
    else if (cutoffType == CutoffType::MIWA)
    {
        const Real fourOverRcut = 4.0 / rcut;
        for (Size i = 0; i < r.size(); ++i)
        {
            if (r[i] > rcut)
            {
                fc[i] = 0.0;
                dfc[i] = 0.0;
            }
            else
            {
                const Real x = fourOverRcut * r[i] - 3.0;
                if (x < -1.0)
                {
                    fc[i] = 1.0;
                    dfc[i] = 0.0;
                }
                else
                {
                    fc[i] = x * (0.25 * x * x - 0.75) + 0.5;
                    dfc[i] = fourOverRcut * 0.75 * (x * x - 1.0);
                }
            }
        }
    }
    else if (cutoffType == CutoffType::POLY2)
    {
        const Real oneOverRcut = 1.0 / rcut;
        for (Size i = 0; i < r.size(); ++i)
        {
            if (r[i] > rcut)
            {
                fc[i] = 0.0;
                dfc[i] = 0.0;
            }
            else
            {
                const Real x = r[i] * oneOverRcut;
                const Real x2 = x * x;
                fc[i] = ((15.0 - 6.0 * x) * x - 10.0) * x * x2 + 1.0;
                dfc[i] = oneOverRcut * x2 * ((60.0 - 30.0 * x) * x - 30.0);
            }
        }
    }
    else if (cutoffType == CutoffType::BUMP)
    {
        const Real oneOverRcut = 1.0 / rcut;
        for (Size i = 0; i < r.size(); ++i)
        {
            if (r[i] > rcut)
            {
                fc[i] = 0.0;
                dfc[i] = 0.0;
            }
            else
            {
                const Real x = r[i] * oneOverRcut;
                const Real temp = 1.0 / (x * x - 1.0);
                const Real temp2 = exp(temp);
                fc[i] = constants::E * temp2;
                dfc[i] = -2.0 * oneOverRcut * constants::E * x * temp * temp * temp2;
            }
        }
    }
    else if (cutoffType == CutoffType::JINNOUCHI || cutoffType == CutoffType::WMC)
    {
        throw std::runtime_error("ERROR: Cutoff function type \"" + toString(cutoffType)
                                 + "\" not implemented yet.");
    }

    return;
}
