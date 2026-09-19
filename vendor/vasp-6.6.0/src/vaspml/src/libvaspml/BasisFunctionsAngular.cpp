#include "BasisFunctionsAngular.hpp"

#include "BasisFunctions.hpp"
#include "ParallelEnvironment.hpp"
#include "SphericalHarmonics.hpp"
#include "constants.hpp"

#include <cmath> // for sqrt

using namespace vaspml;
using CT = vaspml::math::CutoffType;

BasisFunctionsAngular::BasisFunctionsAngular() : BasisFunctionsRadialSpline()
{
    lmax = 0;
}

BasisFunctionsAngular::BasisFunctionsAngular(CT                cutOffType_in,
                                             Real              cutOff_in,
                                             Real              widthBroadening_in,
                                             Size              nGrid_in,
                                             Size              maxOrderBessel_in,
                                             Size              nRootsBessel_in,
                                             BasisFunctionType type,
                                             Size              angularFiltering_in,
                                             Real              filterScale_in) :
    BasisFunctionsRadialSpline(cutOffType_in,
                               cutOff_in,
                               widthBroadening_in,
                               nGrid_in,
                               maxOrderBessel_in,
                               nRootsBessel_in,
                               type)
{
    lmax = maxOrderBessel_in;
    angularFiltering = angularFiltering_in;
    filterScale = filterScale_in;
    angularBasis = std::make_shared<math::SphericalHarmonics>(lmax);
    ldim = angularBasis->get_ldim();
    computeAngularFilteringCoefficient();
    compute_totalNumberBasisFunctions();
}

void BasisFunctionsAngular::set_lmax(Size lmax_in)
{
    lmax = lmax_in;
    angularBasis = std::make_shared<math::SphericalHarmonics>(lmax);
    ldim = angularBasis->get_ldim();

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsAngular::computeAngularBasis(const Vec1Real& normed_xyz,
                                                const Vec1Real& norm,
                                                Vec1Real&       ylm,
                                                Vec1Real&       ylmd) const
{
    angularBasis->computeSphericalHarmonicsAndGradients(normed_xyz, norm, ylm, ylmd);

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsAngular::computeAngularBasis(const Real* normed_xyz,
                                                const Real* norm,
                                                const Size& nPairs,
                                                Real*       ylm,
                                                Real*       ylmd) const
{
    angularBasis->computeSphericalHarmonicsAndGradients(normed_xyz, norm, nPairs, ylm, ylmd);

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsAngular::computeAngularBasis(const Real& normed_x,
                                                const Real& normed_y,
                                                const Real& normed_z,
                                                const Real& norm,
                                                const Size& l,
                                                const Int&  m,
                                                Real&       ylm,
                                                Real&       ylmd_x,
                                                Real&       ylmd_y,
                                                Real&       ylmd_z) const
{
    angularBasis->computeSphericalHarmonicAndGradient(normed_x,
                                                      normed_y,
                                                      normed_z,
                                                      norm,
                                                      l,
                                                      m,
                                                      ylm,
                                                      ylmd_x,
                                                      ylmd_y,
                                                      ylmd_z);

    return;
}

void BasisFunctionsAngular::computeAngularFilteringCoefficient()
{
    Real pi8 = (Real)2 * constants::PI2 * constants::PI2;
    filteringFactors.resize(lmax + 1);
    if (angularFiltering == 1)
    {
        for (Size ll = 0; ll < lmax + 1; ll++)
        {
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1)) / (Real)(2 * ll + 1);
        }
    }
    else if (angularFiltering == 2)
    {
        for (Size ll = 0; ll < lmax + 1; ll++)
        {
            Real factor = (Real)(1 + filterScale * ll * (ll + 1) * ll * (ll + 1));
            factor *= factor;
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1)) / factor;
        }
    }
    else
    {
        for (Size ll = 0; ll < lmax + 1; ll++)
        {
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1));
        }
    }

    return;
}

const Vec1Real& BasisFunctionsAngular::get_filteringFactors() const
{
    return filteringFactors;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctionsAngular::get_ldim() const
{
    return ldim;
}

void BasisFunctionsAngular::compute_totalNumberBasisFunctions()
{
    totalNumberBasisFunctions = 0;
    for (Size ll = 0; ll < lmax + 1; ll++)
    {
        for (Size nRoot = 0; nRoot < nValidRoots[ll]; nRoot++)
        {
            for (Size mm = 0; mm < 2 * ll + 1; mm++) { totalNumberBasisFunctions++; }
        }
    }

    return;
}
