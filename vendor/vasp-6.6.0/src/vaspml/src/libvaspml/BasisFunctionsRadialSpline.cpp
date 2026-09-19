#include "BasisFunctionsRadialSpline.hpp"

#include "ParallelEnvironment.hpp"
#include "constants.hpp"
#include "cutoff.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm> // for fill, min
#include <cmath>     // for sqrt, pow
#include <iostream>  // for basic_ostream, char_traits, opera...
#include <stdexcept> // for runtime_error

using namespace vaspml;
using CT = vaspml::math::CutoffType;

BasisFunctionsRadialSpline::BasisFunctionsRadialSpline() : BasisFunctions()
{
    maxOrderBessel = 0;
    nRootsBessel = 0;
    nGrid = 0;
}

BasisFunctionsRadialSpline::BasisFunctionsRadialSpline(CT                cutoffType_in,
                                                       Real              cutoff_in,
                                                       Real              widthBroadening_in,
                                                       Size              nGrid_in,
                                                       Size              maxOrderBessel_in,
                                                       Size              nRootsBessel_in,
                                                       BasisFunctionType type) :
    BasisFunctions(cutoffType_in, cutoff_in, widthBroadening_in, type)
{
    nGrid = nGrid_in;
    gridSpacing = (cutoff) / (Real)nGrid;
    maxOrderBessel_store = maxOrderBessel_in;
    maxOrderBessel = maxOrderBessel_store + 1;
    nRootsBessel = nRootsBessel_in;
    update();
}

Real BasisFunctionsRadialSpline::get_gridSpacing() const
{
    return gridSpacing;
}

Int BasisFunctionsRadialSpline::get_nGrid() const
{
    return nGrid;
}

Size BasisFunctionsRadialSpline::get_maxOrder() const
{
    return maxOrderBessel_store;
}

Int BasisFunctionsRadialSpline::get_nRootsBessel() const
{
    return nRootsBessel;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctionsRadialSpline::get_totalNumberBasisFunctions() const
{
    return totalNumberBasisFunctions;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctionsRadialSpline::get_totalNumberBesselFunctions() const
{
    return totalNumberBesselFunctions;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctionsRadialSpline::get_ldim() const
{
    return 1;
}

const Vec1Size& BasisFunctionsRadialSpline::get_nValidRoots() const
{
    return nValidRoots;
}

void BasisFunctionsRadialSpline::set_nGrid(const Size in)
{
    nGrid = in;

    return;
}

void BasisFunctionsRadialSpline::set_maxOrderBessel(const Size in)
{
    maxOrderBessel = in + 1;
    maxOrderBessel_store = in;

    return;
}

void BasisFunctionsRadialSpline::set_nRootsBessel(const Size in)
{
    nRootsBessel = in;

    return;
}

void BasisFunctionsRadialSpline::set_gridSpacing()
{
    gridSpacing = cutoff / (Real)nGrid;

    return;
}

void BasisFunctionsRadialSpline::makeGrid()
{
    grid.resize(nGrid);
    for (Size i = 0; i < nGrid; i++) { grid[i] = gridSpacing * (i + 1.0); }

    return;
}

void BasisFunctionsRadialSpline::printValuesRadial() const
{
    std::cout << "nGrid          = " << nGrid << std::endl;
    std::cout << "maxOrderBessel = " << maxOrderBessel_store << std::endl;
    std::cout << "nRootsBessel   = " << nRootsBessel << std::endl;

    return;
}

void BasisFunctionsRadialSpline::computeCutoffFunction()
{
    cutoffFunction.resize(nGrid);
    cutoffFunction_derivative.resize(nGrid);

    math::cutoffAndDerivative(grid, cutoffFunction, cutoffFunction_derivative, cutoff, cutoffType);

    return;
}

void BasisFunctionsRadialSpline::update()
{
    Vec3Real sphericalBessel;
    set_gridSpacing();
    makeGrid();
    // compute cutoff function for radial basis
    computeCutoffFunction();

    // temporary array to store roots (zeros) of the basis set
    Vec2Real rootsBesselBasis;
    // compute roots of Bessel functions
    computeSphericalBesselRoots(rootsBesselBasis, maxOrderBessel);
    checkSphericalBesselRootsConvergence(rootsBesselBasis);
    // rescaling roots such that Basis set is zero at cutoff
    rescaleRoots(rootsBesselBasis);
    // compute the spherical Bessel functions, scaled by the roots
    computeSphericalBesselFunctions(sphericalBessel, rootsBesselBasis);
    // normalize spherical Bessel functions
    computeBesselNorms(sphericalBessel);

    Vec3Real projectionRadial;
    Vec2Real projectionRadialDerivativeOrigin;
    computeProjectionRadial(sphericalBessel, projectionRadial);
    // compute derivatives of modified spherical bessel functions
    // at grid origin to compute spline interpolation
    computeProjectionRadialDerivativeOrigin(sphericalBessel, projectionRadialDerivativeOrigin);
    computeSplineCoefficients(projectionRadial, projectionRadialDerivativeOrigin);

    return;
}

void BasisFunctionsRadialSpline::checkSphericalBesselRootsConvergence(Vec2Real& roots)
{
    nValidRoots.resize(maxOrderBessel);
    std::fill(nValidRoots.begin(), nValidRoots.end(), 0);
    totalNumberBasisFunctions = 0;
    totalNumberBesselFunctions = 0;
    for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        for (Size i = 0; i < roots[actOrder].size(); i++)
        {
            if (roots[actOrder][i] <= rootMax)
            {
                nValidRoots[actOrder]++;
                totalNumberBasisFunctions++;
                totalNumberBesselFunctions++;
            }
            else { roots[actOrder][i] = (Real)0; }
        }
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselRoots(Vec2Real& rootsBesselBasis,
                                                             const Int maxOrderBessel)
{
    vector_tools::allocate_vector(rootsBesselBasis, maxOrderBessel, nRootsBessel);
    rootMax = 1.0e12;
    for (Int actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        math::sphericalBesselRoots(rootsBesselBasis[actOrder], actOrder);
        rootMax = std::min(rootMax, rootsBesselBasis[actOrder][nRootsBessel - 1]);
    }

    return;
}

void BasisFunctionsRadialSpline::rescaleRoots(Vec2Real& rootsBesselBasis)
{
    for (auto& r : rootsBesselBasis)
    {
        for (auto& x : r) { x /= cutoff; }
    }
    rootMax /= cutoff;

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctions(Vec3Real&       bessel,
                                                                 const Vec2Real& rootsBesselBasis)
{
    vector_tools::allocate_vector(bessel, maxOrderBessel);
    for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        computeSphericalBesselFunctionsSingleOrder(bessel[actOrder],
                                                   rootsBesselBasis[actOrder],
                                                   actOrder,
                                                   nValidRoots[actOrder]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctionsSingleOrder(
    Vec2Real&       bessel,
    const Vec1Real& roots,
    const Size      actOrderBessel,
    const Size      nBasis)
{
    vector_tools::allocate_vector(bessel, nBasis, nGrid);
    for (Size i = 0; i < nBasis; i++)
    {
        computeSphericalBesselFunctionsSingleRoot(bessel[i], roots[i], actOrderBessel);
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctionsSingleRoot(Vec1Real& bessel,
                                                                           Real      root,
                                                                           const Int actOrderBessel)
{
    for (Size i = 0; i < bessel.size(); i++)
    {
        bessel[i] = math::sphericalBessel(grid[i] * root, actOrderBessel);
    }

    return;
}

void BasisFunctionsRadialSpline::computeBesselNorms(Vec3Real& bessel)
{

    for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        computeBesselNormsSingleOrder(bessel[actOrder], nValidRoots[actOrder]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeBesselNormsSingleOrder(Vec2Real& bessel, Size nBasis)
{
    Vec1Real normalization;
    normalization.resize(nBasis);
    for (Size i = 0; i < nBasis; i++)
    {
        // Originally, in VASPml this integration scheme was used:
        //normalization[i] = math::trapezSphericalL2Integral(bessel[i], grid, gridSpacing);
        // However, for compatibility reasons we use this simpler form instead:
        for (Size j = 0; j < nGrid; j++)
        {
            const double tmp = bessel[i][j] * grid[j];
            normalization[i] += tmp * tmp * gridSpacing;
        }
        normalization[i] = std::sqrt(normalization[i]);
        math::vectorTimesScalarNoCopy(bessel[i], (Real)1 / normalization[i]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeProjectionRadial(const Vec3Real& sphericalBessel,
                                                         Vec3Real&       projectionRadial)
{
    Real prefactor = gridSpacing * constants::PI4 * pow(widthBroadening / constants::PI, (Real)1.5);

    Vec1Real rrj = grid;
    math::squareVector(rrj);
    math::vectorTimesScalarNoCopy(rrj, widthBroadening);

    Vec2Real product;
    vector_tools::allocate_vector(product, maxOrderBessel, nGrid);

    Vec1Real temp_integral(nGrid);

    // allocate modifiedSphericalBessel maxOrderBessel,nRoots[l],nGrid
    vector_tools::allocate_vector(projectionRadial, maxOrderBessel);
    for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        vector_tools::allocate_vector(projectionRadial[actOrder], nValidRoots[actOrder]);
        for (Size root = 0; root < nValidRoots[actOrder]; root++)
        {
            vector_tools::allocate_vector(projectionRadial[actOrder][root], nGrid);
        }
    }

    for (Size i = 0; i < nGrid; i++)
    {
        Vec1Real rr = math::vectorTimesScalar(grid, 2.0 * widthBroadening * grid[i]);
        Real     rri = widthBroadening * grid[i] * grid[i];
        math::expModSphBessel(rri, rr, rrj, product);

        for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
        {
            for (Size j = 0; j < nGrid; j++)
            {
                temp_integral[j] =
                    prefactor * cutoffFunction[i] * product[actOrder][j] * grid[j] * grid[j];
            }
            for (Size root = 0; root < nValidRoots[actOrder]; root++)
            {
                projectionRadial[actOrder][root][i] =
                    math::dotProduct(temp_integral, sphericalBessel[actOrder][root]);
            }
        }
    }

    return;
}

void BasisFunctionsRadialSpline::computeProjectionRadialDerivativeOrigin(
    const Vec3Real& sphericalBessel,
    Vec2Real&       projectionRadialDerivativeOrigin)
{
    // For more details on this routine one should consider the paper
    // https://doi.org/10.1103/PhysRevB.100.014105, equation 19.

    Real prefactor = gridSpacing * constants::PI4 * pow(widthBroadening / constants::PI, 1.5);

    Vec1Real gauss;
    Vec1Real gauss_derivative;

    vector_tools::allocate_vector(gauss, nGrid);
    vector_tools::allocate_vector(gauss_derivative, nGrid);

    Vec2Real modSphBessel;
    Vec2Real modSphBesselDerivative;
    vector_tools::allocate_vector(modSphBessel, maxOrderBessel, nGrid);
    vector_tools::allocate_vector(modSphBesselDerivative, maxOrderBessel, nGrid);

    math::gaussianAndDerivative(grid, gauss, gauss_derivative, widthBroadening);

    // Rescaled grid used to compute the spherical Bessel function and its derivative.
    Vec1Real rescaledGrid = math::vectorTimesScalar(grid, (Real)2 * widthBroadening * grid[0]);

    math::modifiedSphericalBesselAndDerivative(rescaledGrid, modSphBessel, modSphBesselDerivative);

    // Precomputed factors needed for computing the integral.
    Real factor1 = prefactor * cutoffFunction_derivative[0] * gauss[0];
    Real factor2 = prefactor * cutoffFunction[0] * gauss_derivative[0];
    Real factor3 = prefactor * cutoffFunction[0] * gauss[0] * (Real)2 * widthBroadening;
    vector_tools::allocate_vector(projectionRadialDerivativeOrigin, maxOrderBessel);
    for (Size actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        projectionRadialDerivativeOrigin[actOrder].resize(nValidRoots[actOrder]);
        for (Size root = 0; root < nValidRoots[actOrder]; root++)
        {
            Real value = computeProjectionRadialDerivativeOriginSingle(
                factor1,
                factor2,
                factor3,
                gauss,
                sphericalBessel[actOrder][root],
                modSphBessel[actOrder].begin(),
                modSphBesselDerivative[actOrder].begin());
            projectionRadialDerivativeOrigin[actOrder][root] = value;
        }
    }

    return;
}

Real BasisFunctionsRadialSpline::computeProjectionRadialDerivativeOriginSingle(
    const Real               factor1,
    const Real               factor2,
    const Real               factor3,
    const Vec1Real&          gauss,
    const Vec1Real&          spherical_bessel,
    Vec1Real::const_iterator modifiedBessel,
    Vec1Real::const_iterator modifiedBesselDerivative)
{
    Real integral = 0;
    for (Size ir = 0; ir < nGrid; ir++)
    {
        integral += (gauss[ir] * (*(modifiedBessel + ir)) * (factor1 + factor2)
                     + factor3 * gauss[ir] * (*(modifiedBesselDerivative + ir)) * grid[ir])
                  * grid[ir] * grid[ir] * spherical_bessel[ir];
    }

    return integral;
}

void BasisFunctionsRadialSpline::computeSplineCoefficients(const Vec3Real& functions,
                                                           const Vec2Real& derivative)
{
    vector_tools::allocate_vector(splines, functions.size());
    for (Size actOrder = 0; actOrder < functions.size(); actOrder++)
    {
        vector_tools::allocate_vector(splines[actOrder], nValidRoots[actOrder]);
        for (Size root = 0; root < functions[actOrder].size(); root++)
        {
            splines[actOrder][root].setup(grid,
                                          functions[actOrder][root],
                                          derivative[actOrder][root]);
        }
    }

    return;
}

Vec2Real BasisFunctionsRadialSpline::interpolate(const Real x) const
{
    Vec2Real result;
    vector_tools::allocate_vector(result, splines.size());

    for (Size actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        vector_tools::allocate_vector(result[actOrder], nValidRoots[actOrder]);
        for (Size root = 0; root < splines[actOrder].size(); root++)
        {
            result[actOrder][root] = splines[actOrder][root].interpolate(x);
        }
    }

    return result;
}

void BasisFunctionsRadialSpline::interpolate(const Real x, Vec2Real& values) const
{
    for (Size actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        for (Size root = 0; root < splines[actOrder].size(); root++)
        {
            values[actOrder][root] = splines[actOrder][root].interpolate(x);
        }
    }

    return;
}

void BasisFunctionsRadialSpline::interpolate(const Real x,
                                             Vec2Real&  values,
                                             Vec2Real&  derivatives) const
{
    vector_tools::allocate_vector(values, splines.size());
    vector_tools::allocate_vector(derivatives, splines.size());
    for (Size actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        values[actOrder].resize(splines[actOrder].size());
        derivatives[actOrder].resize(splines[actOrder].size());
        for (Size root = 0; root < splines[actOrder].size(); root++)
        {
            splines[actOrder][root].interpolateWithDerivative(x,
                                                              values[actOrder][root],
                                                              derivatives[actOrder][root]);
        }
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::interpolate(const Real x,
                                             Vec1Real&  values,
                                             Vec1Real&  derivatives) const
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (totalNumberBesselFunctions != values.size())
        {
            throw std::runtime_error("ERROR in BasisFunctionsRadialSpline::interpolate: Size of "
                                     "supplied function array incorrect");
        }
    );
    VASPML_DEBUG_L1(
        if (totalNumberBesselFunctions != derivatives.size())
        {
            throw std::runtime_error(
                "ERROR in BasisFunctionsRadialSpline::interpolate: Size of supplied derivative"
                "function array incorrect");
        }
    );
#endif

    Size col = 0;
    Real func;
    Real func_deri;
    for (Size actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        for (Size root = 0; root < splines[actOrder].size(); root++)
        {
            splines[actOrder][root].interpolateWithDerivative(x, func, func_deri);
            values[col] = func;
            derivatives[col] = func_deri;
            col++;
        }
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::interpolate(const Real  x,
                                             const Size& offset,
                                             Vec1Real&   values,
                                             Vec1Real&   derivatives) const
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (totalNumberBesselFunctions != values.size())
        {
            throw std::runtime_error("ERROR in BasisFunctionsRadialSpline::interpolate: Size of "
                                     "supplied function array incorrect");
        }
    );
    VASPML_DEBUG_L1(
        if (totalNumberBesselFunctions != derivatives.size())
        {
            throw std::runtime_error(
                "ERROR in BasisFunctionsRadialSpline::interpolate: Size of supplied derivative"
                "function array incorrect");
        }
    );
#endif

    Size col = offset;
    for (Size actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        for (Size root = 0; root < splines[actOrder].size(); root++)
        {
            splines[actOrder][root].interpolateWithDerivative(x, values[col], derivatives[col]);
            col++;
        }
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::interpolate(const Real& x,
                                             const Size& order,
                                             const Size& root,
                                             Real&       values,
                                             Real&       derivatives) const
{
    splines[order][root].interpolateWithDerivative(x, values, derivatives);
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::computeAngularBasis(const Vec1Real& /* rhat */,
                                                     const Vec1Real& /* rnorm */,
                                                     Vec1Real& ylm,
                                                     Vec1Real& ylmd) const
{
    const Real radFactor = (Real)1 / std::sqrt((Real)4 * constants::PI);
    std::fill(ylm.begin(), ylm.end(), radFactor);
    std::fill(ylmd.begin(), ylmd.end(), (Real)0);

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::computeAngularBasis(const Real* /* rhat */,
                                                     const Real* /* rnorm */,
                                                     const Size nPair,
                                                     Real*      ylm,
                                                     Real*      ylmd) const
{
    const Real radFactor = (Real)1 / std::sqrt((Real)4 * constants::PI);
    for (Size i = 0; i < nPair * totalNumberBasisFunctions; i++)
    {
        ylm[i] = radFactor;
        ylmd[3 * i] = 0.0;
        ylmd[3 * i + 1] = 0.0;
        ylmd[3 * i + 2] = 0.0;
    }

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctionsRadialSpline::computeAngularBasis(const Real&, //rhat_x,
                                                     const Real&, //rhat_y,
                                                     const Real&, //rhat_z,
                                                     const Real&, //rnorm,
                                                     const Size&, //l,
                                                     const Size&, //m,
                                                     Real& ylm,
                                                     Real& ylmd_x,
                                                     Real& ylmd_y,
                                                     Real& ylmd_z) const
{
    const Real radFactor = (Real)1 / std::sqrt((Real)4 * constants::PI);
    ylm = radFactor;
    ylmd_x = 0.0;
    ylmd_y = 0.0;
    ylmd_z = 0.0;

    return;
}
