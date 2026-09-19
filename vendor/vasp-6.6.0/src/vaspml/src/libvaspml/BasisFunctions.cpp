#include "BasisFunctions.hpp"

#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "ParallelEnvironment.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "cutoff.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace vaspml;

using CT = vaspml::math::CutoffType;

BasisFunctions::BasisFunctions() :
    cutoffType(CT::NONE),
    cutoff(0.0),
    widthBroadening(0.0),
    functionType(BasisFunctionType::none)
{}

BasisFunctions::BasisFunctions(const CT                cutoffType_in,
                               const Real              cutoff_in,
                               const Real              widthBroadening_in,
                               const BasisFunctionType type) :
    cutoffType(cutoffType_in),
    cutoff(cutoff_in),
    functionType(type)
{
    set_widthBroadening(widthBroadening_in);
}

void BasisFunctions::printValues()
{
    std::cout << "cutoffType      = " << cutoffType << std::endl;
    std::cout << "cutoff          = " << cutoff << std::endl;
    std::cout << "widthBroadening = " << get_widthBroadening() << std::endl;

    return;
}

// Setters
void BasisFunctions::set_cutoffType(const CT in)
{
    cutoffType = in;
    return;
}

void BasisFunctions::set_cutoff(const Real in)
{
    cutoff = in;

    return;
}

void BasisFunctions::set_widthBroadening(const Real in)
{
    widthBroadening = 0.5 / (in * in);

    return;
}

// Getters
CT BasisFunctions::get_cutoffType() const
{
    return cutoffType;
}

Real BasisFunctions::get_cutoff() const
{
    return cutoff;
}

Real BasisFunctions::get_widthBroadening() const
{
    return std::sqrt(0.5 / widthBroadening);
}

Size BasisFunctions::get_maxOrder() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_maxOrder();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_maxOrder();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }

    return 0;
}

Real BasisFunctions::get_gridSpacing() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_gridSpacing();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_gridSpacing();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }

    return 0;
}

Int BasisFunctions::get_nGrid() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nGrid();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nGrid();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }

    return 0;
}

Int BasisFunctions::get_nRootsBessel() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nRootsBessel();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nRootsBessel();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }

    return 0;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctions::get_totalNumberBasisFunctions() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)
            ->get_totalNumberBasisFunctions();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_totalNumberBasisFunctions();
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif

    return 0;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctions::get_totalNumberBesselFunctions() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)
            ->get_totalNumberBesselFunctions();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_totalNumberBesselFunctions();
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif

    return 0;
}

VASPML_EXEC_SPACE_SPECIFIER
Size BasisFunctions::get_ldim() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_ldim();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_ldim();
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif

    return 0;
}

const Vec1Size& BasisFunctions::get_nValidRoots() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nValidRoots();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nValidRoots();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
        throw std::runtime_error("");
    }
}

void BasisFunctions::set_nGrid(const Size in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_nGrid(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_nGrid(in);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::set_maxOrderBessel(const Size in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_maxOrderBessel(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_maxOrderBessel(in);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::set_nRootsBessel(const Size in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_nRootsBessel(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_nRootsBessel(in);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::printValuesRadial() const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->printValuesRadial();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->printValuesRadial();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::update()
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<BasisFunctionsRadialSpline*>(this)->update();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<BasisFunctionsAngular*>(this)->update();
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

Vec2Real BasisFunctions::interpolate(const Real x) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->interpolate(x);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
        throw std::runtime_error("");
    }
}

void BasisFunctions::interpolate(const Real x, Vec2Real& result) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x, result);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::interpolate(const Real x, Vec2Real& result, Vec2Real& result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result, result_derivative);
    }
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::interpolate(const Real x, Vec1Real& result, Vec1Real& result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result, result_derivative);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::interpolate(const Real  x,
                                 const Size& offset,
                                 Vec1Real&   result,
                                 Vec1Real&   result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          offset,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x,
                                                                     offset,
                                                                     result,
                                                                     result_derivative);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::interpolate(const Real  x,
                                 const Size& order,
                                 const Size& root,
                                 Real&       result,
                                 Real&       result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          order,
                                                                          root,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x,
                                                                     order,
                                                                     root,
                                                                     result,
                                                                     result_derivative);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::computeAngularBasis(const Vec1Real& rhat,
                                         const Vec1Real& rnorm,
                                         Vec1Real&       ylm,
                                         Vec1Real&       ylmd) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->computeAngularBasis(rhat,
                                                                                  rnorm,
                                                                                  ylm,
                                                                                  ylmd);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->computeAngularBasis(rhat,
                                                                             rnorm,
                                                                             ylm,
                                                                             ylmd);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::computeAngularBasis(const Real* rhat,
                                         const Real* rnorm,
                                         const Size  nPair,
                                         Real*       ylm,
                                         Real*       ylmd) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->computeAngularBasis(rhat,
                                                                                  rnorm,
                                                                                  nPair,
                                                                                  ylm,
                                                                                  ylmd);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->computeAngularBasis(rhat,
                                                                             rnorm,
                                                                             nPair,
                                                                             ylm,
                                                                             ylmd);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_EXEC_SPACE_SPECIFIER
void BasisFunctions::computeAngularBasis(const Real& rhat_x,
                                         const Real& rhat_y,
                                         const Real& rhat_z,
                                         const Real& rnorm,
                                         const Size& l,
                                         const Size& m,
                                         Real&       ylm,
                                         Real&       ylmd_x,
                                         Real&       ylmd_y,
                                         Real&       ylmd_z) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)
            ->computeAngularBasis(rhat_x, rhat_y, rhat_z, rnorm, l, m, ylm, ylmd_x, ylmd_y, ylmd_z);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)
            ->computeAngularBasis(rhat_x, rhat_y, rhat_z, rnorm, l, m, ylm, ylmd_x, ylmd_y, ylmd_z);
    }
#ifndef VASPML_PALGO_GPU
    else
    {
        String functionName = __func__;
        global::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}
