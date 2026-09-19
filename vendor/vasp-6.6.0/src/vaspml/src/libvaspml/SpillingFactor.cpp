#include "MlMPI.hpp"

#include "SpillingFactor.hpp"

#include "Linalg.hpp"
#include "Record.hpp"
#include "ShmemArray.hpp"
#include "TypeMap.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace vaspml;

namespace vaspml
{

std::shared_ptr<ShmemArray2DVariableLen<Real>> setup::makeInverseCovMatrixArray(
    const Record&         inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn)
{
    const Vec1Int& numberlocRef = inputParameters.cget<Vec1Int>("numLrc");
    Vec1Int locRefSquared(numberlocRef.size());
    std::transform(numberlocRef.begin(),
                   numberlocRef.end(),
                   locRefSquared.begin(),
                   [](const Int& value) { return value * value; });
    std::shared_ptr<ShmemArray2DVariableLen<Real>> array;
    array = std::make_shared<ShmemArray2DVariableLen<Real>>(locRefSquared, mpiIn);
    const Vec2Real& invKernelMatrix = inputParameters.cget<Vec2Real>("inverse-kernel-matrix");
    for (Size type = 0; type < numberlocRef.size(); type++)
    {
        for (Size locRef = 0; locRef < (Size)locRefSquared[type]; locRef++)
        {
#ifdef use_shmem
            if (!array->get_useMPI()
                || ((*array)["mpiShmem"].get_rank() == 0 && (*array)["mpiInter"].get_rank() == 0))
#endif
                array->set_value(type, locRef, invKernelMatrix[type][locRef]);
        }
    }
#ifdef use_shmem
    (*array)["mpiShmem"].barrier();
    array->distributeInternode();
#endif

    return array;
}

} //namespace vaspml

MapString vaspml::dataMapSpillingFactor()
{
    return MapString{
        {"spillingFactor", "Vec2Real"}
    };
}

SpillingFactor::SpillingFactor(const Record&                 inputParameters,
                               const std::shared_ptr<MlMPI>& mpiIn,
                               const ShRec&                  dataIn) :
    data(assignOrMakeRecord(dataIn, dataMapSpillingFactor())),
    spillingFactor( data->get<Vec2Real>( "spillingFactor" ) ),
    numberLocalRefConfs(std::make_shared<Vec1Int>(inputParameters.cget<Vec1Int>("numLrc"))),
    inverseKernelMatrix(setup::makeInverseCovMatrixArray(inputParameters, mpiIn)),
    isComputed(false),
    statisticsComputed(false),
    averageSpillingFactor(0.0),
    maxSpillingFactor(0.0),
    minSpillingFactor(0.0),
    averageSpillingFactorType(0),
    maxSpillingFactorType(0),
    minSpillingFactorType(0)
{}

void SpillingFactor::computeSpillingFactor(const Vec2Real& kernelMatrix,
                                           const Vec1Int&  nAtomsType,
                                           const TypeMap&  typeMap)
{
    ShmemArray2DVariableLen<Real>& inverseKernelMatrix = (*this->inverseKernelMatrix);
    const Vec1Int&                 numberLocalRefConfs = (*this->numberLocalRefConfs);
    // compute the working array sizes for every atom type of MLFF
    const Vec1Int arraySizes = vector_tools::elementWiseProduct(nAtomsType, numberLocalRefConfs);
    Int           maxSize = vector_tools::get_max(arraySizes);
    Vec1Real      workingArray;
    workingArray.reserve(maxSize);
    vector_tools::allocate_vector(spillingFactor, nAtomsType);
    for (Size type = 0; type < nAtomsType.size(); type++)
    {
        Size typeFF = typeMap.toType(type);
        workingArray.resize(arraySizes[typeFF]);
        linalg::matMul(linalg::Transpose::NOTRANS,
                       linalg::Transpose::NOTRANS,
                       nAtomsType[type],
                       numberLocalRefConfs[typeFF],
                       numberLocalRefConfs[typeFF],
                       linalg::one,
                       kernelMatrix[type],
                       numberLocalRefConfs[typeFF],
                       inverseKernelMatrix.get_slice(typeFF),
                       numberLocalRefConfs[typeFF],
                       (Real)0.0,
                       workingArray,
                       numberLocalRefConfs[typeFF]);
        for (Size atomInType = 0; atomInType < (Size)nAtomsType[type]; atomInType++)
        {
            const Real* workingArrayPtr = &workingArray[atomInType * numberLocalRefConfs[typeFF]];
            const Real* kernelMatrixPtr =
                &kernelMatrix[type][atomInType * numberLocalRefConfs[typeFF]];
            spillingFactor[type][atomInType] =
                std::min(1.0,
                         std::abs(1.0
                                  - linalg::dotProduct(workingArrayPtr,
                                                       kernelMatrixPtr,
                                                       numberLocalRefConfs[typeFF])));
        }
    }
    isComputed = true;
    return;
}

void SpillingFactor::computeStatistics()
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error("ERROR: SpillingFactor::computeStatistics( void )\n"
                                     "Trying to compute statistics of spilling factor, but "
                                     "spilling factor was not computed yet");
        }
    );
    averageSpillingFactorType.resize(spillingFactor.size());
    maxSpillingFactorType.resize(spillingFactor.size());
    minSpillingFactorType.resize(spillingFactor.size());
    for (Size type = 0; type < averageSpillingFactorType.size(); type++)
    {
        averageSpillingFactorType[type] = math::average(spillingFactor[type]);
        minSpillingFactorType[type]     = math::minimum(spillingFactor[type]);
        maxSpillingFactorType[type]     = math::maximum(spillingFactor[type]);
    }
    const Vec1Real sfA = vector_tools::Vec2ToVec1( spillingFactor );
    averageSpillingFactor = math::average(sfA);
    minSpillingFactor     = math::minimum(sfA);
    maxSpillingFactor     = math::maximum(sfA);
    statisticsComputed = true;
    return;
}

const Real& SpillingFactor::get_spillingFactor(const Size type, const Size atom) const
{
    return spillingFactor[type][atom];
}

const Vec1Real& SpillingFactor::get_spillingFactor(const Size type) const
{
    return spillingFactor[type];
}

const Vec2Real& SpillingFactor::get_spillingFactor() const
{
    return spillingFactor;
}

void SpillingFactor::writeToScreen() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: void SpillingFactor::writeToScreen( void ) const.\n"
                "Trying to access spillingFactor array which was not computed yet.\n"
                "Call void SpillingFactor::computeSpillingFactor( const Vec2Real& kernelMatrix, "
                "const Vec1Int& nAtomsType, const TypeMap& typeMap ) first");
        }
    );
    Size atom = 0;
    for (Size type = 0; type < spillingFactor.size(); type++)
    {
        for (Size atomInType = 0; atomInType < spillingFactor[type].size(); atomInType++)
        {
            String output =
                str("%-6zu -%6zu %24.16E ", atom, type, spillingFactor[type][atomInType]);
            std::cout << output << std::endl;
            atom++;
        }
    }
    return;
}

const Real& SpillingFactor::get_averageSpillingFactor() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_averageSpillingFactor( void ) const\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet."
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return averageSpillingFactor;
}

const Real& SpillingFactor::get_maxSpillingFactor() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_maxSpillingFactor( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return maxSpillingFactor;
}

const Real& SpillingFactor::get_minSpillingFactor() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_minSpillingFactor( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return minSpillingFactor;
}

const Vec1Real& SpillingFactor::get_averageSpillingFactorType() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: Vec1Real& SpllingFactor::get_averageSpillingFactorType( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return averageSpillingFactorType;
}

const Vec1Real& SpillingFactor::get_maxSpillingFactorType() const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: Vec1Real& SpllingFactor::get_maxSpillingFactorType( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return maxSpillingFactorType;
}

const Vec1Real& SpillingFactor::get_minSpillingFactorType() const
{
    return minSpillingFactorType;
}

const Real& SpillingFactor::get_averageSpillingFactorType(const Size type) const
{
    return averageSpillingFactorType[type];
}

const Real& SpillingFactor::get_maxSpillingFactorType(const Size type) const
{
    return maxSpillingFactorType[type];
}

const Real& SpillingFactor::get_minSpillingFactorType(const Size type) const
{
    return minSpillingFactorType[type];
}
