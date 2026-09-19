#include "MlMPI.hpp"

#include "Predictor.hpp"

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "Kernel.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "ShmemArray.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm> // for fill, for_each
#include <fstream>   // for basic_ostream, basic_ofstream, basic_ios, endl
#include <iostream>  // for cout
#include <numeric>   // for iota
#include <stdexcept> // for runtime_error

using namespace vaspml;

template<>
SmartEnum<TotalEnergyType>::EnumMap SmartEnum<TotalEnergyType>::mapEnumsNames()
{
    return SmartEnum<TotalEnergyType>::EnumMap{
        {TotalEnergyType::IsolatedAtom,       "IsolatedAtom"         },
        {TotalEnergyType::AverageTrainEnergy, "AverageTrainingEnergy"}
    };
}

MapString vaspml::dataMapPredictor()
{
    return MapString{
        {"energyArray",       "Vec2Real"},
        {"atomicForces",      "Vec1Real"},
        {"totalStressTensor", "Vec1Real"}
    };
}

namespace vaspml
{

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> make_descriptorsRefConfsWeighted(
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::map<String, Real>&      weights,
    const std::shared_ptr<MlMPI>&      mpiIn)
{
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> arrays;
    for (const auto& [key, item] : featureSpaceSize)
    {
        arrays[key] = std::make_shared<ShmemArray2DVariableLen<Real>>(
            weights.at(key) > 0 ? *featureSpaceSize.at(key) : Vec1Int{0},
            mpiIn);
    }

    return arrays;
}

std::map<String, ShVec2Real> makeMaps2D()
{

    std::map<String, ShVec2Real> data;
    for (const String& key : constants::descriptorKeyList)
    {
        data[key] = std::make_shared<Vec2Real>();
    }

    return data;
}

std::map<String, ShVec1Real> makeMaps1D()
{
    std::map<String, ShVec1Real> data;
    for (const String& key : constants::descriptorKeyList)
    {
        data[key] = std::make_shared<Vec1Real>();
    }

    return data;
}

std::map<String, Vec1Size> makeMaps1DSize()
{
    std::map<String, Vec1Size> data;
    for (const String& key : constants::descriptorKeyList) data[key] = Vec1Size();

    return data;
}

std::map<String, ArrayResizing1D> makeResizeMap1D()
{
    std::map<String, ArrayResizing1D> data;
    for (const String& key : constants::descriptorKeyList) data[key] = ArrayResizing1D();

    return data;
}

std::map<String, ArrayResizing2D> makeResizeMap2D()
{
    std::map<String, ArrayResizing2D> data;
    for (const String& key : constants::descriptorKeyList) data[key] = ArrayResizing2D();

    return data;
}

} //namespace vaspml

Predictor::Predictor(ShRec record) :
    data(assignOrMakeRecord(record, dataMapPredictor())),
    energyArray(data->get<Vec2Real>("energyArray")),
    totalStressTensor(data->get<Vec1Real>("totalStressTensor")),
    atomicForces(data->get<Vec1Real>("atomicForces"))
{}

Predictor::Predictor(const Record&                 inputParameters,
                     const std::shared_ptr<MlMPI>& mpiIn,
                     ShRec                         record) :
    data(assignOrMakeRecord(record, dataMapPredictor())),
    numberTypesForceField(inputParameters.cget<Int>("numTypes")),
    numberLocalRefConfs(std::make_shared<Vec1Int>(inputParameters.cget<Vec1Int>("numLrc"))),
    numberDescriptors(setup::make_numberDescriptorsMap(inputParameters)),
    featureSpaceSize(setup::generate_featureSpaceMap(inputParameters)),
    weights(setup::make_weightsMap(inputParameters)),
    regressionCoefficients(setup::make_regressionCoefficients(inputParameters, mpiIn)),
    referenceEnergyPerType(
        std::make_shared<Vec1Real>(inputParameters.cget<Vec1Real>("ML_EATOM_REF"))),
    averageTrainEnergy(inputParameters.cget<Real>("average-energy-per-atom")),
    descriptorsRefConfsWeighted(make_descriptorsRefConfsWeighted(featureSpaceSize, weights, mpiIn)),
    energyArray(data->get<Vec2Real>("energyArray")),
    derivativeMatrix(makeMaps2D()),
    derivativeMatrixDim(makeMaps1DSize()),
    derivativeMatrixSize(makeResizeMap2D()),
    forcePreContract(makeMaps2D()),
    forcePreContractDim(makeMaps1DSize()),
    forcePreContractSize(makeResizeMap2D()),
    pairForces(makeMaps2D()),
    pairForcesDim(makeMaps1DSize()),
    pairForcesSize(makeResizeMap2D()),
    centralForces(makeMaps1D()),
    centralForcesSize(makeResizeMap1D()),
    stressTensor(makeMaps1D()),
    totalStressTensor(data->get<Vec1Real>("totalStressTensor")),
    atomicForces(data->get<Vec1Real>("atomicForces"))
{
    forceArraysComputed = false;

    compute_descriptorsRefConfsWeighted(inputParameters);

    switch (inputParameters.cget<Int>("ML_ISCALE_TOTEN"))
    {
    case 1:
        this->energyType = TotalEnergyType::IsolatedAtom;
        break;
    case 2:
        this->energyType = TotalEnergyType::AverageTrainEnergy;
        break;
    }
}

void Predictor::init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn)
{
    numberTypesForceField = inputParameters.cget<Int>("numTypes");
    numberLocalRefConfs = std::make_shared<Vec1Int>(inputParameters.cget<Vec1Int>("numLrc"));
    numberDescriptors = setup::make_numberDescriptorsMap(inputParameters);
    featureSpaceSize = setup::generate_featureSpaceMap(inputParameters);
    weights = setup::make_weightsMap(inputParameters);
    regressionCoefficients = setup::make_regressionCoefficients(inputParameters, mpiIn);
    referenceEnergyPerType =
        std::make_shared<Vec1Real>(inputParameters.cget<Vec1Real>("ML_EATOM_REF"));
    averageTrainEnergy = inputParameters.cget<Real>("average-energy-per-atom");
    descriptorsRefConfsWeighted =
        make_descriptorsRefConfsWeighted(featureSpaceSize, weights, mpiIn);
    derivativeMatrix = makeMaps2D();
    derivativeMatrixDim = makeMaps1DSize();
    derivativeMatrixSize = makeResizeMap2D();
    forcePreContract = makeMaps2D();
    forcePreContractDim = makeMaps1DSize();
    forcePreContractSize = makeResizeMap2D();
    pairForces = makeMaps2D();
    pairForcesDim = makeMaps1DSize();
    pairForcesSize = makeResizeMap2D();
    centralForces = makeMaps1D();
    centralForcesSize = makeResizeMap1D();
    stressTensor = makeMaps1D();

    forceArraysComputed = false;

    compute_descriptorsRefConfsWeighted(inputParameters);

    switch (inputParameters.cget<Int>("ML_ISCALE_TOTEN"))
    {
    case 1:
        this->energyType = TotalEnergyType::IsolatedAtom;
        break;
    case 2:
        this->energyType = TotalEnergyType::AverageTrainEnergy;
        break;
    }

    return;
}

void Predictor::compute_descriptorsRefConfsWeighted(const Record& inputParameters)
{
    for (String key : constants::descriptorKeyList)
    {
        String          tag = key + "-reference-configs";
        const Vec2Real& regCoeff = inputParameters.cget<Vec2Real>("regression-coeff");
        if (weights[key] > 0)
        {
            const Vec2Real descriptor = inputParameters.cget<Vec2Real>(tag);
            for (Size type = 0; type < descriptor.size(); type++)
            {
                Size featureSpaceCounter = 0;
                for (Size local_reference = 0; local_reference < (Size)(*numberLocalRefConfs)[type];
                     local_reference++)
                {
                    for (Size desc_count = 0; desc_count < (Size)(*numberDescriptors[key])[type];
                         desc_count++)
                    {
                        if (!descriptorsRefConfsWeighted[key]->get_useMPI()
                            || ((*descriptorsRefConfsWeighted[key])["mpiShmem"].get_rank() == 0
                                && (*descriptorsRefConfsWeighted[key])["mpiInter"].get_rank() == 0))
                        {
                            descriptorsRefConfsWeighted[key]->set_value(
                                type,
                                featureSpaceCounter,
                                regCoeff[type][local_reference]
                                    * descriptor[type][featureSpaceCounter]);
                            featureSpaceCounter++;
                        }
                    }
                }
            }
            if (descriptorsRefConfsWeighted[key]->get_useMPI())
            {
                (*descriptorsRefConfsWeighted[key])["mpiShmem"].barrier();
                descriptorsRefConfsWeighted[key]->distributeInternode();
            }
        }
    }

    return;
}

void Predictor::compute_energyArrayDim(const Size numberTypesStruc)
{
    energyArrayDim.resize(numberTypesStruc);
    for (Size typeStruc = 0; typeStruc < numberTypesStruc; typeStruc++)
    {
        Size typeFF = typeMap->toType(typeStruc);
        energyArrayDim[typeStruc] = (*numberLocalRefConfs)[typeFF];
    }

    return;
}

void Predictor::allocate_energyArray(const Size numberTypesStruc)
{
    compute_energyArrayDim(numberTypesStruc);
    if (global::parallel.off()) allocate_energyArrayCPU(numberTypesStruc);
    else allocate_energyArrayGPU(numberTypesStruc);

    return;
}

void Predictor::allocate_energyArrayCPU(const Size numberTypesStruc)
{
    energyArray.resize(numberTypesStruc);
    for (Size typeStruc = 0; typeStruc < numberTypesStruc; typeStruc++)
    {
        energyArray[typeStruc].resize(energyArrayDim[typeStruc]);
        std::fill(energyArray[typeStruc].begin(), energyArray[typeStruc].end(), (Real)0);
    }

    return;
}

void Predictor::compute_derivativeMatrixDim(const Vec1Int& numberAtomType)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        derivativeMatrixDim[key].resize(numberAtomType.size());
        for (Size typeStruc = 0; typeStruc < numberAtomType.size(); typeStruc++)
        {
            Size typeForceField = typeMap->toType(typeStruc);
            derivativeMatrixDim[key][typeStruc] =
                numberAtomType[typeStruc] * (*numberDescriptors[key])[typeForceField];
        }
    }

    return;
}

void Predictor::allocate_derivativeMatrix(const Vec1Int& numberAtomType)
{
    compute_derivativeMatrixDim(numberAtomType);
    if (global::parallel.off()) allocate_derivativeMatrixCPU(numberAtomType);
    else allocate_derivativeMatrixGPU(numberAtomType);

    return;
}

void Predictor::allocate_derivativeMatrixCPU(const Vec1Int& numberAtomType)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        derivativeMatrix[key]->resize(numberAtomType.size());
        for (Size typeStruc = 0; typeStruc < numberAtomType.size(); typeStruc++)
        {
            //Size typeForceField  =  typeMap->toType( typeStruc );
            (*derivativeMatrix[key])[typeStruc].resize(derivativeMatrixDim[key][typeStruc]);
            std::fill((*derivativeMatrix[key])[typeStruc].begin(),
                      (*derivativeMatrix[key])[typeStruc].end(),
                      (Real)0);
        }
    }

    return;
}

void Predictor::compute_forcePreContractDim(const DescriptorCollector& descriptorCollection,
                                            const Size                 numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        forcePreContractDim[key].resize(numberAtoms);
        for (Size nAtom = 0; nAtom < numberAtoms; nAtom++)
        {
            const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
            forcePreContractDim[key][nAtom] = descriptor.get_forcePreContractSize(nAtom);
        }
    }

    return;
}

void Predictor::allocate_forcePreContract(const DescriptorCollector& descriptorCollection,
                                          const Size                 numberAtoms)
{
    compute_forcePreContractDim(descriptorCollection, numberAtoms);
    if (global::parallel.off()) allocate_forcePreContractCPU(numberAtoms);
    else allocate_forcePreContractGPU(numberAtoms);

    return;
}

void Predictor::allocate_forcePreContractCPU(const Size numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        forcePreContract[key]->resize(numberAtoms);
        for (Size nAtom = 0; nAtom < numberAtoms; nAtom++)
        {
            (*forcePreContract[key])[nAtom].resize(forcePreContractDim[key][nAtom]);
            std::fill((*forcePreContract[key])[nAtom].begin(),
                      (*forcePreContract[key])[nAtom].end(),
                      (Real)0);
        }
    }

    return;
}

void Predictor::allocate_centralAtomIndex(const Size numberAtoms)
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }

    return;
}

void Predictor::compute_pairForcesDim(const DescriptorCollector& descriptorCollection)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        std::shared_ptr<const NearestNeighborNSquare> neighborList =
            descriptor.get_neighborList_ptr();
        pairForcesDim[key].resize(neighborList->get_nAtoms());
        for (Size atom = 0; atom < pairForcesDim[key].size(); atom++)
        {
            // factor 3 for xyz
            pairForcesDim[key][atom] = 3 * neighborList->get_size(atom);
        }
    }

    return;
}

void Predictor::allocate_pairForces(const DescriptorCollector& descriptorCollection,
                                    const Size&                numberAtoms)
{
    compute_pairForcesDim(descriptorCollection);
    if (global::parallel.off()) allocate_pairForcesCPU(numberAtoms);
    else allocate_pairForcesGPU(numberAtoms);

    return;
}

void Predictor::allocate_pairForcesCPU(const Size& numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        pairForces[key]->resize(numberAtoms);
        for (Size atom = 0; atom < pairForces[key]->size(); atom++)
        {
            // factor 3 for xyz
            (*pairForces[key])[atom].resize(pairForcesDim[key][atom]);
            std::fill((*pairForces[key])[atom].begin(), (*pairForces[key])[atom].end(), (Real)0);
        }
    }

    return;
}

void Predictor::allocate_centralForces(const Int numberIons)
{
    if (global::parallel.off()) allocate_centralForcesCPU(numberIons);
    else allocate_centralForcesGPU(numberIons);

    return;
}

void Predictor::allocate_centralForcesCPU(const Int numberIons)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        centralForces[key]->resize(3 * numberIons);
        std::fill(centralForces[key]->begin(), centralForces[key]->end(), (Real)0.0);
    }

    return;
}

void Predictor::allocate_stressTensor()
{
    if (global::parallel.off()) allocate_stressTensorCPU();
    else allocate_stressTensorGPU();

    return;
}

void Predictor::allocate_stressTensorCPU()
{
    for (const String& key : constants::descriptorKeyList)
    {
        stressTensor[key]->resize(9);
        std::fill(stressTensor[key]->begin(), stressTensor[key]->end(), (Real)0);
    }
    totalStressTensor.resize(9);
    std::fill(totalStressTensor.begin(), totalStressTensor.end(), (Real)0);

    return;
}

void Predictor::allocate_atomicForces()
{
    if (global::parallel.off()) allocate_atomicForcesCPU();
    else allocate_atomicForcesGPU();

    return;
}

void Predictor::allocate_atomicForces(const Size& nAtoms)
{
    if (global::parallel.off()) allocate_atomicForcesCPU(nAtoms);
    else allocate_atomicForcesGPU(nAtoms);

    return;
}

void Predictor::allocate_atomicForcesCPU()
{
    // all descriptors have to have the same number of central atoms
    // choose one which is initialized
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        atomicForces.resize(centralForces[key]->size());
        break;
    }
    std::fill(atomicForces.begin(), atomicForces.end(), (Real)0);

    return;
}

void Predictor::allocate_atomicForcesCPU(const Size& nAtoms)
{
    atomicForces.resize(3 * nAtoms);
    std::fill(atomicForces.begin(), atomicForces.end(), (Real)0);

    return;
}

void Predictor::compute_tempForceVectorDim(const Vec1Int& nAtomsType)
{
    tempForceVectorDim.resize(nAtomsType.size());
    for (Size typeStruc = 0; typeStruc < nAtomsType.size(); typeStruc++)
        tempForceVectorDim[typeStruc] = nAtomsType[typeStruc];

    return;
}

void Predictor::allocate_tempForceVector(const Vec1Int& nAtomsType)
{
    compute_tempForceVectorDim(nAtomsType);
    if (global::parallel.off()) allocate_tempForceVectorCPU();
    else allocate_tempForceVectorGPU();

    return;
}

void Predictor::allocate_tempForceVectorCPU()
{
    tempForceVector.resize(tempForceVectorDim.size());
    for (Size typeStruc = 0; typeStruc < tempForceVectorDim.size(); typeStruc++)
    {
        tempForceVector[typeStruc].resize(tempForceVectorDim[typeStruc]);
    }

    return;
}

void Predictor::update(const Kernel& kernel)
{
    typeMap = kernel.get_typeMap();
    //nvtxRangePush( "updateEnergy" );
    updateEnergy(kernel);
    //nvtxRangePop();
    //nvtxRangePush( "updateForces" );
    updateForces(kernel);
    //nvtxRangePop();

    return;
}

void Predictor::updateEnergy(const Kernel& kernel)
{
    const Vec1Int& nAtomsType = kernel.get_nAtomsType();
    const Int      numberIons = vector_tools::sum(nAtomsType);
    allocate_energyArray(nAtomsType.size());
    // TODO include the LTOTEN_SYSTEM tag for scaling
    Real scaleFactor = (Real)1.0 / (Real)numberIons;
    compute_referenceEnergyTotal(nAtomsType);
    totalEnergy = (Real)0;
    for (Size typeStruc = 0; typeStruc < nAtomsType.size(); typeStruc++)
    {
        Size typeForceField = typeMap->toType(typeStruc);
        compute_energyArray(kernel, typeStruc, typeForceField, nAtomsType[typeStruc], scaleFactor);
        totalEnergy += ((Real)numberIons) * computeEnergyPerType(typeStruc, typeForceField);
    }
    finalizeEnergyComputation(nAtomsType);

    return;
}

void Predictor::compute_energyArray(const Kernel& kernel,
                                    const Size    typeStruc,
                                    const Size    typeForceField,
                                    const Size    atomsPerType,
                                    const Real    scaleFactor)
{
    for (Size atom = 0; atom < atomsPerType; atom++)
    {
        linalg::scaleVectorPlusVector(scaleFactor,
                                      kernel.get_kernelMatrixAtom(typeStruc, atom),
                                      energyArray[typeStruc],
                                      (*numberLocalRefConfs)[typeForceField]);
    }

    return;
}

Real Predictor::computeEnergyPerType(const Size typeStruc, const Size typeForceField)
{
    return linalg::dotProduct(energyArray[typeStruc],
                              regressionCoefficients->get_slice(typeForceField),
                              (*numberLocalRefConfs)[typeForceField]);
}

void Predictor::compute_referenceEnergyTotal(const Vec1Int& atomsPerType)
{
    referenceEnergyTotal = (Real)0;
    for (Size typeStruc = 0; typeStruc < atomsPerType.size(); typeStruc++)
    {
        Size typeFF = typeMap->toType(typeStruc);
        referenceEnergyTotal += atomsPerType[typeStruc] * (*referenceEnergyPerType)[typeFF];
    }

    return;
}

void Predictor::finalizeEnergyComputation(const Vec1Int& atomsPerType)
{
#ifndef VASPML_USE_FIT_DEBUG
    if (energyType == TotalEnergyType::IsolatedAtom) { totalEnergy += referenceEnergyTotal; }
    else if (energyType == TotalEnergyType::AverageTrainEnergy)
    {
        for (Size typeStruc = 0; typeStruc < atomsPerType.size(); typeStruc++)
        {
            Size typeFF = typeMap->toType(typeStruc);
            totalEnergy +=
                atomsPerType[typeStruc] * ((*referenceEnergyPerType)[typeFF] + averageTrainEnergy);
        }
    }
#endif

    return;
}

void Predictor::updateForces(const Kernel& kernel)
{

    const Vec1Int&             nAtomsType = kernel.get_nAtomsType();
    const Int                  numberIons = vector_tools::sum(nAtomsType);
    const DescriptorCollector& descriptorCollection = kernel.get_descriptorCollection();
    Size                       numberAtoms = vector_tools::sum(nAtomsType);

    allocate_derivativeMatrix(nAtomsType);
    allocate_forcePreContract(kernel.get_descriptorCollection(), numberAtoms);
    allocate_pairForces(descriptorCollection, numberAtoms);
    allocate_centralForces(numberIons);
    allocate_tempForceVector(nAtomsType);

    for (Size typeStruc = 0; typeStruc < nAtomsType.size(); typeStruc++)
    {
        Size typeForceField = typeMap->toType(typeStruc);
        compute_derivativeMatrix(kernel, tempForceVector[typeStruc], typeStruc, typeForceField);
    }

    //write_tempForceVector( "TempForceVector.dat" );
    //write_derivativeMatrix( "DerivativeMatrix.dat" );
    compute_forcePreContract(descriptorCollection, *typeMap);
    //write_forcePreContract( "ForcePreContract.dat" );
    computeForceArrays(descriptorCollection);
    //write_pairForces( "pairForces.dat" );
    //write_centralForces( "centralForces.dat" );
    forceArraysComputed = true;

    return;
}

void Predictor::compute_derivativeMatrix(const Kernel& kernel,
                                         Vec1Real&     tempForceVector,
                                         const Size    typeStruc,
                                         const Size    typeForceField)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        kernel.computeDerivativeKernel((*derivativeMatrix[key])[typeStruc],
                                       tempForceVector,
                                       descriptorsRefConfsWeighted[key]->get_slice(typeForceField),
                                       regressionCoefficients->get_slice(typeForceField),
                                       typeStruc,
                                       typeForceField,
                                       key);
    }

    return;
}

void Predictor::compute_forcePreContract(const DescriptorCollector& descriptorCollection,
                                         const TypeMap&             typeMap)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        descriptor.compute_forcePreContract(*derivativeMatrix[key],
                                            *forcePreContract[key],
                                            typeMap);
    }

    return;
}

void Predictor::computeForceArrays(const DescriptorCollector& descriptorCollection)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        descriptor.computeForceTerms(*forcePreContract[key], *pairForces[key], *centralForces[key]);
    }

    return;
}

void Predictor::compute_atomicForces(const DescriptorCollector& descriptorCollection)
{
    if (forceArraysComputed)
    {
        allocate_atomicForces();
        if (global::parallel.off()) compute_atomicForcesCPU(descriptorCollection);
        else compute_atomicForcesGPU(descriptorCollection);
    }
    else
    {
        throw std::runtime_error(
            "ERROR: void Predictor::compute_atomicForces( const DescriptorCollector& "
            "descriptorCollection )\n"
            "force arrays were not computed. Please call void Predictor::updateForcesAndStress( "
            "const Kernel& kernel ) before.");
    }

    return;
}

void Predictor::compute_atomicForces(const DescriptorCollector& descriptorCollection,
                                     const Vec1Int&             globalIonIndex,
                                     const Int&                 nAtomsGlob)
{
    if (forceArraysComputed)
    {
        allocate_atomicForces(nAtomsGlob);
        //    switch (algoExecution)
        //    {
        //    case ExecutionPolicy::cpuSingleCore:
        // copy to CPU
        //for ( const auto& key : constants::descriptorKeyList )
        //{
        //    const NearestNeighborNSquare neighborList =
        //        *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        //    cudaMemPrefetchAsync( centralForces[ key ]->data(),
        //                          centralForces[ key ]->size() *
        //                           sizeof( Real ), cudaCpuDeviceId );
        //    for ( Size centralAtom = 0; centralAtom < neighborList.get_nAtoms(); centralAtom++ )
        //    {
        //        cudaMemPrefetchAsync( (*pairForces[ key ])[ centralAtom ].data(),
        //                               (*pairForces[ key ])[ centralAtom ].size() *
        //                               sizeof( Real ), cudaCpuDeviceId );
        //    }
        //}
        compute_atomicForcesCPU(descriptorCollection, globalIonIndex);
        //// copy back to GPU
        //for ( const auto& key : constants::descriptorKeyList )
        //{
        //    const NearestNeighborNSquare neighborList =
        //        *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        //    cudaMemPrefetchAsync( centralForces[ key ]->data(),
        //                          centralForces[ key ]->size() *
        //                           sizeof( Real ), 0 );
        //    for ( Size centralAtom = 0; centralAtom < neighborList.get_nAtoms(); centralAtom++ )
        //    {
        //        cudaMemPrefetchAsync( (*pairForces[ key ])[ centralAtom ].data(),
        //                               (*pairForces[ key ])[ centralAtom ].size() *
        //                               sizeof( Real ), 0 );
        //    }
        //}

        //        break;
        //    case ExecutionPolicy::gpuStdLib:
        //        VASPML_TBD_PARALLEL(
        //            compute_atomicForcesGPU( descriptorCollection, globalIonIndex, nAtomsGlob );
        //        );
        //        break;
        //    }
    }
    else
    {
        throw std::runtime_error(
            "ERROR: void Predictor::compute_atomicForces( const DescriptorCollector& "
            "descriptorCollection )\n"
            "force arrays were not computed. Please call void Predictor::updateForcesAndStress( "
            "const Kernel& kernel ) before.");
    }

    return;
}

void Predictor::compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptorCollection.getDescriptor(key).get_weight() <= 0) continue;
        const NearestNeighborNSquare neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec2Real&      pairForces = (*this->pairForces[key]);
        Vec1Real&      centralForces = (*this->centralForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();
        for (Size centralAtom = 0; centralAtom < neighborList.get_nAtoms(); centralAtom++)
        {
            atomicForces[3 * centralAtom] += centralForces[3 * centralAtom];
            atomicForces[3 * centralAtom + 1] += centralForces[3 * centralAtom + 1];
            atomicForces[3 * centralAtom + 2] += centralForces[3 * centralAtom + 2];
            SumFunctorForce functor(atomicForces, pairForces[centralAtom]);
            std::for_each(neighborIndex[centralAtom].cbegin(),
                          neighborIndex[centralAtom].cend(),
                          functor);
        }
    }

    return;
}

void Predictor::compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection,
                                        const Vec1Int&             globalIonIndex)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptorCollection.getDescriptor(key).get_weight() <= 0) continue;
        const NearestNeighborNSquare neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec2Real&      pairForces = (*this->pairForces[key]);
        Vec1Real&      centralForces = (*this->centralForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();
        for (Size centralAtom = 0; centralAtom < neighborList.get_nAtoms(); centralAtom++)
        {
            const Int& globIndx = globalIonIndex[centralAtom];
            atomicForces[3 * globIndx] += centralForces[3 * centralAtom] * constants::FUNIT;
            atomicForces[3 * globIndx + 1] += centralForces[3 * centralAtom + 1] * constants::FUNIT;
            atomicForces[3 * globIndx + 2] += centralForces[3 * centralAtom + 2] * constants::FUNIT;

            SumFunctorForce functor(atomicForces, pairForces[centralAtom]);
            std::for_each(neighborIndex[centralAtom].cbegin(),
                          neighborIndex[centralAtom].cend(),
                          functor);
        }
    }

    return;
}

void Predictor::compute_stressTensor(const DescriptorCollector& descriptorCollection,
                                     const Real&                volume)
{
    if (forceArraysComputed)
    {
        allocate_stressTensor();
        if (global::parallel.off()) compute_stressTensorCPU(descriptorCollection, volume);
        else compute_stressTensorGPU(descriptorCollection, volume);
    }
    else
    {
        global::tutor.bug("ERROR: void Predictor::compute_atomicForces( const "
                          "DescriptorCollector& descriptorCollection )\n"
                          "force arrays were not computed. Please call void "
                          "Predictor::updateForcesAndStress( const Kernel& kernel ) before.");
    }

    return;
}

void Predictor::compute_stressTensorCPU(const DescriptorCollector& descriptorCollection,
                                        const Real&                volume)
{
    const Real& inverseVolume = (Real)-1.0 / volume;
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptorCollection.getDescriptor(key).get_weight() <= 0) continue;
        const NearestNeighborNSquare& neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec1Real&      stressTensor = (*this->stressTensor[key]);
        Vec2Real&      pairForces = (*this->pairForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();
        const Vec2Real connectionVector = neighborList.get_connectionVector();
        const Size     nAtoms = neighborList.get_nAtoms();

        for (Size centralAtom = 0; centralAtom < nAtoms; centralAtom++)
        {
            SumFunctorStress functor(stressTensor,
                                     pairForces[centralAtom],
                                     connectionVector[centralAtom],
                                     inverseVolume);
            std::for_each(neighborIndex[centralAtom].cbegin(),
                          neighborIndex[centralAtom].cend(),
                          functor);
        }
        // symmetrize
        stressTensor[3] = stressTensor[1];
        stressTensor[6] = stressTensor[2];
        stressTensor[7] = stressTensor[5];

        totalStressTensor[0] += stressTensor[0];
        totalStressTensor[1] += stressTensor[1];
        totalStressTensor[2] += stressTensor[2];
        totalStressTensor[4] += stressTensor[4];
        totalStressTensor[5] += stressTensor[5];
        totalStressTensor[8] += stressTensor[8];
    }
    // symmetrize total stress
    totalStressTensor[3] = totalStressTensor[1];
    totalStressTensor[6] = totalStressTensor[2];
    totalStressTensor[7] = totalStressTensor[5];

    return;
}

const Real& Predictor::get_pairForcesX(const String& key,
                                       const Size    centralAtom,
                                       const Size    neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesX( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesX( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );

    return (*pairForces.at(key))[centralAtom][3 * neighborAtom];
}

const Real& Predictor::get_pairForcesY(const String& key,
                                       const Size    centralAtom,
                                       const Size    neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesY( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesY( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );

    return (*pairForces.at(key))[centralAtom][3 * neighborAtom + 1];
}

const Real& Predictor::get_pairForcesZ(const String& key,
                                       const Size    centralAtom,
                                       const Size    neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesZ( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesZ( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );

    return (*pairForces.at(key))[centralAtom][3 * neighborAtom + 2];
}

const std::tuple<const Real&, const Real&, const Real&>
Predictor::get_pairForces(const String& key, const Size centralAtom, const Size neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForces( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForces( const String& key, const "
                "Size centralAtom, const Size neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );

    return std::tie(get_pairForcesX(key, centralAtom, neighborAtom),
                    get_pairForcesY(key, centralAtom, neighborAtom),
                    get_pairForcesZ(key, centralAtom, neighborAtom));
}

const std::map<String, ShVec2Real>& Predictor::get_pairForces() const
{
    return pairForces;
}

const Real& Predictor::get_centralForcesX(const String& key, const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesX( const "
                                     "String& key, const Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return (*centralForces.at(key))[3 * centralAtom];
}

const Real& Predictor::get_centralForcesY(const String& key, const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesY( const "
                                     "String& key, const Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return (*centralForces.at(key))[3 * centralAtom + 1];
}

const Real& Predictor::get_centralForcesZ(const String& key, const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesZ( const "
                                     "String& key, const Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return (*centralForces.at(key))[3 * centralAtom + 2];
}

const std::tuple<const Real&, const Real&, const Real&> Predictor::get_centralForces(
    const String& key,
    const Size    centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForces( const "
                                     "String& key, const Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return std::tie(get_centralForcesX(key, centralAtom),
                    get_centralForcesY(key, centralAtom),
                    get_centralForcesZ(key, centralAtom));
}

const std::map<String, ShVec1Real>& Predictor::get_centralForces() const
{
    return centralForces;
}

const Real& Predictor::get_atomicForcesX(const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces.size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesX( const "
                                     "Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return atomicForces[3 * centralAtom];
}

const Real& Predictor::get_atomicForcesY(const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces.size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesY( const "
                                     "Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return atomicForces[3 * centralAtom + 1];
}

const Real& Predictor::get_atomicForcesZ(const Size centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces.size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesZ( const "
                                     "Size centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );

    return atomicForces[3 * centralAtom + 2];
}

const Vec1Real& Predictor::get_atomicForces() const
{
    return atomicForces;
}

const Real& Predictor::get_stressTensor(const String& key, const Size indx0, const Size indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_stressTensor( const String& key, const Size indx0, "
                "const Size indx1 )const\n"
                "Index indx0 is larger than 3 ( only 3 space dimensions available )\n");
        }
        if (indx1 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_stressTensor( const String& key, const Size indx0, "
                "const Size indx1 )const\n"
                "Index indx1 is larger than 3 ( only 3 space dimensions available )\n");
        }
    );

    return (*stressTensor.at(key))[indx0 * 3 + indx1];
}

const Real& Predictor::get_totalStressTensor(const Size indx0, const Size indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_totalStressTensor( const Size indx0, const "
                "Size indx1 )const\n"
                "Index indx0 is larger than 3 ( only 3 space dimensions available )\n");
        }
        if (indx1 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_totalStressTensor( const Size indx0, const "
                "Size indx1 )const\n"
                "Index indx1 is larger than 3 ( only 3 space dimensions available )\n");
        }
    );

    return totalStressTensor[indx0 * 3 + indx1];
}

const Vec1Real& Predictor::get_totalStressTensor() const
{
    return totalStressTensor;
}

const Real& Predictor::get_totalEnergy() const
{
    return totalEnergy;
}

void Predictor::write_atomicForceToScreen() const
{
    for (Size atom = 0; atom < atomicForces.size(); atom += 3)
    {
        std::cout << str("%-24.16E %-24.16E %-24.16E\n",
                         atomicForces[atom] * constants::FUNIT,
                         atomicForces[atom + 1] * constants::FUNIT,
                         atomicForces[atom + 2] * constants::FUNIT);
    }

    return;
}

void Predictor::write_tempForceVector(const String& fname) const
{
    file_io::writeRealVec2D(tempForceVector, fname);

    return;
}

void Predictor::write_derivativeMatrix(const String& baseName) const
{
    for (const String& key : constants::descriptorKeyList)
    {
        file_io::writeRealVec2D(*derivativeMatrix.at(key), baseName + key);
    }

    return;
}

void Predictor::write_forcePreContract(const String& baseName) const
{
    for (const String& key : constants::descriptorKeyList)
    {
        file_io::writeRealVec2D(*forcePreContract.at(key), baseName + key);
    }

    return;
}

void Predictor::write_pairForces(const String& baseName) const
{
    for (const String& key : constants::descriptorKeyList)
    {
        file_io::writeRealVec2D(*pairForces.at(key), baseName + key);
    }

    return;
}

void Predictor::write_centralForces(const String& baseName) const
{
    for (const String& key : constants::descriptorKeyList)
    {
        auto file = file_io::openFileO(baseName + key);
        for (const auto& x : *centralForces.at(key)) file << str("%24.16E ", x) << std::endl;
        file.close();
    }

    return;
}
