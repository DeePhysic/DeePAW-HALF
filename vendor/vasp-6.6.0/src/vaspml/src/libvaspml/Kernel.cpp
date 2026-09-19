#include "Kernel.hpp"

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "Frame.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "ShmemArray.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm> // for fill, for_each
#include <fstream>   // for basic_ostream, basic_ofstream
#include <stdexcept> // for runtime_error


using namespace vaspml;
using namespace vaspml::setup;

namespace vaspml
{

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> generate_ShmemArrayMap(
    const std::map<String, Real>&      weights,
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::shared_ptr<MlMPI>&      mpiIn)
{

    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> arrays;
    for (auto& [key, item] : featureSpaceSize)
    {
        arrays[key] = std::make_shared<ShmemArray2DVariableLen<Real>>(
            weights.at(key) > 0 ? *featureSpaceSize.at(key) : Vec1Int{0},
            mpiIn);
    }

    return arrays;
}

} //namespace vaspml

MapString vaspml::dataMapKernel()
{
    return MapString{
        {"kernelMatrix", "Vec2Real"}
    };
}

Kernel::Kernel(ShRec record) :
    data(assignOrMakeRecord(record, dataMapKernel())),
    kernelMatrix(data->get<Vec2Real>("kernelMatrix"))
{}

Kernel::Kernel(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn, ShRec record) :
    data(assignOrMakeRecord(record, dataMapKernel())),
    mpi(mpiIn),
    weights(setup::make_weightsMap(inputParameters)),
    numberLocalRefConfs(std::make_shared<Vec1Int>(inputParameters.cget<Vec1Int>("numLrc"))),
    numberDescriptors(setup::make_numberDescriptorsMap(inputParameters)),
    featureSpaceSize(setup::generate_featureSpaceMap(inputParameters)),
    descriptorsRefConfs(generate_ShmemArrayMap(weights, featureSpaceSize, mpiIn)),
    descriptorCollection(
        std::make_shared<DescriptorCollector>(0, std::map<String, ShVec2Real>(), nullptr)),
    kernelMatrix(data->get<Vec2Real>("kernelMatrix"))
{
    // filling descriptors into shared memory arrays
    fillDescriptor(inputParameters);
}

Kernel::Kernel(const Record&                 incar,
               const Record&                 mlab,
               const std::shared_ptr<MlMPI>& mpiIn,
               ShRec                         record) :
    data(assignOrMakeRecord(record, dataMapKernel())),
    mpi(mpiIn),
    weights(setup::make_weightsMap(incar)),
    numberLocalRefConfs(std::make_shared<Vec1Int>(mlab.cget<Vec1Int>("numLrc"))),

    //descriptorsRefConfs(generate_ShmemArrayMap(weights, featureSpaceSize, mpiIn)),
    descriptorCollection(
        std::make_shared<DescriptorCollector>(0, std::map<String, ShVec2Real>(), nullptr)),
    kernelMatrix(data->get<Vec2Real>("kernelMatrix"))
{}

void Kernel::init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn)
{
    weights = setup::make_weightsMap(inputParameters);
    mpi = mpiIn;
    numberLocalRefConfs = std::make_shared<Vec1Int>(inputParameters.cget<Vec1Int>("numLrc"));
    numberDescriptors = setup::make_numberDescriptorsMap(inputParameters);
    featureSpaceSize = setup::generate_featureSpaceMap(inputParameters);
    descriptorsRefConfs = generate_ShmemArrayMap(weights, featureSpaceSize, mpi);
    descriptorCollection =
        std::make_shared<DescriptorCollector>(0, std::map<String, ShVec2Real>(), nullptr);
    fillDescriptor(inputParameters);
}

void Kernel::fillDescriptor(const Record& inputParameters)
{
    // filling descriptors into shared memory arrays
    for (const String& key : constants::descriptorKeyList)
    {
        String tag = key + "-reference-configs";
        if (weights[key] > 0)
            //descriptorsRefConfs[key]->set_value(inputParameters[tag].dcget<ShVec2Real>());
            descriptorsRefConfs[key]->set_value(inputParameters.cget<Vec2Real>(tag));
    }
}

void Kernel::updateKernel(const Frame& frame)
{
    const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors =
        frame.get_descriptors();
    VASPML_DEBUG_L1(
        for (Size i = 0; i < constants::descriptorKeyList.size() - 1; i++)
        {
            for (Size j = i + 1; j < constants::descriptorKeyList.size(); j++)
            {
                if (structureDescriptors.at(constants::descriptorKeyList[i])->get_weight() > 0
                    and structureDescriptors.at(constants::descriptorKeyList[j])->get_weight() > 0)
                {
                    if (structureDescriptors.at(constants::descriptorKeyList[i])->get_nAtoms()
                        != structureDescriptors.at(constants::descriptorKeyList[j])->get_nAtoms())
                    {
                        throw std::runtime_error(
                            "ERROR: void Kernel::update( const Frame& frame ) \n"
                            "the number of atoms does not agree in your descriptors");
                    }
                }
            }
        }
    );

    // determine number of atoms per type. At same time check if there was a valid descriptor set.
    [[maybe_unused]] bool descFound = false;
    for (const auto& [key, desc] : structureDescriptors)
    {
        if (desc->get_weight() > 0)
        {
            nAtomsPerType = desc->get_nAtomsType();
            descFound = true;
            break;
        }
    }
    VASPML_DEBUG_L1(
        if (!descFound)
        {
            global::tutor.error("ERROR: " + flf(VASPML_FLF) + " no valid descriptor found!");
        }
    );

    typeMap = frame.get_typeMap();
    allocate_kernelMatrix(nAtomsPerType);
    //nvtxRangePush( "setUpDescriptorCollector" );
    setUpDescriptorCollector(structureDescriptors);
    //nvtxRangePop();
    //nvtxRangePush( "compute_kernelMatrix" );
    compute_kernelMatrix();
    //nvtxRangePop();
}

void Kernel::compute_kernelMatrixDim(const Vec1Int& nAtomPerType)
{
    kernelMatrixDim.resize(typeMap->countStructureTypes());
    for (Size strucType = 0; strucType < typeMap->countStructureTypes(); strucType++)
    {
        Int forceFieldType = typeMap->toType(strucType);
        kernelMatrixDim[strucType] =
            nAtomPerType[strucType] * (*numberLocalRefConfs)[forceFieldType];
    }
}

void Kernel::allocate_kernelMatrix(const Vec1Int& nAtomPerType)
{
    compute_kernelMatrixDim(nAtomPerType);
    if (global::parallel.off())
    {
        kernelMatrix.resize(typeMap->countStructureTypes());
        for (Size strucType = 0; strucType < kernelMatrix.size(); strucType++)
        {
            kernelMatrix[strucType].resize(kernelMatrixDim[strucType], (Real)0);
            std::fill(kernelMatrix[strucType].begin(), kernelMatrix[strucType].end(), (Real)0);
        }
    }
    else { resizeKernelMatrixGPU(nAtomPerType); }
}

void Kernel::resizeKernelMatrixGPU(const Vec1Int& nAtomPerType)
{
    bool resize = kernelMatrixSize.checkResize1Dim(nAtomPerType.size());
    if (resize) { kernelMatrixSize.resizeArray1Dim(kernelMatrix, kernelMatrixDim); }
    else
    {
        resize = kernelMatrixSize.checkResize2Dim(kernelMatrixDim);
        if (resize) kernelMatrixSize.resizeArray2Dim(kernelMatrix, kernelMatrixDim);
    }
    // fill with zeros
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          kernelMatrix.begin(),
                          kernelMatrix.begin() + kernelMatrixSize.act1Dim,
                          [](Vec1Real& slice)
                          { std::fill(slice.begin(), slice.end(), (Real)0.0); });
        },
        __func__);
}

void Kernel::setUpDescriptorCollector(
    const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors)
{
    descriptorCollection->setDescriptorMap(structureDescriptors);
    descriptorCollection->updateCollector();
    descriptorCollection->rearrangeSHS2Body(*typeMap);
}

void Kernel::compute_kernelMatrix()
{
    // unpack variables
    const std::map<String, ShVec2Real>& descriptorStructure =
        descriptorCollection->get_descriptorsNormalized();
    Vec1Int& numberLocalRefConfs = (*this->numberLocalRefConfs);

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec1Int& numberDescriptors = (*this->numberDescriptors[key]);
        // extract descriptor of structure
        const Vec2Real& descriptor = (*descriptorStructure.at(key));
        for (Size strucType = 0; strucType < nAtomsPerType.size(); strucType++)
        {
            Int            forceFieldType = typeMap->toType(strucType);
            constexpr Real one = (Real)1;
            linalg::matMul(linalg::Transpose::NOTRANS,
                           linalg::Transpose::TRANS,
                           nAtomsPerType[strucType],            // m
                           numberLocalRefConfs[forceFieldType], // n
                           numberDescriptors[forceFieldType],   // k
                           one,
                           descriptor[strucType],             // A\in\mathbb{R}^{Natom/times/Ndesc}
                           numberDescriptors[forceFieldType], // k
                           descriptorsRefConfs[key]->get_slice(
                               forceFieldType), //B^{T}\in\mathbb{R}^{Nref/times/Ndesc}
                           numberDescriptors[forceFieldType], // k
                           one,
                           kernelMatrix[strucType],
                           numberLocalRefConfs[forceFieldType]);
        }
    }
}

const std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>>&
Kernel::get_descriptorsRefConfs() const
{
    return descriptorsRefConfs;
}

ShVec1Int Kernel::get_numberLocalRefConfs() const
{
    return numberLocalRefConfs;
}

const std::map<String, ShVec1Int>& Kernel::get_numberDescriptors() const
{
    return numberDescriptors;
}

const std::map<String, ShVec1Int>& Kernel::get_featureSpaceSize() const
{
    return featureSpaceSize;
}

const std::map<String, Real>& Kernel::get_weights() const
{
    return weights;
}

const Vec2Real& Kernel::get_kernelMatrix() const
{
    return kernelMatrix;
}

const Vec1Int& Kernel::get_nAtomsType() const
{
    return nAtomsPerType;
}

std::shared_ptr<const TypeMap> Kernel::get_typeMap() const
{
    return typeMap;
}

const Real* Kernel::get_kernelMatrixAtom(Size type, Size atomInType) const
{
    VASPML_DEBUG_L1(
        if (type >= nAtomsPerType.size())
        {
            throw std::runtime_error("ERROR: const Real* Kernel::get_kernelMatrixAtom( Size type,"
                                     "Size atomInType )\n       type index "
                                     + std::to_string(type) + " out of bounds\n");
        }
    );
    VASPML_DEBUG_L1(
        if (atomInType >= (Size)nAtomsPerType[type])
        {
            throw std::runtime_error("ERROR: const Real* Kernel::get_kernelMatrixAtom( Size type,"
                                     "Size atomInType )\n       atomInType index "
                                     + std::to_string(atomInType) + " out of bounds\n");
        }
    );

    Real* ptr;
    Size  index = atomInType * (*numberLocalRefConfs)[typeMap->toType(type)];
    ptr = &(kernelMatrix)[type][index];
    return ptr;
}

const DescriptorCollector& Kernel::get_descriptorCollection() const
{
    return *descriptorCollection;
}

const std::shared_ptr<DescriptorCollector>& Kernel::get_descriptorCollectionPtr() const
{
    return descriptorCollection;
}

void Kernel::write_kernelMatrix(const String& fname) const
{
    auto file = file_io::openFileO(fname);
    if (global::parallel.off())
    {
        for (Size i = 0; i < kernelMatrix.size(); i++)
            for (Size j = 0; j < kernelMatrix[i].size(); j++)
                file << str("%24.16E ", kernelMatrix[i][j]) << std::endl;
    }
    else
    {
        for (Size i = 0; i < kernelMatrixDim.size(); i++)
            for (Size j = 0; j < kernelMatrixDim[i]; j++)
                file << str("%24.16E ", kernelMatrix[i][j]) << std::endl;
    }
    file.close();
}

void Kernel::set_descriptorsRefConfs(
    const std::shared_ptr<DescriptorCollector>& referenceConfigurations)
{
    const auto& normedDescriptors = referenceConfigurations->get_descriptorsNormalized();
    for (auto& [key, array] : descriptorsRefConfs)
    {
        for (Size type = 0; type < array->get_size0(); type++)
        {
            for (Size locRef = 0; locRef < array->get_size1(type); locRef++)
            {
                array->set_value(type, locRef, (*normedDescriptors.at(key))[type][locRef]);
            }
        }
        //array->writeToFile(  key + "RefConfs.dat" );
    }
}

void Kernel::generate_featureSpaceMap()
{
    const Vec1Int& locRef = *numberLocalRefConfs;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] > 0)
            featureSpaceSize[key] = std::make_shared<Vec1Int>(
                math::elementwiseProduct(locRef, *(numberDescriptors[key])));
    }
}

void Kernel::set_numberDescriptors(const std::map<String, ShVec1Int>& numberDescriptors)
{
    this->numberDescriptors = numberDescriptors;
}

void Kernel::make_ShmemArrayMap()
{
    descriptorsRefConfs = generate_ShmemArrayMap(weights, featureSpaceSize, mpi);
}
