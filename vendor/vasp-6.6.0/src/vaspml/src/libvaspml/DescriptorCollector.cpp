#include "DescriptorCollector.hpp"

#include "Descriptor.hpp"
#include "DescriptorSHS2.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm> // for fill, for_each, transform
#include <cmath>     // for sqrt
#include <fstream>   // for basic_ofstream, basic_ostream
#include <numeric>   // for iota

using namespace vaspml;

template<>
SmartEnum<DescriptorStorage>::EnumMap SmartEnum<DescriptorStorage>::mapEnumsNames()
{
    return SmartEnum<DescriptorStorage>::EnumMap{
        {DescriptorStorage::Type,        "TypeOrder"       },
        {DescriptorStorage::CentralAtom, "CentralAtomOrder"}
    };
}

MapString vaspml::dataMapDescriptorCollector()
{
    return MapString{
        {"normsAtom", "Vec1Real"}
    };
}

DescriptorCollector::DescriptorCollector(Int                                 storageOrder,
                                         const std::map<String, ShVec2Real>& descriptorsNormalized,
                                         ShRec                               record) :
    data(assignOrMakeRecord(record, dataMapDescriptorCollector())),
    normsAtom(data->get<Vec1Real>("normsAtom")),
    descriptorKeyList(constants::descriptorKeyList)
{
    switch (storageOrder)
    {
    case 0:
        this->storageOrder = DescriptorStorage::Type;
        break;
    case 1:
        this->storageOrder = DescriptorStorage::CentralAtom;
        break;
    default:
        global::tutor.bug("ERROR: " + flf(VASPML_FLF)
                          + " A storage Order was chose which does not exist > 1.");
        break;
    }
    if (descriptorsNormalized.empty())
    {
        this->descriptorsNormalized["SHS2-2-body"] = std::make_shared<Vec2Real>();
        this->descriptorsNormalized["SHS3-3-body"] = std::make_shared<Vec2Real>();
    }
    else { this->descriptorsNormalized = descriptorsNormalized; }
    for (const String& key : constants::descriptorKeyList)
    {
        length_descriptorsNormalized[key] = std::make_shared<Vec1Int>();
        descriptorsNormalizedSize[key] = ArrayResizing2D();
        length_descriptorsNormalizedSize[key] = ArrayResizing1D();
    }
}

void DescriptorCollector::setDescriptorMap(
    const std::map<String, std::shared_ptr<Descriptor>>& descriptors)
{
    this->descriptors = descriptors;
}

void DescriptorCollector::updateCollector()
{
    VASPML_DEBUG_L1(
        for (Size i = 0; i < constants::descriptorKeyList.size() - 1; i++)
        {
            for (Size j = i + 1; j < constants::descriptorKeyList.size(); j++)
            {
                if (descriptors.at(constants::descriptorKeyList[i])->get_weight() > 0
                    && descriptors.at(constants::descriptorKeyList[j])->get_weight() > 0)
                {
                    if (descriptors.at(constants::descriptorKeyList[i])->get_nAtoms()
                        != descriptors.at(constants::descriptorKeyList[j])->get_nAtoms())
                    {
                        global::tutor.bug("ERROR: DescriptorCollector::updateCollector(): "
                                          "the number of atoms does not agree in your descriptors");
                    }
                }
            }
        }
    );
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptors.at(key)->get_weight() > 0)
        {
            numberAtoms = descriptors.at(key)->get_nAtoms();
            break;
        }
    }
    allocateArrays();
    normalizeDescriptors();
}

void DescriptorCollector::normalizeDescriptors()
{
    if (global::parallel.off()) normalizeDescriptorsCPU();
    else normalizeDescriptorsGPU();
}

void DescriptorCollector::normalizeDescriptorsCPU()
{
    for (Size i = 0; i < (Size)numberAtoms; i++)
    {
        Real norm = computeNormSingleAtom(i);
        if (storageOrder == DescriptorStorage::Type)
        {
            normalizeDescriptorsSingleAtomTypeOrder(i, norm);
        }
        else if (storageOrder == DescriptorStorage::CentralAtom)
        {
            normalizeDescriptorsSingleAtomCentralOrder(i, norm);
        }
        normsAtom[i] = norm;
    }
}

void DescriptorCollector::normalizeDescriptorsGPU()
{
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          centralAtomIndex.cbegin(),
                          centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                          [&](const Size& atom) mutable
                          {
                              Real norm = computeNormSingleAtomGPU(atom);
                              normalizeDescriptorsSingleAtomTypeOrder(atom, norm);
                              normsAtom[atom] = norm;
                          });
        },
        __func__);
}

Real DescriptorCollector::computeNormSingleAtomGPU(const Size& atomIndex)
{
    Real normTotal = 0;
    for (const String& key : descriptorKeyList)
    {
        Real weight = descriptors.at(key)->get_weight();
        if (weight > 0)
        {
            //const Vec1Real& data = descriptors.at( key ) -> get_descriptor( atomIndex );
            // TODO: Shouldn't this use the linalg:: version (on GPU)?
            Real normDescriptor = math::l2Norm(descriptors.at(key)->get_descriptor(atomIndex),
                                               descriptors.at(key)->get_sizeDescriptor(atomIndex));
            normTotal += weight * normDescriptor * normDescriptor;
        }
    }

    return sqrt(normTotal);
}

Real DescriptorCollector::computeNormSingleAtom(Size atomIndex)
{
    Real normTotal = 0;
    for (const String& key : descriptorKeyList)
    {
        Real weight = descriptors.at(key)->get_weight();
        if (weight > 0)
        {
            // const Vec1Real& data = descriptors.at( key ) -> get_descriptor( atomIndex );
            Real normDescriptor =
                linalg::l2Norm(descriptors.at(key)->get_descriptor(atomIndex),
                               descriptors.at(key)->get_sizeDescriptor(atomIndex));
            normTotal += weight * normDescriptor * normDescriptor;
        }
    }

    return sqrt(normTotal);
}

void DescriptorCollector::normalizeDescriptorsSingleAtomTypeOrder(const Size& atomIndex,
                                                                  const Real  norm)
{
    for (const String& key : descriptorKeyList)
    {
        Real weight = descriptors.at(key)->get_weight();
        if (weight <= (Real)0) continue;

        const Int& type = descriptors.at(key)->get_typeIndexCentral(atomIndex);
        const Int& sizeDesc = descriptors.at(key)->get_sizeDescriptor(atomIndex);
        const Int& atomIndexPerType =
            descriptors.at(key)->get_neighborList().get_centralAtomIndexPerType()[atomIndex];
        if (descriptors.at(key)->get_isNormalized())
        {
            if (norm > constants::EPS_TOL)
            {
                // note gpu needs .at( key ) access for std::maps
                Real            factor = std::sqrt(weight) / norm;
                const Vec1Real& desc = descriptors.at(key)->get_descriptor(atomIndex);
                Vec1Real&       descriptorsNormalized =
                    (*this->descriptorsNormalized.at(key))[type]; // race condition
                Int& length_descriptorsNormalized =
                    (*this->length_descriptorsNormalized.at(key))[type];
                std::transform(VASPML_POLICY_SEQ
                               desc.cbegin(),
                               desc.cend(),
                               descriptorsNormalized.begin() + sizeDesc * atomIndexPerType,
                               [&factor](const Real& x) { return x * factor; });
                length_descriptorsNormalized += desc.size();
            }
        }
    }
}

void DescriptorCollector::normalizeDescriptorsSingleAtomCentralOrder(Size       atomIndex,
                                                                     const Real norm)
{
    for (const String& key : descriptorKeyList)
    {
        if (descriptors.at(key)->get_isNormalized())
        {
            Real weight = descriptors.at(key)->get_weight();
            if (weight > (Real)0)
            {
                if (norm > constants::EPS_TOL)
                {
                    Real            factor = sqrt(weight) / norm;
                    const Vec1Real& desc = descriptors.at(key)->get_descriptor(atomIndex);
                    // rescale descriptors
                    Vec1Real& descriptorsNormalized =
                        (*this->descriptorsNormalized[key])[atomIndex];
                    Int& length_descriptorsNormalized =
                        (*this->length_descriptorsNormalized[key])[atomIndex];
                    std::for_each( // par_unseq,
                        desc.cbegin(),
                        desc.cend(),
                        [&](const Real& desc)
                        {
                            descriptorsNormalized[length_descriptorsNormalized] = desc * factor;
                            length_descriptorsNormalized++;
                        });
                }
            }
        }
    }
}

void DescriptorCollector::allocateArrays()
{
    if (global::parallel.off()) { allocateArraysCPU(); }
    else { allocateArraysGPU(); }
}

void DescriptorCollector::allocateArraysCPU()
{
    switch (storageOrder)
    {
    case DescriptorStorage::Type:
        for (const String& key : constants::descriptorKeyList)
        {
            if (descriptors.at(key)->get_weight() <= 0) continue;
            const Vec1Int& nAtomsPerType = descriptors.at(key)->get_nAtomsType();
            length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            std::fill(length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->end(),
                      (Int)0.0);
            descriptorsNormalized[key]->resize(nAtomsPerType.size());
            // number of descriptors is same for all atoms for certain descriptor
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            for (Size type = 0; type < nAtomsPerType.size(); type++)
            {
                (*descriptorsNormalized[key])[type].resize(nAtomsPerType[type] * numberDescriptors);
                std::fill((*descriptorsNormalized[key])[type].begin(),
                          (*descriptorsNormalized[key])[type].end(),
                          (Real)0.0);
            }
        }
        break;
    case DescriptorStorage::CentralAtom:
        for (const String& key : constants::descriptorKeyList)
        {
            if (descriptors.at(key)->get_weight() <= 0) continue;
            const Size& nAtoms = descriptors.at(key)->get_nAtoms();
            descriptorsNormalized[key]->resize(nAtoms);
            length_descriptorsNormalized[key]->resize(nAtoms);
            std::fill(length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->end(),
                      (Int)0.0);
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            for (Size atom = 0; atom < nAtoms; atom++)
            {
                // number of descriptors is same for all atoms for certain descriptor
                (*descriptorsNormalized[key])[atom].resize(numberDescriptors);
            }
        }
        break;
    }
    normsAtom.resize(numberAtoms);
}

void DescriptorCollector::allocateArraysGPU()
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        normsAtom.resize(numberAtoms);
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(),
                  centralAtomIndex.begin() + centralAtomIndexSize.actDim,
                  0);
    }

    switch (storageOrder)
    {
    case DescriptorStorage::Type:
        for (const String& key : constants::descriptorKeyList)
        {
            if (descriptors.at(key)->get_weight() <= 0) continue;
            const Vec1Int& nAtomsPerType = descriptors.at(key)->get_nAtomsType();
            length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            // zero filling has to be shifted to gpu
            // number of descriptors is same for all atoms for certain descriptor
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            resize = descriptorsNormalizedSize[key].checkResize1Dim(nAtomsPerType.size());
            if (resize)
            {
                descriptorsNormalizedSize[key].resizeArray1Dim(*(descriptorsNormalized[key]),
                                                               nAtomsPerType,
                                                               numberDescriptors);
                length_descriptorsNormalizedSize[key].maxDim = nAtomsPerType.size();
                length_descriptorsNormalizedSize[key].actDim = nAtomsPerType.size();
                length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            }
            else
            {
                length_descriptorsNormalizedSize[key].actDim = nAtomsPerType.size();
                resize = descriptorsNormalizedSize[key].checkResize2Dim(nAtomsPerType,
                                                                        numberDescriptors);
                if (resize)
                {
                    descriptorsNormalizedSize[key].resizeArray2Dim(*(descriptorsNormalized[key]),
                                                                   nAtomsPerType,
                                                                   numberDescriptors);
                }
            }
            global::parallel.run(
                [&]([[maybe_unused]] const auto& policy)
                {
                    std::fill(VASPML_POLICY
                              length_descriptorsNormalized[key]->begin(),
                              length_descriptorsNormalized[key]->begin()
                                  + length_descriptorsNormalizedSize[key].actDim,
                              (Int)0);
                },
                __func__);
        }
        break;
    case DescriptorStorage::CentralAtom:
        for (const String& key : constants::descriptorKeyList)
        {
            if (descriptors.at(key)->get_weight() <= 0) continue;
            const Size& nAtoms = descriptors.at(key)->get_nAtoms();
            resize = descriptorsNormalizedSize[key].checkResize1Dim(nAtoms);
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            if (resize)
            {
                descriptorsNormalizedSize[key].resizeArray1Dim(*(descriptorsNormalized[key]),
                                                               numberDescriptors);
                length_descriptorsNormalizedSize[key].maxDim = nAtoms;
                length_descriptorsNormalizedSize[key].actDim = nAtoms;
            }
            else
            {
                resize = descriptorsNormalizedSize[key].checkResize2Dim(numberDescriptors);
                length_descriptorsNormalizedSize[key].actDim = nAtoms;
                if (resize)
                {
                    descriptorsNormalizedSize[key].resizeArray2Dim((*descriptorsNormalized[key]),
                                                                   numberDescriptors);
                }
            }
            length_descriptorsNormalized[key]->resize(nAtoms);
            global::parallel.run(
                [&]([[maybe_unused]] const auto& policy)
                {
                    std::fill(VASPML_POLICY
                              length_descriptorsNormalized[key]->begin(),
                              length_descriptorsNormalized[key]->begin()
                                  + length_descriptorsNormalizedSize[key].actDim,
                              (Int)0);
                },
                __func__);
        }
        break;
    }
}

const std::map<String, ShVec2Real>& DescriptorCollector::get_descriptorsNormalized() const
{
    return descriptorsNormalized;
}

Real DescriptorCollector::get_normAtom(Size atomIndx) const
{
    return normsAtom[atomIndx];
}

const Vec1Real& DescriptorCollector::get_normAtom() const
{
    return normsAtom;
}

void DescriptorCollector::rearrangeSHS2Body(const TypeMap& typeMap)
{
    if (typeMap.countStructureTypes() == typeMap.countForceFieldTypes()) return;
    if (storageOrder == DescriptorStorage::Type)
    {
        if (global::parallel.off()) rearrangeSHS2BodyTypeOrderCPU(typeMap);
        else rearrangeSHS2BodyTypeOrderGPU(typeMap);
    }
    else if (storageOrder == DescriptorStorage::CentralAtom)
    {
        rearrangeSHS2BodyAtomOrderCPU(typeMap);
    }
}

void DescriptorCollector::rearrangeSHS2BodyTypeOrderCPU(const TypeMap& typeMap)
{
    // need to make copy here otherwise data will be overwritten before use
    const Vec2Real descriptorSHS2body = *descriptorsNormalized["SHS2-2-body"];
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    const Vec1Int& nAtomsPerType = descriptors.at("SHS2-2-body")->get_nAtomsType();
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);

    Size centralType = 0;
    std::for_each(descriptorsNormalized.begin(),
                  descriptorsNormalized.begin() + typeMap.countStructureTypes(),
                  [&](Vec1Real& slice)
                  {
                      slice.resize(nAtomsPerType[centralType] * typeMap.countForceFieldTypes()
                                   * numberDescriptors);
                      for (Size centralAtom = 0; centralAtom < (Size)nAtomsPerType[centralType];
                           centralAtom++)
                      {
                          for (Size neighborTypeStruc = 0;
                               neighborTypeStruc < typeMap.countStructureTypes();
                               neighborTypeStruc++)
                          {
                              Size neighborTypeFF = typeMap.toType(neighborTypeStruc);
                              Size neighborOffsetFF = neighborTypeFF * numberDescriptors
                                                    + centralAtom * numberDescriptors
                                                          * typeMap.countForceFieldTypes();
                              Size neighborOffsetStruc =
                                  neighborTypeStruc * numberDescriptors
                                  + centralAtom * numberDescriptors * typeMap.countStructureTypes();
                              for (Size desc = 0; desc < (Size)numberDescriptors; desc++)
                              {
                                  slice[neighborOffsetFF + desc] =
                                      descriptorSHS2body[centralType][neighborOffsetStruc + desc];
                              }
                          }
                      }
                      centralType++;
                  });
}

void DescriptorCollector::rearrangeSHS2BodyTypeOrderGPU(const TypeMap& typeMap)
{
    // need to make copy here otherwise data will be overwritten before use
    const Vec2Real descriptorSHS2body = *descriptorsNormalized["SHS2-2-body"];
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);
    const Vec1Int& typeIndexCentral =
        descriptors.at("SHS2-2-body")->get_neighborList().get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType =
        descriptors.at("SHS2-2-body")->get_neighborList().get_centralAtomIndexPerType();

    // allocate descriptorsNormalized for SHS2-body
    allocateDescriptorSHS2BodyNormalized(typeMap);

    // Now do the actual calculations over numberAtoms
    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
        [&](const Size& centralAtom) mutable
        {
            Size centralType = typeIndexCentral[centralAtom];
            Size centralAtomInType = centralAtomIndexPerType[centralAtom];
            for (Size neighborTypeStruc = 0; neighborTypeStruc < typeMap.countStructureTypes();
                 neighborTypeStruc++)
            {
                Size neighborTypeFF = typeMap.toType(neighborTypeStruc);
                Size neighborOffsetFF =
                    neighborTypeFF * numberDescriptors
                    + centralAtomInType * numberDescriptors * typeMap.countForceFieldTypes();
                Size neighborOffsetStruc =
                    neighborTypeStruc * numberDescriptors
                    + centralAtomInType * numberDescriptors * typeMap.countStructureTypes();
                for (Size desc = 0; desc < (Size)numberDescriptors; desc++)
                {
                    descriptorsNormalized[centralType][neighborOffsetFF + desc] =
                        descriptorSHS2body[centralType][neighborOffsetStruc + desc];
                }
            }
        });
}

void DescriptorCollector::rearrangeSHS2BodyAtomOrderCPU(const TypeMap& typeMap)
{
    const Vec2Real descriptorSHS2body = *descriptorsNormalized["SHS2-2-body"];
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    Vec2Real& descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);
    Size      centralAtom = 0;
    std::for_each(descriptorsNormalized.begin(),
                  descriptorsNormalized.end(),
                  [&](Vec1Real& slice)
                  {
                      slice.resize(typeMap.countForceFieldTypes() * numberDescriptors);
                      for (Size neighborTypeStruc = 0;
                           neighborTypeStruc < typeMap.countStructureTypes();
                           neighborTypeStruc++)
                      {
                          Size neighborTypeFF = typeMap.toType(neighborTypeStruc);
                          Size neighborOffsetFF = neighborTypeFF * numberDescriptors;
                          Size neighborOffsetStruc = neighborTypeStruc * numberDescriptors;
                          for (Size desc = 0; desc < (Size)numberDescriptors; desc++)
                          {
                              slice[neighborOffsetFF + desc] =
                                  descriptorSHS2body[centralAtom][neighborOffsetStruc + desc];
                          }
                      }
                      centralAtom++;
                  });
}

void DescriptorCollector::allocateDescriptorSHS2BodyNormalized(const TypeMap& typeMap)
{
    const Vec1Int& nAtomsPerType = descriptors.at("SHS2-2-body")->get_nAtomsType();
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    Size centralType = 0;
    std::for_each(descriptorsNormalized.begin(),
                  descriptorsNormalized.begin() + typeMap.countStructureTypes(),
                  [&](Vec1Real& slice)
                  {
                      slice.resize(nAtomsPerType[centralType] * typeMap.countForceFieldTypes()
                                   * numberDescriptors);
                      std::fill(slice.begin(), slice.end(), (Real)0);
                      centralType++;
                  });
}

VASPML_EXEC_SPACE_SPECIFIER
const Descriptor& DescriptorCollector::getDescriptor(const String& key) const
{
    return (*descriptors.at(key));
}

void DescriptorCollector::writeDescriptorCollector() const
{
    for (const auto& [key, data] : descriptorsNormalized)
    {
        auto file = file_io::openFileO("DescriptorCollector_" + key + ".dat");
        for (const auto& x : *data)
        {
            for (const auto& y : x) file << str("%24.16E ", y) << std::endl;
        }
        file.close();
    }
}

const NearestNeighborNSquare& DescriptorCollector::get_neighborList(const String& key) const
{
    return (*descriptors.at(key)).get_neighborList();
}

const std::map<String, std::shared_ptr<Descriptor>>& DescriptorCollector::get_descriptors() const
{
    return descriptors;
}

const Real* DescriptorCollector::get_descriptorNormalizedPtrTO(const Int     type,
                                                               const Int     atom,
                                                               const String& key) const
{
    VASPML_DEBUG_L1(
        if (this->storageOrder != DescriptorStorage::Type)
        {
            global::tutor.bug("ERROR: DescriptorCollector::get_descriptorNormalizedPtrTO(): can "
                              "only be called if storageOrder is TypeOrder");
        }
        if (type >= (Int)descriptorsNormalized.at(key)->size())
        {
            global::tutor.bug("ERROR: DescriptorCollector::get_descriptorNormalizedPtrTO(): type "
                              "index out of bounds");
        }
        if (atom >= numberAtoms)
        {
            global::tutor.bug("ERROR: DescriptorCollector::get_descriptorNormalizedPtrTO(): atom "
                              "index out of bounds");
        }
        if (descriptorsNormalized.find(key) != descriptorsNormalized.cend())
        {
            global::tutor.bug("ERROR: DescriptorCollector::get_descriptorNormalizedPtrTO(): can "
                              "not find descriptor with key "
                              + key);
        }
    );
    Size descSize = descriptors.at(key)->get_sizeDescriptor(0);
    return &((*descriptorsNormalized.at(key))[type][atom * descSize]);
}

std::shared_ptr<const NearestNeighborNSquare> DescriptorCollector::getMaxNeighborList() const
{
    std::map<String, std::shared_ptr<const NearestNeighborNSquare>> neighborMap;
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptors.at(key)->get_weight() <= 0) continue;
        neighborMap[key] = descriptors.at(key)->get_neighborList_ptr();
    }
    return neighbor_list::getMaxNeighborListPtr(neighborMap);
}

void DescriptorCollector::computeBracketTerms(
    const Int                                            kAtom,
    const std::map<String, std::shared_ptr<Descriptor>>& descriptors,
    const NearestNeighborNSquare&                        maxNeighborList)
{
    this->descriptors = descriptors;
    compute_innerBracketSlowDerivativeHead(kAtom, maxNeighborList);
    compute_innerBracketSlowDerivative(kAtom, maxNeighborList);
    compute_outerBracketSlowDerivative();
}

void DescriptorCollector::compute_innerBracketSlowDerivative(
    const Int                     kAtom,
    const NearestNeighborNSquare& maxNeighborList)
{
    allocate_innerBracketSlowDerivative(maxNeighborList.get_size(kAtom));
    // insert code here
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptors[key]->get_weight() <= 0) continue;
    }
}

void DescriptorCollector::compute_innerBracketSlowDerivativeHead(
    const Int                     kAtom,
    const NearestNeighborNSquare& maxNeighborList)
{
    allocate_innerBracketSlowDerivativeHead(maxNeighborList.get_numberTypes());
    const Int kType = maxNeighborList.get_typeIndexCentral(kAtom);
    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptors[key]->get_weight() <= 0.0) continue;
        const Vec2Real& headTerm = descriptors[key]->get_refitHeadTerm();
        const Size      size = descriptors[key]->get_sizeDescriptor(kAtom);
        const Real*     head0 = &headTerm[kType][0];
        const Real*     headX = &headTerm[kType][size];
        const Real*     headY = &headTerm[kType][2 * size];
        const Real*     headZ = &headTerm[kType][3 * size];
        innerBracketSlowDerivativeHead[kType][0] += linalg::dotProduct(head0, headX, size);
        innerBracketSlowDerivativeHead[kType][1] += linalg::dotProduct(head0, headY, size);
        innerBracketSlowDerivativeHead[kType][2] += linalg::dotProduct(head0, headZ, size);
    }
}

void DescriptorCollector::compute_outerBracketSlowDerivative()
{
    //allocate_outerBracketSlowDerivative( maxNeighborList.get_size( kAtom ) );
    allocate_outerBracketSlowDerivative();
}

void DescriptorCollector::allocate_innerBracketSlowDerivative(const Int numberNeighbors)
{
    innerBracketSlowDerivative.resize(3 * numberNeighbors);
    std::fill(innerBracketSlowDerivative.begin(), innerBracketSlowDerivative.end(), (Real)0);
}

void DescriptorCollector::allocate_innerBracketSlowDerivativeHead(const Int numberTypes)
{
    innerBracketSlowDerivativeHead.resize(numberTypes);
    std::for_each(innerBracketSlowDerivativeHead.begin(),
                  innerBracketSlowDerivativeHead.end(),
                  [](Vec1Real& slice)
                  {
                      slice.resize(3);
                      std::fill(slice.begin(), slice.end(), (Real)0);
                  });
}

void DescriptorCollector::allocate_outerBracketSlowDerivative()
{}
