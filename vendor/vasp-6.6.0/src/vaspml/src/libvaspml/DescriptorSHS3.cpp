#include "DescriptorSHS3.hpp"

#include "DescriptorSHS2.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm>  // for for_each, fill, max_element, tran...
#include <cmath>      // for sqrt
#include <fstream>    // for basic_ofstream, basic_ios, basic_...
#include <functional> // for multiplies, bind, _1, _Placeholder
#include <stdexcept>  // for runtime_error

using namespace vaspml;

MapString vaspml::dataMapDescriptorSHS3()
{
    MapString m = {
        {"spectrumVasp", "Vec2Real"},
        {"spectrumPair", "Vec2Real"}
    };

    for (const auto& dataMapEntry : dataMapDescriptor()) m.insert(dataMapEntry);

    return m;
}

DescriptorSHS3::DescriptorSHS3(const Vec2Int&       descriptorList,
                               const bool           angularFilterOn,
                               const Int            angularFilterType,
                               const Real           angularFilterScaling,
                               const Vec1Size&      nRoots,
                               const Int            maxOrder,
                               const Real           weight,
                               const bool           isNormalized,
                               const DescriptorType descriptorType,
                               ShRec                record) :
    Descriptor(weight,
               isNormalized,
               descriptorType,
               assignOrMakeRecord(record, dataMapDescriptorSHS3())),
    spectrumVasp(data->get<Vec2Real>("spectrumVasp")),
    spectrumPair(data->get<Vec2Real>("spectrumPair"))
{
    this->angularFilterOn = angularFilterOn;
    this->angularFilterType = angularFilterType;
    this->angularFilterScaling = angularFilterScaling;
    descriptorSHS2 = nullptr;
    isSparseforcePreContractReady = false;
    numberElementsStructure = 0;
    numberElementsMLFF = descriptorList.size();
    make_angularFilter(maxOrder);
    make_sparseList(descriptorList, nRoots, maxOrder);
}

DescriptorSHS3::DescriptorSHS3(const Size&          nTypes,
                               const bool           angularFilterOn,
                               const Int            angularFilterType,
                               const Real           angularFilterScaling,
                               const Vec1Size&      nRoots,
                               const Int            maxOrder,
                               const Real           weight,
                               const bool           isNormalized,
                               const DescriptorType descriptorType,
                               ShRec                record) :
    Descriptor(weight,
               isNormalized,
               descriptorType,
               assignOrMakeRecord(record, dataMapDescriptorSHS3())),
    spectrumVasp(data->get<Vec2Real>("spectrumVasp")),
    spectrumPair(data->get<Vec2Real>("spectrumPair"))
{
    this->angularFilterOn = angularFilterOn;
    this->angularFilterType = angularFilterType;
    this->angularFilterScaling = angularFilterScaling;
    descriptorSHS2 = nullptr;
    isSparseforcePreContractReady = false;
    numberElementsStructure = 0;
    numberElementsMLFF = nTypes;
    make_angularFilter(maxOrder);
    make_sparseListRefit(nTypes, nRoots, maxOrder);
}

DescriptorSHS3::DescriptorSHS3(ShRec record) :
    Descriptor(assignOrMakeRecord(record, dataMapDescriptorSHS3())),
    spectrumVasp(data->get<Vec2Real>("spectrumVasp")),
    spectrumPair(data->get<Vec2Real>("spectrumPair"))
{}

void DescriptorSHS3::computeSHS3(const std::shared_ptr<DescriptorSHS2>& descriptorSHS2,
                                 const TypeMap&                         typeMap)
{

    if (weight <= 0) return;
    set_parameters(descriptorSHS2, typeMap);
    // actual calculations
    if (global::parallel.off()) { computeSHS3CPU(typeMap); }
    else
    {
        //nvtxRangePush( "computeSHS3GPULookUp" );
        make_lookUpSpectrumVasp(typeMap);
        make_lookUpForcePreContract(typeMap);
        computeSHS3GPULookUp();
        //nvtxRangePop();
    }
}

void DescriptorSHS3::set_parameters(const std::shared_ptr<DescriptorSHS2>& descriptorSHS2,
                                    const TypeMap&                         typeMap)
{
    this->descriptorSHS2 = descriptorSHS2;
    set_neighborList(this->descriptorSHS2->get_neighborList_ptr());
    numberAtoms = descriptorSHS2->get_nAtoms();
    numberElementsStructure = typeMap.countStructureTypes();
    allocateArrays(typeMap);

    if (!isSparseforcePreContractReady) prepareSparseforcePreContract(typeMap);
}

void DescriptorSHS3::computeSHS3CPU(const TypeMap& /* typeMap */)
{
    for (Size atom = 0; atom < (Size)numberAtoms; atom++)
    {
        computeSHS3SingleAtom(*descriptorSHS2, atom);
    }
}

VASPML_EXEC_SPACE_SPECIFIER
void DescriptorSHS3::computeSHS3SingleAtom(const DescriptorSHS2& descriptorSHS2, const Size atom)
{

    const Int centralType = neighborList->get_typeIndexCentral(atom);

    const Int       centralType_ff = typeMapLoc.toType(centralType);
    const Vec1Real& clnm_vasp = descriptorSHS2.get_clnmVasp(atom);
    Vec1Real&       spectrumVasp = this->spectrumVasp[atom];

    Size nDesc = 0;
    std::for_each(type1List[centralType_ff].cbegin(),
                  type1List[centralType_ff].cend(),
                  [&](const Int& type1_ff) mutable
                  {
                      const Int type2_ff = type2List[centralType_ff][nDesc];
                      const Int type1 = typeMapLoc.toSubType(type1_ff);
                      const Int type2 = typeMapLoc.toSubType(type2_ff);
                      if (type1 >= 0 && type2 >= 0)
                      {
                          const Int order = angularList[centralType_ff][nDesc];
                          const Int n0 = n0List[centralType_ff][nDesc];
                          const Int n1 = n1List[centralType_ff][nDesc];
                          Real      descTemp = 0;
                          Size      indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size)0);
                          Size      indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size)0);
                          for (Size m = 0; m < (Size)2 * order + 1; m++)
                          {
                              descTemp += clnm_vasp[indx0] * clnm_vasp[indx1];
                              indx0++;
                              indx1++;
                          }
                          spectrumVasp[nDesc] = weightFactor[centralType_ff][nDesc] * descTemp;
                      }
                      nDesc++;
                  });
}

void DescriptorSHS3::make_sparseList(const Vec2Int&  descriptorList,
                                     const Vec1Size& nRoots,
                                     const Int       maxOrder)
{
    make_sparseListInference(descriptorList, nRoots, maxOrder);
}

void DescriptorSHS3::make_sparseListInference(const Vec2Int&  descriptorList,
                                              const Vec1Size& nRoots,
                                              const Int       maxOrder)
{
    type0List.resize(descriptorList.size());
    type1List.resize(descriptorList.size());
    type2List.resize(descriptorList.size());
    angularList.resize(descriptorList.size());
    n0List.resize(descriptorList.size());
    n1List.resize(descriptorList.size());
    weightFactor.resize(descriptorList.size());

    for (Size type0 = 0; type0 < descriptorList.size(); type0++)
    {
        Int descCounter = 0;
        for (Size type1 = 0; type1 < descriptorList.size(); type1++)
        {
            for (Size type2 = 0; type2 < descriptorList.size(); type2++)
            {
                for (Size orderL = 0; orderL < (Size)maxOrder + 1; orderL++)
                {
                    for (Size order0 = 0; order0 < nRoots[orderL]; order0++)
                    {
                        for (Size order1 = order0; order1 < nRoots[orderL]; order1++)
                        {
                            descCounter++;
                            // check if descriptor list is sparsified
                            for (Size nDesc = 0; nDesc < descriptorList[type0].size(); nDesc++)
                            {
                                if (descCounter == descriptorList[type0][nDesc])
                                {
                                    type0List[type0].push_back(type0);
                                    type1List[type0].push_back(type1);
                                    type2List[type0].push_back(type2);
                                    angularList[type0].push_back(orderL);
                                    n0List[type0].push_back(order0);
                                    n1List[type0].push_back(order1);
                                    if (order0 == order1)
                                    {
                                        weightFactor[type0].push_back((Real)1.0
                                                                      * angularFilter[orderL]);
                                    }
                                    else
                                    {
                                        weightFactor[type0].push_back(constants::SQRT2
                                                                      * angularFilter[orderL]);
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void DescriptorSHS3::make_sparseListRefit(const Size&     nTypesSup,
                                          const Vec1Size& nRoots,
                                          const Int       maxOrder)
{
    type0List.clear();
    type0List.resize(nTypesSup);
    type1List.clear();
    type1List.resize(nTypesSup);
    type2List.clear();
    type2List.resize(nTypesSup);
    angularList.clear();
    angularList.resize(nTypesSup);
    n0List.clear();
    n0List.resize(nTypesSup);
    n1List.clear();
    n1List.resize(nTypesSup);
    weightFactor.clear();
    weightFactor.resize(nTypesSup);
    data->add("descriptorList", "Vec2Int");
    Vec2Int& descriptorList = data->get<Vec2Int>("descriptorList");
    descriptorList.resize(nTypesSup);
    for (Size type0 = 0; type0 < nTypesSup; type0++)
    {
        Int descCounter = 0;
        for (Size type1 = 0; type1 < nTypesSup; type1++)
        {
            for (Size type2 = 0; type2 < nTypesSup; type2++)
            {
                for (Size orderL = 0; orderL < (Size)maxOrder + 1; orderL++)
                {
                    for (Size order0 = 0; order0 < nRoots[orderL]; order0++)
                    {
                        for (Size order1 = order0; order1 < nRoots[orderL]; order1++)
                        {
                            type0List[type0].push_back(type0);
                            type1List[type0].push_back(type1);
                            type2List[type0].push_back(type2);
                            angularList[type0].push_back(orderL);
                            n0List[type0].push_back(order0);
                            n1List[type0].push_back(order1);
                            if (order0 == order1)
                            {
                                weightFactor[type0].push_back((Real)1.0 * angularFilter[orderL]);
                            }
                            else
                            {
                                weightFactor[type0].push_back(constants::SQRT2
                                                              * angularFilter[orderL]);
                            }
                            descriptorList[type0].push_back(descCounter);
                            descCounter++;
                        }
                    }
                }
            }
        }
    }
}

void DescriptorSHS3::make_angularFilter(const Int maxOrder)
{
    angularFilter.resize(maxOrder + 1);
    Real factor = constants::PI * constants::PI * (Real)8.0;

    for (Size order = 0; order < (Size)maxOrder + 1; order++)
    {
        if (angularFilterOn)
        {
            if (angularFilterType == 1)
            {
                angularFilter[order] =
                    std::sqrt(factor / (Real)(2 * order + 1)) / std::sqrt((Real)(2 * order + 1));
            }
            else if (angularFilterType == 2)
            {
                Real factor2 =
                    angularFilterScaling * (Real)((order * (order + 1)) * (order * (order + 1)));
                factor2 = ((Real)1 + factor2);
                factor2 *= factor2;
                angularFilter[order] = std::sqrt(factor / ((Real)(2 * order + 1))) / factor2;
            }
        }
        else { angularFilter[order] = std::sqrt(factor / (Real)(2 * order + 1)); }
    }
}

void DescriptorSHS3::allocateArrays(const TypeMap& typeMap)
{
    if (global::parallel.off())
    {
        typeMapLoc = typeMap;
        const auto& nList = descriptorSHS2->get_neighborList();
        spectrumVasp.resize(numberAtoms);
        for (Size atom = 0; atom < (Size)numberAtoms; atom++)
        {
            const Int& type = nList.get_typeIndexCentral(atom);
            const Int  ff_type = typeMap.toType(type);
            spectrumVasp[atom].resize(type0List[ff_type].size());
        }
    }
    else { resizeArraysGPU(typeMap); }
}

const Vec1Real& DescriptorSHS3::get_SHS3Atom_vasp(const Int centralAtom) const
{
    return spectrumVasp[centralAtom];
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& DescriptorSHS3::get_descriptor(const Size centralAtom) const
{
    return spectrumVasp[centralAtom];
}

VASPML_EXEC_SPACE_SPECIFIER
Int DescriptorSHS3::get_sizeDescriptor(const Size centralAtom) const
{
    return spectrumVasp[centralAtom].size();
}

Size DescriptorSHS3::get_numberElementsMLFF() const
{
    return numberElementsMLFF;
}
Size DescriptorSHS3::get_numberElementsStructure() const
{
    return numberElementsStructure;
}
bool DescriptorSHS3::get_angularFilterOn() const
{
    return angularFilterOn;
}
Int DescriptorSHS3::get_angularFilterType() const
{
    return angularFilterType;
}
Real DescriptorSHS3::get_angularFilterScaling() const
{
    return angularFilterScaling;
}

const Vec2Real& DescriptorSHS3::get_spectrumVasp() const
{
    return spectrumVasp;
}
const Vec1Real& DescriptorSHS3::get_spectrumVasp(const Size atomIdx) const
{
    return spectrumVasp[atomIdx];
}

void DescriptorSHS3::rescale_descriptor(const Size centralAtom, const Real scaleFactor)
{
    std::transform( // par_unseq,
        spectrumVasp[centralAtom].begin(),
        spectrumVasp[centralAtom].end(),
        spectrumVasp[centralAtom].begin(),
        std::bind(std::multiplies<Real>(), std::placeholders::_1, scaleFactor));
}

void DescriptorSHS3::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                              Vec2Real&       forcePreContract,
                                              const TypeMap&  typeMap) const
{
    VASPML_DEBUG_L2(
        if (!isSparseforcePreContractReady)
        {
            throw std::runtime_error("ERROR: You are trying to use compute_forcePreContract(const "
                                     "ShVec2Real & derivativeMatrix, ShVec2Real& forcePreContract, "
                                     "const TypeMap& typeMap) const.\n"
                                     "Sparse index arrays were not computed yet. Make sure that "
                                     "prepareSparseforcePreContract( const TypeMap& typeMap ) is "
                                     "called before");
        }
    );
    if (global::parallel.off())
    {
        compute_forcePreContractSparseSingleCPU(derivativeMatrix, forcePreContract, typeMap);
    }
    else { compute_forcePreContractSparse_stdpar(derivativeMatrix, forcePreContract); }
}

void DescriptorSHS3::compute_forcePreContractSparseSingleCPU(const Vec2Real& derivativeMatrix,
                                                             Vec2Real&       forcePreContract,
                                                             const TypeMap&  typeMap) const
{
    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();
    Size           centralAtom = 0;
    std::for_each(
        forcePreContract.begin(),
        forcePreContract.begin() + numberAtoms,
        [&](Vec1Real& slice)
        {
            Size            typeStruc = typeIndexCentral[centralAtom];
            Size            typeStrucFF = typeMap.toType(typeStruc);
            Size            centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(centralAtom);
            Size            centralAtomPerType = centralAtomIndexPerType[centralAtom];
            Size            nCombi = 0;
            std::for_each(
                preFactorPreContract[typeStrucFF].cbegin(),
                preFactorPreContract[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    Size index_derivativeMatrix = sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                    Size shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                    Size lnmCombiIndex = sparseMap_central[typeStrucFF][nCombi];
                    slice[lnmCombiIndex] +=
                        preFac
                        * derivativeMatrix[typeStruc][centralAtomShift * centralAtomPerType
                                                      + index_derivativeMatrix]
                        * descriptor2Body[shs2Index];
                    nCombi++;
                });

            centralAtom++;
        });
}

Size DescriptorSHS3::get_forcePreContractSize(Size atomIndex) const
{
    return descriptorSHS2->get_forcePreContractSize(atomIndex);
}

void DescriptorSHS3::computeForceTerms(const Vec2Real& forcePreContract,
                                       Vec2Real&       pairForces,
                                       Vec1Real&       centralForces) const
{
    descriptorSHS2->computeForceTerms(forcePreContract, pairForces, centralForces);
}

void DescriptorSHS3::prepareSparseforcePreContract(const TypeMap& typeMap)
{
    sparseMap_derivativeMatrix.resize(n0List.size());
    preFactorPreContract.resize(n0List.size());
    sparseMap_SHS2.resize(n0List.size());
    sparseMap_central.resize(n0List.size());
    // loop over central atom type
    for (Size centralTypeFF = 0; centralTypeFF < n0List.size(); centralTypeFF++)
    {
        // loop over central atom neighbor types. These are neighbor types as in the clnm
        for (Size centralNeighborTypeFF = 0; centralNeighborTypeFF < n0List.size();
             centralNeighborTypeFF++)
        {
            Int centralNeighborType = typeMap.toSubType(centralNeighborTypeFF);
            if (centralNeighborType < 0) continue;
            // loop over the l index of the central clnm.
            for (Size l0 = 0; l0 < descriptorSHS2->get_maxOrder() + 1; l0++)
            {
                // loop over the the radial index of the clnm of the central atom
                for (Size nRadial0 = 0; nRadial0 < descriptorSHS2->get_nRoots()[l0]; nRadial0++)
                {
                    // loop over the m index of the clnm of the central atom.
                    // One has to filter prefactors and index maps with respect to all possible
                    // combinations of the central atom clnm
                    for (Size m0 = 0; m0 < 2 * l0 + 1; m0++)
                    {
                        Size nDesc0 =
                            descriptorSHS2->get_Index((Size)centralNeighborType, l0, nRadial0, m0);
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        // loop over all possible descriptor combinations
                        // p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}
                        for (Size nDesc1 = 0; nDesc1 < n0List[centralTypeFF].size(); nDesc1++)
                        {
                            Size neighborType1FF = type1List[centralTypeFF][nDesc1];
                            Size neighborType2FF = type2List[centralTypeFF][nDesc1];
                            if (l0 != angularList[centralTypeFF][nDesc1]) continue;
                            Int neighborType1 = typeMap.toSubType(neighborType1FF);
                            Int neighborType2 = typeMap.toSubType(neighborType2FF);
                            if (neighborType1 < 0 || neighborType2 < 0) continue;
                            Size nRadial1 = n0List[centralTypeFF][nDesc1];
                            Size nRadial2 = n1List[centralTypeFF][nDesc1];
                            Size nDesc2 =
                                descriptorSHS2->get_Index((Size)neighborType1, l0, nRadial1, m0);
                            Size nDesc3 =
                                descriptorSHS2->get_Index((Size)neighborType2, l0, nRadial2, m0);
                            // check if the derivative
                            // \frac{p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{c^{iJ_{3}}_{ln_{3}m}} is non
                            // zero.
                            if (nDesc0 == nDesc2 || nDesc0 == nDesc3)
                            {
                                sparseMap_central[centralTypeFF].push_back(nDesc0);
                                sparseMap_derivativeMatrix[centralTypeFF].push_back(nDesc1);
                                if (nRadial1 == nRadial2 && neighborType1 == neighborType2)
                                {
                                    preFactorPreContract[centralTypeFF].push_back(
                                        (Real)2.0 * angularFilter[l0]);
                                }
                                else
                                {
                                    if (nRadial1 == nRadial2 && neighborType1 != neighborType2)
                                    {
                                        preFactorPreContract[centralTypeFF].push_back(
                                            angularFilter[l0]);
                                    }
                                    else
                                    {
                                        preFactorPreContract[centralTypeFF].push_back(
                                            std::sqrt((Real)2.0) * angularFilter[l0]);
                                    }
                                }
                                // store the entries of the descriptorSHS2 which have to be
                                // multiplied by the derivativeMatrix and summed over. note always
                                // the entry which is not the one in target forcePreContract is
                                // stored.
                                if (nDesc0 == nDesc2)
                                {
                                    sparseMap_SHS2[centralTypeFF].push_back(nDesc3);
                                }
                                else { sparseMap_SHS2[centralTypeFF].push_back(nDesc2); }
                            }
                        }
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    }
                }
            }
        }
        descriptorSize.push_back(n0List[centralTypeFF].size());
    }
    isSparseforcePreContractReady = true;
}

void DescriptorSHS3::write_spectrumVasp(const String& fname) const
{
    file_io::writeRealVec2D(spectrumVasp, fname);
}

void DescriptorSHS3::setParameters(const Real           weight,
                                   const bool           isNormalized,
                                   const DescriptorType descriptorType,
                                   const Size           numberElementsMLFF,
                                   const Size           numberElementsStructure,
                                   const bool           angularFilterOn,
                                   const Int            angularFilterType,
                                   const Real           angularFilterScaling)
{
    this->weight = weight;
    this->isNormalized = isNormalized;
    this->descriptorType = descriptorType;
    this->numberElementsMLFF = numberElementsMLFF;
    this->numberElementsStructure = numberElementsStructure;
    this->angularFilterOn = angularFilterOn;
    this->angularFilterType = angularFilterType;
    this->angularFilterScaling = angularFilterScaling;
}

void DescriptorSHS3::setData(const Vec2Real& spectrumVasp)
{
    this->spectrumVasp = spectrumVasp;
    //    this->spectrumPair = spectrumPair;
}

void DescriptorSHS3::extendData(const Vec2Real& spectrumVasp)
{
    this->spectrumVasp.insert(this->spectrumVasp.end(), spectrumVasp.cbegin(), spectrumVasp.cend());
    //this->spectrumPair.insert( this->spectrumPair.end(),
    //                   spectrumPair.cbegin(), spectrumPair.cend());
}

void DescriptorSHS3::addElement(const Vec1Real& spectrumVasp)
{
    this->spectrumVasp.push_back(spectrumVasp);
    //    this->spectrumPair.push_back( spectrumPair );
}

const Vec2Int& DescriptorSHS3::get_descriptorList() const
{
    return data->get<Vec2Int>("descriptorList");
}

void DescriptorSHS3::compute_RefitDescDerivatives(const Int                     kAtom,
                                                  const NearestNeighborNSquare& maxNeighborList)
{
    compute_dPSHead(kAtom);
    compute_dPS(kAtom, maxNeighborList);
}

void DescriptorSHS3::compute_dPS(const Int kAtom, const NearestNeighborNSquare& maxNeighborList)
{
    allocate_dPS(maxNeighborList);
    const DescriptorSHS2& descriptorSHS2 = *(this->descriptorSHS2);
    Int                   kType = neighborList->get_typeIndexCentral(kAtom);

    const Vec1Int& typeIndex = maxNeighborList.get_typeIndex(kAtom);
    const Vec1Int& globalIndex = maxNeighborList.get_globalIndex(kAtom);

    //std::ofstream indexFile("Index.dat", std::ios::app);

    const Vec2Real& clnmVasp = descriptorSHS2.get_clnmVasp();
    //const Vec2Real& clnmVaspDerivative  =  descriptorSHS2.get_clnmDerivativeCentralVasp();
    //const Vec2Real& clnmPairDerivative  =  descriptorSHS2.get_clnmPairDerivative();
    // IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0))
    for (Size centralAtom = 0; centralAtom < maxNeighborList.get_size(kAtom); centralAtom++)
    {
        const Int centralGlobal = globalIndex[centralAtom];
        const Int centralType = typeIndex[centralAtom];

        //        indexFile << kType << "   " << centralType << std::endl;

        Int nDesc = 0;
        std::for_each(type1List[centralType].cbegin(),
                      type1List[centralType].cend(),
                      [&](const Int& type1_ff) mutable
                      {
                          const Int type2_ff = type2List[centralType][nDesc];
                          const Int type1 = typeMapLoc.toSubType(type1_ff);
                          const Int type2 = typeMapLoc.toSubType(type2_ff);
                          if (type1 >= 0 && type2 >= 0)
                          {
                              const Int order = angularList[centralType][nDesc];
                              const Int n0 = n0List[centralType][nDesc];
                              const Int n1 = n1List[centralType][nDesc];
                              Real      descTemp = 0;
                              Size      indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size)0);
                              Size      indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size)0);
                              for (Size m = 0; m < (Size)2 * order + 1; m++)
                              {
                                  descTemp += clnmVasp[centralGlobal][indx0]
                                            * clnmVasp[centralGlobal][indx1];
                                  indx0++;
                                  indx1++;
                              }
                              dPS[centralAtom][nDesc] = weightFactor[centralType][nDesc] * descTemp;
                          }
                          nDesc++;
                      });

        if (centralAtom < neighborList->get_size(kAtom))
        {
            for (Int xyz = 0; xyz < 3; xyz++)
            {
                const Int xyzOffset =
                    type1List[centralType].size() * xyz + type1List[centralType].size();
                nDesc = 0;
                //std::for_each(type1List[centralType].cbegin(),
                //       type1List[centralType].cend(),
                //       [&](const Int& type1_ff) mutable
                for (Size i = 0; i < type1List[centralType].size(); i++)
                {
                    const Int type1_ff = type1List[centralType][i];

                    const Int type2_ff = type2List[centralType][nDesc];
                    const Int type1 = typeMapLoc.toSubType(type1_ff);
                    const Int type2 = typeMapLoc.toSubType(type2_ff);

                    if (type1 >= 0 && type2 >= 0)
                    {
                        const Int order = angularList[centralType][nDesc];
                        //const Int n0 = n0List[centralType][nDesc];
                        const Int n1 = n1List[centralType][nDesc];

                        Real descTerm0 = 0;

                        //Size indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size)0);
                        Size indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size)0);

                        //Int dIndx0 =
                        //    descriptorSHS2.get_IndexDir(centralAtom, xyz, order, n0, (Size)0);
                        //Int dIndx1 =
                        //    descriptorSHS2.get_IndexDir(centralAtom, xyz, order, n1, (Size)0);
                        //
                        // clnmPairDerivative[kAtom 0..4][ nachbar1* x lnm , y lnm zlnm, nachbar 2,
                        // ]
                        //
                        // pseudo code
                        // neighborTyp
                        // if ( centralType == type1 )

                        for (Size m = 0; m < (Size)2 * order + 1; m++)
                        {
                            //if ( kType == type1 )
                            //{
                            //descTerm0 += clnmPairDerivative[kAtom][dIndx0] *
                            //clnmVasp[centralGlobal][indx1]; descTerm0 +=
                            //clnmPairDerivative[kAtom][dIndx0];
                            descTerm0 += clnmVasp[centralGlobal][indx1];
                            //}
                            if (kType == type2)
                            {
                                //descTerm0 += clnmPairDerivative[kAtom][dIndx1] *
                                //clnmVasp[centralGlobal][indx0]; descTerm0 +=
                                //clnmPairDerivative[kAtom][dIndx1]; descTerm0 +=
                                //clnmVasp[centralGlobal][indx0];
                            }
                            //indx0++;
                            indx1++;
                            //dIndx0++;
                            //dIndx1++;
                        }
                        dPS[centralAtom][xyzOffset + nDesc] =
                            weightFactor[centralType][nDesc] * descTerm0;
                    }
                    nDesc++;
                }
            }
        }
    }

    std::ofstream outputFile("dPS.dat", std::ios::app);
    for (Size centralAtom = 0; centralAtom < maxNeighborList.get_size(kAtom); centralAtom++)
    {
        if (centralAtom < neighborList->get_size(kAtom))
        {
            for (Size j = type1List[kType].size(); j < dPS[centralAtom].size(); j++)
            //for ( Int j = 0; j <  type1List[kType].size(); j++ )
            {
                outputFile << dPS[centralAtom][j] << std::endl;
            }
        }
    }
}

void DescriptorSHS3::compute_dPSHead(const Int kAtom)
{
    allocate_dPSHead();
    Int                   kType = neighborList->get_typeIndexCentral(kAtom);
    const DescriptorSHS2& descriptorSHS2 = *(this->descriptorSHS2);
    const Vec2Real&       clnmVasp = descriptorSHS2.get_clnmVasp();
    const Vec2Real&       clnmVaspDerivative = descriptorSHS2.get_clnmDerivativeCentralVasp();

    Int nDesc = 0;
    std::for_each(type1List[kType].cbegin(),
                  type1List[kType].cend(),
                  [&](const Int& type1_ff) mutable
                  {
                      const Int type2_ff = type2List[kType][nDesc];
                      const Int type1 = typeMapLoc.toSubType(type1_ff);
                      const Int type2 = typeMapLoc.toSubType(type2_ff);
                      if (type1 >= 0 && type2 >= 0)
                      {
                          const Int order = angularList[kType][nDesc];
                          const Int n0 = n0List[kType][nDesc];
                          const Int n1 = n1List[kType][nDesc];
                          Real      descTemp = 0;
                          Size      indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size)0);
                          Size      indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size)0);
                          for (Size m = 0; m < (Size)2 * order + 1; m++)
                          {
                              descTemp += clnmVasp[kAtom][indx0] * clnmVasp[kAtom][indx1];
                              indx0++;
                              indx1++;
                          }
                          dPSHead[kType][nDesc] = weightFactor[kType][nDesc] * descTemp;
                      }
                      nDesc++;
                  });

    for (Int xyz = 0; xyz < 3; xyz++)
    {
        Int xyzOffset = type1List[kType].size() * xyz + type1List[kType].size();
        nDesc = 0;
        std::for_each(
            type1List[kType].cbegin(),
            type1List[kType].cend(),
            [&](const Int& type1_ff) mutable
            {
                const Int type2_ff = type2List[kType][nDesc];
                const Int type1 = typeMapLoc.toSubType(type1_ff);
                const Int type2 = typeMapLoc.toSubType(type2_ff);

                if (type1 >= 0 && type2 >= 0)
                {
                    const Int order = angularList[kType][nDesc];
                    const Int n0 = n0List[kType][nDesc];
                    const Int n1 = n1List[kType][nDesc];

                    Real descTerm0 = 0;
                    Real descTerm1 = 0;

                    Size indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size)0);
                    Size indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size)0);

                    Size indx0Deri = descriptorSHS2.get_IndexDir(type1, xyz, order, n0, (Size)0);
                    Size indx1Deri = descriptorSHS2.get_IndexDir(type2, xyz, order, n1, (Size)0);
                    for (Size m = 0; m < (Size)2 * order + 1; m++)
                    {
                        descTerm0 += clnmVaspDerivative[kAtom][indx0Deri] * clnmVasp[kAtom][indx1];
                        descTerm1 += clnmVasp[kAtom][indx0] * clnmVaspDerivative[kAtom][indx1Deri];
                        indx0++;
                        indx1++;
                        indx0Deri++;
                        indx1Deri++;
                    }
                    dPSHead[kType][nDesc + xyzOffset] =
                        weightFactor[kType][nDesc] * (descTerm0 + descTerm1);
                }
                nDesc++;
            });
    }

    //std::ofstream outputFile( "dPSHead.dat", std::ios::app );
    ////std::ofstream outputFile( "DescSHS3.dat" );
    //for ( Int j = 0; j <  dPSHead[kType].size(); j++ )
    //for ( Int j = 0; j <  type1List[kType].size(); j++ )
    //{
    //    //std::cout << << dC00Head[kType][j]<< std::endl;
    //    outputFile << dPSHead[kType][j] << std::endl;
    //}
    //outputFile.close();
    //    //std::cout << std::endl;
    //    //outputFile << "EMPTY " << std::endl;
    //    //throw;
}

void DescriptorSHS3::allocate_dPS(const NearestNeighborNSquare& maxNeighborList)
{
    if (global::parallel.off())
    {
        const Vec1Int& numberNeighbors = maxNeighborList.get_size();
        Int maxNeighbors = *std::max_element(numberNeighbors.cbegin(), numberNeighbors.cend());
        dPS.resize(maxNeighbors);
        std::for_each(dPS.begin(),
                      dPS.end(),
                      [&](Vec1Real& slice)
                      {
                          //slice.resize( neighborList->get_numberTypes()*type0List[0].size() * 4 );
                          slice.resize(type0List[0].size() * 4);
                          std::fill(slice.begin(), slice.end(), 0);
                      });
    }
    else global::tutor.bug("ERROR: " + flf(VASPML_FLF) + " GPU version not implemented!");
}

void DescriptorSHS3::allocate_dPSHead()
{
    if (global::parallel.off())
    {
        dPSHead.resize(neighborList->get_numberTypes());
        std::for_each(dPSHead.begin(),
                      dPSHead.end(),
                      [&](Vec1Real& slice)
                      {
                          slice.resize(type0List[0].size() * 4);
                          // TODO init to zero
                          std::fill(slice.begin(), slice.end(), 1);
                      });
    }
    else global::tutor.bug("ERROR: " + flf(VASPML_FLF) + " GPU version not implemented!");
}

Int DescriptorSHS3::get_totalNumberBasisFunctions() const
{
    return descriptorSHS2->get_totalNumberBasisFunctions();
}

const Vec2Real& DescriptorSHS3::get_refitHeadTerm() const
{
    return dPSHead;
}

namespace vaspml
{

void copyNeighborList(const DescriptorSHS3& descA, DescriptorSHS3& descB)
{
    descB.set_neighborList(descA.get_neighborList_ptr());
}

void copyParameters(const DescriptorSHS3& descA, DescriptorSHS3& descB)
{
    descB.setParameters(descA.get_weight(),
                        descA.get_isNormalized(),
                        descA.get_descriptorType(),
                        descA.get_numberElementsMLFF(),
                        descA.get_numberElementsStructure(),
                        descA.get_angularFilterOn(),
                        descA.get_angularFilterType(),
                        descA.get_angularFilterScaling());
}

void copyData(const DescriptorSHS3& descA, DescriptorSHS3& descB)
{
    descB.setData(descA.get_spectrumVasp());
}

void extendData(const DescriptorSHS3& descA, DescriptorSHS3& descB)
{
    descB.extendData(descA.get_spectrumVasp());
}

void addElementData(const DescriptorSHS3& descA, DescriptorSHS3& descB, const Size atomIdx)
{
    descB.addElement(descA.get_spectrumVasp(atomIdx));
}

} //namespace vaspml
