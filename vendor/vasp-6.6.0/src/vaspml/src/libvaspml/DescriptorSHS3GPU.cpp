#include "DescriptorSHS3.hpp" // IWYU pragma: associated

#include "DescriptorSHS2.hpp"
#include "ParallelEnvironment.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm> // for for_each

using namespace vaspml;

LookUpSparseMatrixspectrumVasp::LookUpSparseMatrixspectrumVasp()
{
    nAtomsOld = 0;
}

Size LookUpSparseMatrixspectrumVasp::countElements(const Size&     nAtoms,
                                                   const Vec2Size& descriptorList,
                                                   const TypeMap&  typeMap,
                                                   const Vec1Int&  typeIndexCentral)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Int type = typeMap.toType(typeIndexCentral[atom]);
        for (Size n = 0; n < descriptorList[type].size(); n++) { counter++; }
    }
    return counter;
}

void LookUpSparseMatrixspectrumVasp::update(const Size&     nAtoms,
                                            const Vec2Size& descriptorList,
                                            const TypeMap&  typeMap,
                                            const Vec1Int&  typeIndexCentral)
{
    // TODO here I should in principle compare the number of atoms
    // per type

    if (nAtoms != nAtomsOld)
    {
        nAtomsOld = nAtoms;
        Size nElements = countElements(nAtoms, descriptorList, typeMap, typeIndexCentral);
        mainIndex.resize(nElements);
        typeCentralFF.resize(nElements);
        nDesc.resize(nElements);
        centralIndex.resize(nElements);
        refill(nAtoms, descriptorList, typeMap, typeIndexCentral);
    }
}

void LookUpSparseMatrixspectrumVasp::refill(const Size&     nAtoms,
                                            const Vec2Size& descriptorList,
                                            const TypeMap&  typeMap,
                                            const Vec1Int&  typeIndexCentral)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Int type = typeMap.toType(typeIndexCentral[atom]);
        for (Size n = 0; n < descriptorList[type].size(); n++)
        {
            mainIndex[counter] = counter;
            typeCentralFF[counter] = type;
            nDesc[counter] = n;
            centralIndex[counter] = atom;
            counter++;
        }
    }
}

void LookUpForcePreContract::init(const Vec2Real& preFactorPreContract,
                                  const Vec2Size& sparseMap_derivativeMatrix,
                                  const Vec2Size& sparseMap_SHS2,
                                  const Vec2Size& sparseMap_central,
                                  const TypeMap&  typeMap)
{
    if (initialized) return;
    shs2Index.resize(typeMap.countForceFieldTypes());
    indexDerivativeMatrix.resize(typeMap.countForceFieldTypes());
    unique_lnmCombiIndex.resize(typeMap.countForceFieldTypes());
    shs2Start.resize(typeMap.countForceFieldTypes());
    shs2End.resize(typeMap.countForceFieldTypes());
    preFactor.resize(typeMap.countForceFieldTypes());
    for (Size typeFF = 0; typeFF < typeMap.countForceFieldTypes(); typeFF++)
    {
        for (Size nDesc = 0; nDesc < preFactorPreContract[typeFF].size(); nDesc++)
        {
            shs2Index[typeFF].push_back(sparseMap_SHS2[typeFF][nDesc]);
            indexDerivativeMatrix[typeFF].push_back(sparseMap_derivativeMatrix[typeFF][nDesc]);
            preFactor[typeFF].push_back(preFactorPreContract[typeFF][nDesc]);
        }

        unique_lnmCombiIndex[typeFF] = vector_tools::get_unique(sparseMap_central[typeFF]);
        shs2Start[typeFF].resize(unique_lnmCombiIndex[typeFF].size());
        shs2End[typeFF].resize(unique_lnmCombiIndex[typeFF].size());
        for (Size lnm = 0; lnm < unique_lnmCombiIndex[typeFF].size(); lnm++)
        {
            Size counter = 0;
            for (Size nDesc = 0; nDesc < preFactorPreContract[typeFF].size(); nDesc++)
            {
                if (unique_lnmCombiIndex[typeFF][lnm] == sparseMap_central[typeFF][nDesc]
                    and counter == 0)
                {
                    shs2Start[typeFF][lnm] = nDesc;
                    counter++;
                }
                if (unique_lnmCombiIndex[typeFF][lnm] < sparseMap_central[typeFF][nDesc])
                {
                    shs2End[typeFF][lnm] = nDesc;
                    break;
                }
                if (nDesc == preFactorPreContract[typeFF].size() - 1)
                {
                    shs2End[typeFF][lnm] = nDesc + 1;
                }
                if (unique_lnmCombiIndex[typeFF][lnm] == sparseMap_central[typeFF][nDesc])
                    counter++;
            }
        }
    }
    initialized = true;
}

Size LookUpForcePreContract::countElements(const Size&    nAtoms,
                                           const Vec1Int& typeIndexCentral,
                                           const TypeMap& typeMap)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size typeFF = typeMap.toType(typeIndexCentral[atom]);
        for (Size nDesc = 0; nDesc < unique_lnmCombiIndex[typeFF].size(); nDesc++) { counter++; }
    }
    return counter;
}

void LookUpForcePreContract::update(const Size&    nAtoms,
                                    const Vec1Int& typeIndexCentralIn,
                                    const Vec1Int& centralAtomIndexPerTypeIn,
                                    const TypeMap& typeMap)
{
    if (nAtoms != nAtomsOld)
    {
        nAtomsOld = nAtoms;
        Size nElements = countElements(nAtoms, typeIndexCentralIn, typeMap);
        mainIndex.resize(nElements);
        centralIndex.resize(nElements);
        centralIndexPerType.resize(nElements);
        typeIndexCentral.resize(nElements);
        typeIndexCentralFF.resize(nElements);
        lnmCombiCentral.resize(nElements);
        refill(nAtoms, typeIndexCentralIn, centralAtomIndexPerTypeIn, typeMap);
    }
}

void LookUpForcePreContract::refill(const Size&    nAtoms,
                                    const Vec1Int& typeIndexCentralIn,
                                    const Vec1Int& centralAtomIndexPerTypeIn,
                                    const TypeMap& typeMap)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size type = typeIndexCentralIn[atom];
        Size typeFF = typeMap.toType(type);
        for (Size ndesc = 0; ndesc < unique_lnmCombiIndex[typeFF].size(); ndesc++)
        {
            mainIndex[counter] = counter;
            centralIndex[counter] = atom;
            typeIndexCentral[counter] = type;
            typeIndexCentralFF[counter] = typeFF;
            lnmCombiCentral[counter] = ndesc;
            centralIndexPerType[counter] = centralAtomIndexPerTypeIn[atom];
            counter++;
        }
    }
}

void DescriptorSHS3::resizeArraysGPU(const TypeMap& typeMap)
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        const auto& nList = descriptorSHS2->get_neighborList();
        spectrumVaspSize.act1Dim = numberAtoms;
        spectrumVaspSize.max1Dim = numberAtoms;
        Vec1Int sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = type0List[typeMap.toType(size)].size(); });
        spectrumVaspSize.resizeArray1Dim(spectrumVasp, sizeArray);
        typeMapLoc = typeMap;
    }
    else
    {
        spectrumVaspSize.act1Dim = numberAtoms;
        const auto& nList = descriptorSHS2->get_neighborList();
        Vec1Int     sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = descriptorSize[typeMap.toType(size)]; });
        resize = spectrumVaspSize.checkResize2Dim(sizeArray);
        if (resize) { spectrumVaspSize.resizeArray2Dim(spectrumVasp, sizeArray); }
    }
}

void DescriptorSHS3::make_lookUpTables(const TypeMap& typeMap)
{
    make_lookUpSpectrumVasp(typeMap);
    make_lookUpForcePreContract(typeMap);
}

void DescriptorSHS3::make_lookUpSpectrumVasp(const TypeMap& typeMap)
{
    // TODO HERE ONE SHOULD CHECK THE number of atoms per type.
    lookUpSpectrumVasp.update(numberAtoms, n0List, typeMap, neighborList->get_typeIndexCentral());
}

void DescriptorSHS3::make_lookUpForcePreContract(const TypeMap& typeMap)
{
    lookUpForcePreContract.init(preFactorPreContract,
                                sparseMap_derivativeMatrix,
                                sparseMap_SHS2,
                                sparseMap_central,
                                typeMap);
    lookUpForcePreContract.update(numberAtoms,
                                  neighborList->get_typeIndexCentral(),
                                  neighborList->get_centralAtomIndexPerType(),
                                  typeMap);
}

VASPML_EXEC_SPACE_SPECIFIER
void DescriptorSHS3::computeSHS3SingleAtomLookUp(const Size&     centralType_ff,
                                                 const Size&     nDesc,
                                                 const Vec1Real& clnmVasp,
                                                 Real&           spectrumVasp)
{
    const Int type1_ff = type1List[centralType_ff][nDesc];
    const Int type2_ff = type2List[centralType_ff][nDesc];
    const Int type1 = typeMapLoc.toSubType(type1_ff);
    const Int type2 = typeMapLoc.toSubType(type2_ff);
    if (type1 >= 0 && type2 >= 0)
    {
        const Int order = angularList[centralType_ff][nDesc];
        const Int n0 = n0List[centralType_ff][nDesc];
        const Int n1 = n1List[centralType_ff][nDesc];
        Real      descTemp = 0;
        Size      indx0 = descriptorSHS2->get_Index(type1, order, n0, (Size)0);
        Size      indx1 = descriptorSHS2->get_Index(type2, order, n1, (Size)0);

        for (Size m = 0; m < (Size)2 * order + 1; m++)
        {
            descTemp += clnmVasp[indx0] * clnmVasp[indx1];
            indx0++;
            indx1++;
        }
        spectrumVasp = weightFactor[centralType_ff][nDesc] * descTemp;
    }
}

void DescriptorSHS3::computeSHS3GPULookUp()
{
    const Vec2Real& clnmVasp = descriptorSHS2->get_descriptor();
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          lookUpSpectrumVasp.mainIndex.cbegin(),
                          lookUpSpectrumVasp.mainIndex.cend(),
                          [&](const Size& indx)
                          {
                              computeSHS3SingleAtomLookUp(
                                  lookUpSpectrumVasp.typeCentralFF[indx],
                                  lookUpSpectrumVasp.nDesc[indx],
                                  clnmVasp[lookUpSpectrumVasp.centralIndex[indx]],
                                  spectrumVasp[lookUpSpectrumVasp.centralIndex[indx]]
                                              [lookUpSpectrumVasp.nDesc[indx]]);
                          });
        },
        __func__);
}

void DescriptorSHS3::compute_forcePreContractSparse_stdpar(const Vec2Real& derivativeMatrix,
                                                           Vec2Real&       forcePreContract) const
{
    global::parallel.run(
        [&]([[maybe_unused]] auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpForcePreContract.mainIndex.cbegin(),
                lookUpForcePreContract.mainIndex.cend(),
                [&](const Size& indx)
                {
                    const Size&     centralAtom = lookUpForcePreContract.centralIndex[indx];
                    const Size&     type = lookUpForcePreContract.typeIndexCentral[indx];
                    const Size&     typeFF = lookUpForcePreContract.typeIndexCentralFF[indx];
                    const Size&     lnmCombiCentral = lookUpForcePreContract.lnmCombiCentral[indx];
                    const Size&     centralAtomShift = n0List[typeFF].size();
                    const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(centralAtom);
                    const Size&     centralAtomPerType =
                        lookUpForcePreContract.centralIndexPerType[indx];

                    for (Size contInd = lookUpForcePreContract.shs2Start[typeFF][lnmCombiCentral];
                         contInd < lookUpForcePreContract.shs2End[typeFF][lnmCombiCentral];
                         contInd++)
                    {
                        const Size& index_derivativeMatrix =
                            sparseMap_derivativeMatrix[typeFF][contInd];
                        const Size& shs2Index = sparseMap_SHS2[typeFF][contInd];
                        const Real& preFac = preFactorPreContract[typeFF][contInd];
                        forcePreContract[centralAtom][lnmCombiCentral] +=
                            preFac
                            * derivativeMatrix[type][centralAtomShift * centralAtomPerType
                                                     + index_derivativeMatrix]
                            * descriptor2Body[shs2Index];
                    }
                });
        },
        __func__);
}
