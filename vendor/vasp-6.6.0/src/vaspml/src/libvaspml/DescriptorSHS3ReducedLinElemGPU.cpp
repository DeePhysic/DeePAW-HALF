#include "DescriptorSHS3ReducedLinElem.hpp" // IWYU pragma: associated

#include "DescriptorSHS2.hpp"
#include "Linalg.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"

#include <algorithm> // for for_each
#include <numeric>   // for iota

using namespace vaspml;

void DescriptorSHS3ReducedLinElem::computeSHS3GPU(const TypeMap& typeMap)
{
    Vec1Size centralAtomIndex(numberAtoms);
    std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);

    std::for_each( // par_unseq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + numberAtoms,
        [&](const Size& atom) { computeSHS3SingleAtom(typeMap, atom); });
}

void DescriptorSHS3ReducedLinElem::resizeArraysGPU(const TypeMap& typeMap)
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        auto nList = descriptorSHS2->get_neighborList();
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
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
        auto    nList = descriptorSHS2->get_neighborList();
        Vec1Int sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = descriptorSize[typeMap.toType(size)]; });
        resize = spectrumVaspSize.checkResize2Dim(sizeArray);
        if (resize) { spectrumVaspSize.resizeArray2Dim(spectrumVasp, sizeArray); }
    }
}

void DescriptorSHS3ReducedLinElem::computeSHS3clnmNoElementGPU()
{
    const Int basisSetSize = (Size)descriptorSHS2->get_basisSetSize();

    // One dimension is number of atoms,
    // other is descriptor without the types
    for (Size atom = 0; atom < (Size)numberAtoms; atom++)
    {
        // Get whole descriptor (l, n, m) for each atom
        const Vec1Real& cnlmVasp = descriptorSHS2->get_descriptor(atom);
        for (Size type = 0; type < numberElementsStructure; type++)
        {
            linalg::scaleVectorPlusVector((Real)1.0,
                                          &cnlmVasp[type * basisSetSize],
                                          clnmVaspNoElement[atom],
                                          basisSetSize);
        }
    }
}

void DescriptorSHS3ReducedLinElem::compute_forcePreContractSparseGPU(
    const Vec2Real& derivativeMatrix,
    Vec2Real&       forcePreContract,
    const TypeMap&  typeMap) const
{
    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();

    Vec1Size centralAtomIndex(numberAtoms);
    std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);

    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + numberAtoms,
        [this,
         derivativeMatrix,
         forcePreContract,
         typeMap,
         typeIndexCentral,
         centralAtomIndexPerType](const Size& derivativeAtom) mutable
        {
            Size            typeStruc = typeIndexCentral[derivativeAtom];
            Size            typeStrucFF = typeMap.toType(typeStruc);
            Size            centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(derivativeAtom);
            Size            centralAtomPerType = centralAtomIndexPerType[derivativeAtom];
            Size            centralProductTmp = centralAtomShift * centralAtomPerType;
            // first contraction
            Size nCombi = 0;
            std::for_each(
                preFactorPreContractB[typeStrucFF].cbegin(),
                preFactorPreContractB[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    // calculate sum over all element types for Cnlms
                    Real presumTemp = 0;
                    for (Size neighborType2 = 0; neighborType2 < numberElementsStructure;
                         neighborType2++)
                    {
                        presumTemp +=
                            descriptor2Body[sparseMap_SHS2_B[typeStrucFF][neighborType2][nCombi]];
                    }
                    Size index_derivativeMatrix =
                        centralProductTmp + sparseMap_derivativeMatrix_B[typeStrucFF][nCombi];
                    Size centralIndex = sparseMap_central_B[typeStrucFF][nCombi];
                    forcePreContract[derivativeAtom][centralIndex] +=
                        preFac * derivativeMatrix[typeStruc][index_derivativeMatrix] * presumTemp;
                    nCombi++;
                });
            // second contraction
            nCombi = 0;
            std::for_each(preFactorPreContract[typeStrucFF].cbegin(),
                          preFactorPreContract[typeStrucFF].cend(),
                          [&](const Real& preFac)
                          {
                              Size index_derivativeMatrix =
                                  centralProductTmp
                                  + sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                              Size centralIndex = sparseMap_central[typeStrucFF][nCombi];
                              Size shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                              forcePreContract[derivativeAtom][centralIndex] +=
                                  preFac * derivativeMatrix[typeStruc][index_derivativeMatrix]
                                  * descriptor2Body[shs2Index];
                              nCombi++;
                          });
        });
}
