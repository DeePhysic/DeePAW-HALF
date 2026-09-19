#include "DescriptorSHS2.hpp" // IWYU pragma: associated

#include "BasisFunctions.hpp"
#include "ParallelEnvironment.hpp"
#include "Timer.hpp"
#include "TypeMap.hpp"
#include "math.hpp"
#include "nearest_neighbor.hpp"

#include <algorithm> // for for_each, fill
#include <numeric>   // for iota

using namespace vaspml;

Size LookUpClnmVasp::countElements(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++)
            {
                for (Size m = 0; m < 2 * l + 1; m++) { counter++; }
            }
        }
    }
    return counter;
}

void LookUpClnmVasp::update(const Size&      nAtoms,
                            const Size&      lmax,
                            const Vec1Size&  nRoots,
                            ArrayResizing1D& lookUpSize)
{
    Size nElements = countElements(nAtoms, lmax, nRoots);
    bool resize = lookUpSize.checkResize(nElements);
    if (resize)
    {
        mainIndex.resize(nElements);
        centralAtom.resize(nElements);
        ll.resize(nElements);
        nRoot.resize(nElements);
        mm.resize(nElements);
        refill(nAtoms, lmax, nRoots);
    }
}

void LookUpClnmVasp::refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots)
{

    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++)
            {
                for (Size m = 0; m < 2 * l + 1; m++)
                {
                    mainIndex[counter] = counter;
                    centralAtom[counter] = atom;
                    ll[counter] = l;
                    nRoot[counter] = n;
                    mm[counter] = m;
                    counter++;
                }
            }
        }
    }
}

Size ForcePreContractLookUp::countElements(const Size& nAtoms,
                                           const Size& nTypeStruc,
                                           const Size& basisSetSize)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size neighborType = 0; neighborType < nTypeStruc; neighborType++)
        {
            for (Size nDesc = 0; nDesc < basisSetSize; nDesc++) { counter++; }
        }
    }
    return counter;
}

void ForcePreContractLookUp::update(const Size&           nAtoms,
                                    const Size&           nTypesStruc,
                                    const Size&           basisSetSize,
                                    const TypeMap&        typeMap,
                                    const Vec1Int&        typeIndexCentral,
                                    const Vec1Int&        atomIndexInType,
                                    const DescriptorSHS2& descriptorSHS2)
{
    if (nAtoms != nAtomsOld)
    {
        nAtomsOld = nAtoms;
        Size nElements = countElements(nAtoms, nTypesStruc, basisSetSize);
        mainIndex.resize(nElements);
        centralAtom.resize(nElements);
        centralType.resize(nElements);
        atomIndexPerType.resize(nElements);
        neighborTypeStruc.resize(nElements);
        matrixIndex.resize(nElements);
        forcePreIndx.resize(nElements);
        refill(nAtoms,
               nTypesStruc,
               basisSetSize,
               typeMap,
               typeIndexCentral,
               atomIndexInType,
               descriptorSHS2);
    }
}

void ForcePreContractLookUp::refill(const Size&           nAtoms,
                                    const Size&           nTypesStruc,
                                    const Size&           basisSetSize,
                                    const TypeMap&        typeMap,
                                    const Vec1Int&        typeIndexCentral,
                                    const Vec1Int&        atomIndexInType,
                                    const DescriptorSHS2& descriptorSHS2)
{
    Size atomShift = basisSetSize * typeMap.countForceFieldTypes();
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size neighborType = 0; neighborType < nTypesStruc; neighborType++)
        {
            Size neighborTypeFF = typeMap.toType(neighborType);
            for (Size nDesc = 0; nDesc < basisSetSize; nDesc++)
            {
                mainIndex[counter] = counter;
                centralAtom[counter] = atom;
                centralType[counter] = typeIndexCentral[atom];
                atomIndexPerType[counter] = atomIndexInType[atom];
                neighborTypeStruc[counter] = neighborType;
                matrixIndex[counter] =
                    atomIndexInType[atom] * atomShift + neighborTypeFF * basisSetSize + nDesc;
                forcePreIndx[counter] = descriptorSHS2.get_Index(neighborType, 0, nDesc, 0);
                counter++;
            }
        }
    }
}

Size LookUpRadialSplineNN::countElements(const Size&     nAtoms,
                                         const Size&     lmax,
                                         const Vec1Size& nRoots)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++) { counter++; }
        }
    }
    return counter;
}

void LookUpRadialSplineNN::update(const Size&      nAtoms,
                                  const Size&      lmax,
                                  const Vec1Size&  nRoots,
                                  ArrayResizing1D& lookUpSize)
{
    bool resize = (nAtoms != nAtomsOld);
    if (resize)
    {
        Size nElements = countElements(nAtoms, lmax, nRoots);
        lookUpSize.actDim = nElements;
        mainIndex.resize(nElements);
        centralAtom.resize(nElements);
        ll.resize(nElements);
        nRoot.resize(nElements);
        combiIndex.resize(nElements);
        refill(nAtoms, lmax, nRoots);
    }
    nAtomsOld = nAtoms;
}

void LookUpRadialSplineNN::refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size combi = 0;
        offset = 0;
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++)
            {
                mainIndex[counter] = counter;
                centralAtom[counter] = atom;
                ll[counter] = l;
                nRoot[counter] = n;
                combiIndex[counter] = combi;
                counter++;
                combi++;
                offset++;
            }
        }
    }
}

Size LookUpYlmNN::countElements(const Size& nAtoms, const Size& lmax)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size l = 0; l < lmax; l++)
        {
            for (Size m = 0; m < 2 * l + 1; m++) { counter++; }
        }
    }
    return counter;
}

void LookUpYlmNN::update(const Size& nAtoms, const Size& lmax, ArrayResizing1D& lookUpSize)
{
    bool resize = (nAtoms != nAtomsOld);
    if (resize)
    {
        Size nElements = countElements(nAtoms, lmax);
        lookUpSize.actDim = nElements;
        mainIndex.resize(nElements);
        centralAtom.resize(nElements);
        ll.resize(nElements);
        mm.resize(nElements);
        ylmEntry.resize(nElements);
        ylmdEntry.resize(nElements);
        refill(nAtoms, lmax);
    }
    nAtomsOld = nAtoms;
}

void LookUpYlmNN::refill(const Size& nAtoms, const Size& lmax)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size ylmIndex = 0;
        Size ylmdIndex = 0;
        offset = 0;
        for (Size l = 0; l < lmax; l++)
        {
            for (Size m = 0; m < 2 * l + 1; m++)
            {
                mainIndex[counter] = counter;
                centralAtom[counter] = atom;
                ll[counter] = l;
                mm[counter] = m;
                ylmEntry[counter] = ylmIndex;
                ylmdEntry[counter] = ylmdIndex;
                counter++;
                ylmIndex++;
                ylmdIndex += 3;
                offset++;
            }
        }
    }
}

Size LookUpclnmPairNN::countElements(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++)
            {
                for (Size m = 0; m < 2 * l + 1; m++) { counter++; }
            }
        }
    }
    return counter;
}

void LookUpclnmPairNN::update(const Size&      nAtoms,
                              const Size&      lmax,
                              const Vec1Size&  nRoots,
                              ArrayResizing1D& lookUpSize)
{
    bool resize = (nAtomsOld != nAtoms);
    if (resize)
    {
        Size nElements = countElements(nAtoms, lmax, nRoots);
        lookUpSize.actDim = nElements;
        mainIndex.resize(nElements);
        centralAtom.resize(nElements);
        ll.resize(nElements);
        nRoot.resize(nElements);
        mm.resize(nElements);
        clnmPairEntry.resize(nElements);
        clnmPairDOffset.resize(nElements);
        radialEntry.resize(nElements);
        ylmEntry.resize(nElements);
        ylmdEntry.resize(nElements);
        refill(nAtoms, lmax, nRoots);
    }
    nAtomsOld = nAtoms;
}

void LookUpclnmPairNN::refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size nbasis = 0;
        Size besselOffset = 0;
        Size ylm_index_offset = 0;
        Size ylm_index_xyz_offset = 0;
        Size OFFSET = 0;
        offsetRadial = 0;
        offsetYLM = 0;
        offsetCLNM = 0;
        Size nbasisD = 0;
        for (Size l = 0; l < lmax; l++)
        {
            for (Size n = 0; n < nRoots[l]; n++)
            {
                Size ylm_index = ylm_index_offset;
                Size ylm_index_xyz = ylm_index_xyz_offset;
                for (Size m = 0; m < 2 * l + 1; m++)
                {
                    mainIndex[counter] = counter;
                    centralAtom[counter] = atom;
                    ll[counter] = l;
                    nRoot[counter] = n;
                    mm[counter] = m;
                    clnmPairEntry[counter] = nbasis;
                    clnmPairDOffset[counter] = OFFSET + nbasisD;
                    radialEntry[counter] = besselOffset + n;
                    ylmEntry[counter] = ylm_index;
                    ylmdEntry[counter] = ylm_index_xyz;
                    counter++;
                    nbasis++;
                    nbasisD += 1;
                    ylm_index++;
                    ylm_index_xyz += 3;
                    offsetCLNM++;
                }
            }
            ylm_index_offset += 2 * l + 1;
            ylm_index_xyz_offset += 6 * l + 3;
            besselOffset += nRoots[l];
            offsetRadial += nRoots[l];
            offsetYLM += 2 * l + 1;
        }
        OFFSET += 3 * nbasisD;
    }
}

void LookUpForceTerms::update(const Size& nAtoms)
{
    bool resize = (nAtomsOld != nAtoms);
    if (resize)
    {
        Size nElements = 3 * nAtoms;
        mainIndex.resize(nElements);
        centralAtomIndex.resize(nElements);
        xyzIndex.resize(nElements);
        refill(nAtoms);
    }
}

void LookUpForceTerms::refill(const Size& nAtoms)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size xyz = 0; xyz < 3; xyz++)
        {
            mainIndex[counter] = counter;
            centralAtomIndex[counter] = atom;
            xyzIndex[counter] = xyz;
            counter++;
        }
    }
}

void DescriptorSHS2::resizePairArraysGPU(const Size& ldim)
{
    const Size&    nAtoms = neighborList->get_numberAtoms();
    const Vec1Int& numberNeighbors = neighborList->get_size();
    bool           resize = clnmPairSize.checkResize1Dim(nAtoms);
    // number of atoms increased. Reset all arrays set new max Size
    if (resize)
    {
        clnmPairSize.resizeArray1Dim(clnmPair, numberNeighbors, basisSetSize);
        clnmPairDerivativeSize.act1Dim = nAtoms;
        clnmPairDerivativeSize.max1Dim = nAtoms;
        clnmPairDerivativeSize.resizeArray1Dim(clnmPairDerivative,
                                               numberNeighbors,
                                               3 * basisSetSize);

        clnmPairDerivativeSize.resizeArray1Dim(clnmPairDerivativeNO,
                                               numberNeighbors,
                                               3 * basisSetSize);
        ylmSize.act1Dim = nAtoms;
        ylmSize.max1Dim = nAtoms;
        ylmSize.resizeArray1Dim(ylm, numberNeighbors, ldim);
        ylmDerivativeSize.act1Dim = nAtoms;
        ylmDerivativeSize.max1Dim = nAtoms;
        ylmDerivativeSize.resizeArray1Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        centralAtomIndexSize.checkResize(nAtoms);
        centralAtomIndex.resize(nAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }
    else
    {
        // set new actual sizes of first dimension.
        clnmPairDerivativeSize.act1Dim = nAtoms;
        ylmSize.act1Dim = nAtoms;
        ylmDerivativeSize.act1Dim = nAtoms;
        // check second dimension which means number of neighbors changed or
        // basis set size; Check already valid for all arrays
        resize = clnmPairSize.checkResize2Dim(numberNeighbors, basisSetSize);
        if (resize)
        {
            clnmPairSize.resizeArray2Dim(clnmPair, numberNeighbors, basisSetSize);
            clnmPairDerivativeSize.resizeArray2Dim(clnmPairDerivative,
                                                   numberNeighbors,
                                                   3 * basisSetSize);
            clnmPairDerivativeSize.resizeArray2Dim(clnmPairDerivativeNO,
                                                   numberNeighbors,
                                                   3 * basisSetSize);
            ylmSize.resizeArray2Dim(ylm, numberNeighbors, ldim);
            ylmDerivativeSize.resizeArray2Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        }
    }
}

void DescriptorSHS2::resizeArraysVaspFormatGPU()
{
    const Size  nAtoms = neighborList->get_nAtoms();
    const Size& numberTypes = neighborList->get_numberTypes();
    bool        resize = clnmVaspSize.checkResize1Dim(nAtoms);
    // number of atoms increased. Reset all arrays set new max Size
    if (resize)
    {
        clnmVaspSize.resizeArray1Dim(clnmVasp, basisSetSize * numberTypes);
        clnmDerivativeCentralVaspSize.act1Dim = nAtoms;
        clnmDerivativeCentralVaspSize.max1Dim = nAtoms;
        clnmDerivativeCentralVaspSize.resizeArray1Dim(clnmDerivativeCentralVasp,
                                                      3 * basisSetSize * numberTypes);
    }
    else
    {
        clnmDerivativeCentralVaspSize.act1Dim = nAtoms;
        resize = clnmVaspSize.checkResize2Dim(basisSetSize * numberTypes);
        if (resize)
        {

            clnmVaspSize.resizeArray2Dim(clnmVasp, basisSetSize * numberTypes);
            clnmDerivativeCentralVaspSize.resizeArray1Dim(clnmDerivativeCentralVasp,
                                                          3 * basisSetSize * numberTypes);
        }
    }
    // send data back to GPU if it was resized. And rezero
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                clnmVasp.begin(),
                clnmVasp.begin() + clnmVaspSize.act1Dim,
                [&](Vec1Real& slice)
                { std::fill(slice.begin(), slice.begin() + clnmVaspSize.actSizeRec, (Real)0); });
        },
        __func__);
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          clnmDerivativeCentralVasp.begin(),
                          clnmDerivativeCentralVasp.begin() + clnmDerivativeCentralVaspSize.act1Dim,
                          [&](Vec1Real& slice)
                          {
                              std::fill(slice.begin(),
                                        slice.begin() + clnmDerivativeCentralVaspSize.actSizeRec,
                                        (Real)0);
                          });
        },
        __func__);
}

void DescriptorSHS2::resize_radialBasisFunction()
{
    bool resize = radialBasisFunctionSize.checkResize1Dim(nAtoms);
    if (resize)
    {
        radialBasisFunctionSize.act1Dim = nAtoms;
        radialBasisFunctionSize.max1Dim = nAtoms;
        const Vec1Int& numberNeighbors = neighborList->get_size();
        radialBasisFunctionSize.resizeArray1Dim(radialBasisFunction, numberNeighbors, nRadialBasis);
        radialBasisFunctionSize.resizeArray1Dim(radialBasisFunctionDerivative,
                                                numberNeighbors,
                                                nRadialBasis);
    }
    else
    {
        const Vec1Int& numberNeighbors = neighborList->get_size();
        resize = radialBasisFunctionSize.checkResize2Dim(numberNeighbors, nRadialBasis);
        if (resize)
        {
            radialBasisFunctionSize.resizeArray2Dim(radialBasisFunction,
                                                    numberNeighbors,
                                                    nRadialBasis);
            radialBasisFunctionSize.resizeArray2Dim(radialBasisFunctionDerivative,
                                                    numberNeighbors,
                                                    nRadialBasis);
        }
    }
}

void DescriptorSHS2::resize_typeIndexArrayGPU(const Size& numberTypes)
{
    bool resize = typeIndexArraySize.checkResize(numberTypes);
    if (resize)
    {
        typeIndexArray.resize(numberTypes);
        std::iota(typeIndexArray.begin(), typeIndexArray.end(), 0);
    }
}

VASPML_EXEC_SPACE_SPECIFIER
void DescriptorSHS2::updateSinglePairCoefficients_stdpar(const Size& atom,
                                                         const Size& radIndx,
                                                         const Size& ylmIndx,
                                                         const Size& ylmDIndx,
                                                         const Real& normed_x,
                                                         const Real& normed_y,
                                                         const Real& normed_z,
                                                         Real&       clnmPair,
                                                         Real&       clnmPairDerivativeX,
                                                         Real&       clnmPairDerivativeY,
                                                         Real&       clnmPairDerivativeZ)
{
    clnmPair = radialBasisFunction[atom][radIndx] * ylm[atom][ylmIndx];
    clnmPairDerivativeX =
        ylm[atom][ylmIndx] * radialBasisFunctionDerivative[atom][radIndx] * normed_x
        + ylmDerivative[atom][ylmDIndx] * radialBasisFunction[atom][radIndx];

    clnmPairDerivativeY =
        ylm[atom][ylmIndx] * radialBasisFunctionDerivative[atom][radIndx] * normed_y
        + ylmDerivative[atom][ylmDIndx + 1] * radialBasisFunction[atom][radIndx];

    clnmPairDerivativeZ =
        ylm[atom][ylmIndx] * radialBasisFunctionDerivative[atom][radIndx] * normed_z
        + ylmDerivative[atom][ylmDIndx + 2] * radialBasisFunction[atom][radIndx];
}

void DescriptorSHS2::arrayResizing_stdpar()
{
    // this should be removed
    nRoots = basisFunctions->get_nValidRoots();
    // plus 1 to be able to write loops with < smbol
    maxOrder = basisFunctions->get_maxOrder() + 1;
    basisSetSize = compute_BasisSetSize(maxOrder, nRoots);
    // allocate pair coefficient arrays
    nRadialBasis = basisFunctions->get_totalNumberBesselFunctions();
    computeOffsets();
    // allocate arrays for radial Basis function
    nAtoms = neighborList->get_numberAtoms();
    allocatePairCoefficientArrays(basisFunctions->get_ldim());
    resize_radialBasisFunction();
    allocateArraysVaspFormat();
    resize_centralAtomIndex();
    resize_radialBasisFunction();
    resize_ylm(basisFunctions->get_ldim());
    allocate_typeIndexArray(neighborList->get_numberTypes());
    lookUpClnmVasp.update(nAtoms, maxOrder, nRoots, lookUpClnmVaspSize);
}

void DescriptorSHS2::resize_LookUpTablesForceCompute(const TypeMap& typeMap)
{
    resize_forcePreContractLookUp(typeMap);
    resize_lookUpForceTerms();
}

void DescriptorSHS2::resize_forcePreContractLookUp(const TypeMap& typeMap)
{
    if (!global::parallel.off())
    {
        forcePreContractLookUp.update(nAtoms,
                                      typeMap.countStructureTypes(),
                                      basisSetSize,
                                      typeMap,
                                      neighborList->get_typeIndexCentral(),
                                      neighborList->get_centralAtomIndexPerType(),
                                      *this);
    }
}

void DescriptorSHS2::resize_lookUpForceTerms()
{
    if (!global::parallel.off()) lookUpForceTerms.update(nAtoms);
}

void DescriptorSHS2::compute_forcePreContract_stdpar(const Vec2Real& derivativeMatrix,
                                                     Vec2Real&       forcePreContract) const
{
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          forcePreContractLookUp.mainIndex.cbegin(),
                          forcePreContractLookUp.mainIndex.cend(),
                          [&](const Size& indx)
                          {
                              const Size& centralAtom = forcePreContractLookUp.centralAtom[indx];
                              const Size& forcePreIndx = forcePreContractLookUp.forcePreIndx[indx];
                              const Size& centralType = forcePreContractLookUp.centralType[indx];
                              const Size& matrixIndex = forcePreContractLookUp.matrixIndex[indx];
                              forcePreContract[centralAtom][forcePreIndx] =
                                  derivativeMatrix[centralType][matrixIndex];
                          });
        },
        __func__);
}

void DescriptorSHS2::updatePairCoefficients_stdparNN()
{
    const Vec2Real& distance = neighborList->get_distances();
    const Vec2Real& connectVectorNormalized = neighborList->get_connectionVectorNormalized();

    //nvtxRangePush( "updatePairCoefficientsGPU" );
    VASPML_PROFILING_START("SPLINES");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpRadialSplineNN.mainIndex.begin(),
                lookUpRadialSplineNN.mainIndex.begin() + lookUpRadialSplineNNSize.actDim,
                [&](const Size& indx)
                {
                    const Size&     centralAtom = lookUpRadialSplineNN.centralAtom[indx];
                    const Size&     ll = lookUpRadialSplineNN.ll[indx];
                    const Size&     nRoot = lookUpRadialSplineNN.nRoot[indx];
                    const Size&     combiIndex = lookUpRadialSplineNN.combiIndex[indx];
                    const Vec1Real& distance = neighborList->get_distances(centralAtom);
                    for (Size neighborAtom = 0; neighborAtom < distance.size(); neighborAtom++)
                    {
                        Size entry = lookUpRadialSplineNN.offset * neighborAtom + combiIndex;
                        basisFunctions->interpolate(
                            distance[neighborAtom],
                            ll,
                            nRoot,
                            radialBasisFunction[centralAtom][entry],
                            radialBasisFunctionDerivative[centralAtom][entry]);
                    }
                });
        },
        __func__);
    VASPML_PROFILING_STOP("SPLINES");

    VASPML_PROFILING_START("YLMS");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          lookUpYlmNN.mainIndex.begin(),
                          lookUpYlmNN.mainIndex.begin() + lookUpYlmNNSize.actDim,
                          [&](const Size& indx)
                          {
                              const Size& centralAtom = lookUpYlmNN.centralAtom[indx];
                              const Size& ll = lookUpYlmNN.ll[indx];
                              const Size& mm = lookUpYlmNN.mm[indx];
                              Int         m = mm - ll;
                              for (Size neighborAtom = 0;
                                   neighborAtom < distance[centralAtom].size();
                                   neighborAtom++)
                              {
                                  Size entryYLM = lookUpYlmNN.ylmEntry[indx]
                                                + neighborAtom * lookUpYlmNN.offset;
                                  Size entryYLMd = lookUpYlmNN.ylmdEntry[indx]
                                                 + 3 * neighborAtom * lookUpYlmNN.offset;
                                  basisFunctions->computeAngularBasis(
                                      connectVectorNormalized[centralAtom][3 * neighborAtom],
                                      connectVectorNormalized[centralAtom][3 * neighborAtom + 1],
                                      connectVectorNormalized[centralAtom][3 * neighborAtom + 2],
                                      distance[centralAtom][neighborAtom],
                                      ll,
                                      m,
                                      ylm[centralAtom][entryYLM],
                                      ylmDerivative[centralAtom][entryYLMd],
                                      ylmDerivative[centralAtom][entryYLMd + 1],
                                      ylmDerivative[centralAtom][entryYLMd + 2]);
                              }
                          });
        },
        __func__);
    VASPML_PROFILING_STOP("YLMS");

    VASPML_PROFILING_START("CLNMS");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpclnmPairNN.mainIndex.begin(),
                lookUpclnmPairNN.mainIndex.begin() + lookUpclnmPairNNSize.actDim,
                [&](const Size& indx)
                {
                    const Size& centralAtom = lookUpclnmPairNN.centralAtom[indx];
                    for (Size neighborAtom = 0; neighborAtom < distance[centralAtom].size();
                         neighborAtom++)
                    {
                        const Size& radialEntry = lookUpclnmPairNN.radialEntry[indx]
                                                + neighborAtom * lookUpclnmPairNN.offsetRadial;
                        const Size connectionOffset = 3 * neighborAtom;
                        const Size clnmPairEntry = lookUpclnmPairNN.clnmPairEntry[indx]
                                                 + neighborAtom * lookUpclnmPairNN.offsetCLNM;
                        const Size clnmPairDOffset = lookUpclnmPairNN.clnmPairDOffset[indx]
                                                   + 3 * neighborAtom * lookUpclnmPairNN.offsetCLNM;
                        const Size ylmEntry = lookUpclnmPairNN.ylmEntry[indx]
                                            + neighborAtom * lookUpclnmPairNN.offsetYLM;
                        const Size ylmdEntry = lookUpclnmPairNN.ylmdEntry[indx]
                                             + 3 * neighborAtom * lookUpclnmPairNN.offsetYLM;
                        updateSinglePairCoefficients_stdpar(
                            centralAtom,
                            radialEntry,
                            ylmEntry,
                            ylmdEntry,
                            connectVectorNormalized[centralAtom][connectionOffset],
                            connectVectorNormalized[centralAtom][connectionOffset + 1],
                            connectVectorNormalized[centralAtom][connectionOffset + 2],
                            clnmPair[centralAtom][clnmPairEntry],
                            clnmPairDerivative[centralAtom][clnmPairDOffset],
                            clnmPairDerivative[centralAtom][clnmPairDOffset + basisSetSize],
                            clnmPairDerivative[centralAtom][clnmPairDOffset + 2 * basisSetSize]);
                    }
                });
        },
        __func__);
    VASPML_PROFILING_STOP("CLNMS");
    //nvtxRangePop();
}

void DescriptorSHS2::updatePairCoefficients_stdparNO()
{
    const Vec2Real& distance = neighborList->get_distances();
    const Vec2Real& connectVectorNormalized = neighborList->get_connectionVectorNormalized();

    //nvtxRangePush( "updatePairCoefficientsGPU" );
    VASPML_PROFILING_START("SPLINES");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          lookUpRadialSplineNN.mainIndex.begin(),
                          lookUpRadialSplineNN.mainIndex.begin() + lookUpRadialSplineNNSize.actDim,
                          [&](const Size& indx)
                          {
                              const Size& centralAtom = lookUpRadialSplineNN.centralAtom[indx];
                              const Size& ll = lookUpRadialSplineNN.ll[indx];
                              const Size& nRoot = lookUpRadialSplineNN.nRoot[indx];
                              const Size& combiIndex = lookUpRadialSplineNN.combiIndex[indx];
                              Size        entry = combiIndex * distance[centralAtom].size();
                              for (Size neighborAtom = 0;
                                   neighborAtom < distance[centralAtom].size();
                                   neighborAtom++)
                              {
                                  basisFunctions->interpolate(
                                      distance[centralAtom][neighborAtom],
                                      ll,
                                      nRoot,
                                      radialBasisFunction[centralAtom][entry],
                                      radialBasisFunctionDerivative[centralAtom][entry]);
                                  entry++;
                              }
                          });
        },
        __func__);
    VASPML_PROFILING_STOP("SPLINES");
    VASPML_PROFILING_START("YLMS");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpYlmNN.mainIndex.begin(),
                lookUpYlmNN.mainIndex.begin() + lookUpYlmNNSize.actDim,
                [&](const Size& indx)
                {
                    const Size& centralAtom = lookUpYlmNN.centralAtom[indx];
                    const Size& ll = lookUpYlmNN.ll[indx];
                    const Size& mm = lookUpYlmNN.mm[indx];
                    Int         m = mm - ll;
                    Size entryYLM = lookUpYlmNN.ylmEntry[indx] * distance[centralAtom].size();
                    Size entryYLMd = lookUpYlmNN.ylmdEntry[indx] * distance[centralAtom].size();
                    for (Size neighborAtom = 0; neighborAtom < distance[centralAtom].size();
                         neighborAtom++)
                    {
                        basisFunctions->computeAngularBasis(
                            connectVectorNormalized[centralAtom][3 * neighborAtom],
                            connectVectorNormalized[centralAtom][3 * neighborAtom + 1],
                            connectVectorNormalized[centralAtom][3 * neighborAtom + 2],
                            distance[centralAtom][neighborAtom],
                            ll,
                            m,
                            ylm[centralAtom][entryYLM],
                            ylmDerivative[centralAtom][entryYLMd],
                            ylmDerivative[centralAtom][entryYLMd + 1],
                            ylmDerivative[centralAtom][entryYLMd + 2]);
                        entryYLM++;
                        entryYLMd += 3;
                    }
                });
        },
        __func__);
    VASPML_PROFILING_STOP("YLMS");
    VASPML_PROFILING_START("CLNMS");
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpclnmPairNN.mainIndex.begin(),
                lookUpclnmPairNN.mainIndex.begin() + lookUpclnmPairNNSize.actDim,
                [&](const Size& indx)
                {
                    const Size& centralAtom = lookUpclnmPairNN.centralAtom[indx];
                    // cache friendlyness just turn entry around
                    Size radialEntry =
                        lookUpclnmPairNN.radialEntry[indx] * distance[centralAtom].size();
                    Size clnmPairEntry =
                        lookUpclnmPairNN.clnmPairEntry[indx] * distance[centralAtom].size();
                    Size clnmPairDOffset =
                        3 * lookUpclnmPairNN.clnmPairDOffset[indx] * distance[centralAtom].size();
                    Size ylmEntry = lookUpclnmPairNN.ylmEntry[indx] * distance[centralAtom].size();
                    Size ylmdEntry =
                        3 * lookUpclnmPairNN.ylmEntry[indx] * distance[centralAtom].size();
                    Size connectionOffset = 0;
                    for (Size neighborAtom = 0; neighborAtom < distance[centralAtom].size();
                         neighborAtom++)
                    {
                        updateSinglePairCoefficients_stdpar(
                            centralAtom,
                            radialEntry,
                            ylmEntry,
                            ylmdEntry,
                            connectVectorNormalized[centralAtom][connectionOffset],
                            connectVectorNormalized[centralAtom][connectionOffset + 1],
                            connectVectorNormalized[centralAtom][connectionOffset + 2],
                            clnmPair[centralAtom][clnmPairEntry],
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffset],
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffset + 1],
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffset + 2]);
                        radialEntry++;
                        clnmPairEntry++;
                        clnmPairDOffset += 3;
                        ylmEntry++;
                        ylmdEntry += 3;
                        connectionOffset += 3;
                    }
                });
        },
        __func__);
    VASPML_PROFILING_STOP("CLNMS");
    //nvtxRangePop();
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficientsNO_stdpar()
{
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpClnmVasp.mainIndex.cbegin(),
                lookUpClnmVasp.mainIndex.cbegin() + lookUpClnmVaspSize.actDim,
                [&](const Size& index) mutable
                {
                    Size indexAtom = lookUpClnmVasp.centralAtom[index];
                    Size ll = lookUpClnmVasp.ll[index];
                    Size mm = lookUpClnmVasp.mm[index];
                    Size nRoot = lookUpClnmVasp.nRoot[index];

                    auto& clnmPair = this->clnmPair[indexAtom];
                    auto& clnmPairDerivativeNO = this->clnmPairDerivativeNO[indexAtom];
                    auto& clnmVasp = this->clnmVasp[indexAtom];
                    auto& clnmDerivativeCentralVasp = this->clnmDerivativeCentralVasp[indexAtom];
                    const Vec1Int& neighborType = neighborList->get_typeIndex(indexAtom);
                    Size           clnmIndex = lOffset[ll] + nRoot * (2 * ll + 1) + mm;

                    Size clnmPairEntry = neighborList->get_size(indexAtom) * clnmIndex;
                    Size clnmPairDEntry = 3 * neighborList->get_size(indexAtom) * clnmIndex;

                    for (Size neighbor = 0; neighbor < neighborList->get_size(indexAtom);
                         neighbor++)
                    {
                        Size typeOffset = neighborType[neighbor] * basisSetSize;
                        Size typeOffset3x = 3 * neighborType[neighbor] * basisSetSize;

                        clnmVasp[typeOffset + clnmIndex] += clnmPair[clnmPairEntry];

                        clnmDerivativeCentralVasp[typeOffset3x + clnmIndex] -=
                            clnmPairDerivativeNO[clnmPairDEntry];

                        clnmDerivativeCentralVasp[typeOffset3x + basisSetSize + clnmIndex] -=
                            clnmPairDerivativeNO[clnmPairDEntry + 1];

                        clnmDerivativeCentralVasp[typeOffset3x + 2 * basisSetSize + clnmIndex] -=
                            clnmPairDerivativeNO[clnmPairDEntry + 2];
                        clnmPairEntry++;
                        clnmPairDEntry += 3;
                    }
                });
        },
        __func__);
    vaspFormatComputed = true;
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficientsNN_stdpar()
{
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpClnmVasp.mainIndex.cbegin(),
                lookUpClnmVasp.mainIndex.cbegin() + lookUpClnmVaspSize.actDim,
                [&](const Size& index) mutable
                {
                    Size indexAtom = lookUpClnmVasp.centralAtom[index];
                    Size ll = lookUpClnmVasp.ll[index];
                    Size mm = lookUpClnmVasp.mm[index];
                    Size nRoot = lookUpClnmVasp.nRoot[index];

                    auto& clnmPair = this->clnmPair[indexAtom];
                    auto& clnmPairDerivative = this->clnmPairDerivative[indexAtom];
                    auto& clnmVasp = this->clnmVasp[indexAtom];
                    auto& clnmDerivativeCentralVasp = this->clnmDerivativeCentralVasp[indexAtom];

                    const Vec1Int& neighborType = neighborList->get_typeIndex(indexAtom);
                    Size           clnmIndex = lOffset[ll] + nRoot * (2 * ll + 1) + mm;

                    for (Size neighbor = 0; neighbor < neighborList->get_size(indexAtom);
                         neighbor++)
                    {
                        Size typeOffset = neighborType[neighbor] * basisSetSize;
                        Size typeOffset3x = 3 * neighborType[neighbor] * basisSetSize;

                        clnmVasp[typeOffset + clnmIndex] +=
                            clnmPair[neighbor * basisSetSize + clnmIndex];

                        clnmDerivativeCentralVasp[typeOffset3x + clnmIndex] -=
                            clnmPairDerivative[3 * neighbor * basisSetSize + clnmIndex];

                        clnmDerivativeCentralVasp[typeOffset3x + basisSetSize + clnmIndex] -=
                            clnmPairDerivative[3 * neighbor * basisSetSize + basisSetSize
                                               + clnmIndex];

                        clnmDerivativeCentralVasp[typeOffset3x + 2 * basisSetSize + clnmIndex] -=
                            clnmPairDerivative[3 * neighbor * basisSetSize + 2 * basisSetSize
                                               + clnmIndex];
                    }
                });
        },
        __func__);
    vaspFormatComputed = true;
}

void DescriptorSHS2::computeForceTerms_stdpar(const Vec2Real& forcePreContract,
                                              Vec2Real&       pairForces,
                                              Vec1Real&       centralForces) const
{
    const Vec2Int& neighborTypes = neighborList->get_typeIndex();

    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpForceTerms.mainIndex.cbegin(),
                lookUpForceTerms.mainIndex.cend(),
                [&](const Size& indx)
                {
                    const Size& atom = lookUpForceTerms.centralAtomIndex[indx];
                    const Size& xyz = lookUpForceTerms.xyzIndex[indx];
                    std::for_each(
                        typeIndexArray.cbegin(),
                        typeIndexArray.cbegin() + typeIndexArraySize.actDim,
                        [&](const Size& neighborType)
                        {
                            const Real* forcePreContractPtr =
                                &forcePreContract[atom][basisSetSize * neighborType];
                            const Real* shecPtr =
                                &clnmDerivativeCentralVasp[atom][basisSetSize * neighborType * 3
                                                                 + xyz * basisSetSize];
                            centralForces[3 * atom + xyz] -=
                                math::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                        });
                    // TODO: Rewrite this without lambda capture or find out why compiler is
                    // complaining.
                    Size neighborAtom = 0;
                    std::for_each(neighborTypes[atom].cbegin(),
                                  neighborTypes[atom].cend(),
                                  [&](const Int& neighborType)
                                  {
                                      const Real* forcePreContractPtr =
                                          &forcePreContract[atom][basisSetSize * neighborType];
                                      const Real* shecPtr =
                                          &clnmPairDerivative[atom][basisSetSize * neighborAtom * 3
                                                                    + xyz * basisSetSize];
                                      Real dotProduct = math::dotProduct(shecPtr,
                                                                         forcePreContractPtr,
                                                                         basisSetSize);
                                      pairForces[atom][3 * neighborAtom + xyz] = -dotProduct;
                                      neighborAtom++;
                                  });
                });
        },
        __func__);
}

void DescriptorSHS2::reorder_clnmPairDerivative_stdpar()
{
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                lookUpclnmPairNN.mainIndex.begin(),
                lookUpclnmPairNN.mainIndex.begin() + lookUpclnmPairNNSize.actDim,
                [&](const Size& indx)
                {
                    const Size& centralAtom = lookUpclnmPairNN.centralAtom[indx];
                    Size        clnmPairDOffsetNO = 3 * lookUpclnmPairNN.clnmPairDOffset[indx]
                                           * neighborList->get_size(centralAtom);
                    for (Size neighborAtom = 0; neighborAtom < neighborList->get_size(centralAtom);
                         neighborAtom++)
                    {
                        const Size clnmPairDOffset = lookUpclnmPairNN.clnmPairDOffset[indx]
                                                   + 3 * neighborAtom * lookUpclnmPairNN.offsetCLNM;
                        clnmPairDerivative[centralAtom][clnmPairDOffset] =
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffsetNO];
                        clnmPairDerivative[centralAtom][clnmPairDOffset + basisSetSize] =
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffsetNO + 1];
                        clnmPairDerivative[centralAtom][clnmPairDOffset + 2 * basisSetSize] =
                            clnmPairDerivativeNO[centralAtom][clnmPairDOffsetNO + 2];

                        clnmPairDOffsetNO += 3;
                    }
                });
        },
        __func__);
}
