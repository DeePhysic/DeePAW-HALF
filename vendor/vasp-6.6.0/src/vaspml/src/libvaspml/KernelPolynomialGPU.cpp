#include "KernelPolynomial.hpp" // IWYU pragma: associated

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "ParallelEnvironment.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "math.hpp"
#include "nearest_neighbor.hpp"

#include <algorithm>

using namespace vaspml;

/*================================================================================================+
 |
 | SECTION Lookup Table stuff first.
 |
 +================================================================================================*/

Size KernelPolynomialLookUp::countElements(const Size&    nAtoms,
                                           const Vec1Int& typeCentral,
                                           const TypeMap& typeMap,
                                           const Vec1Int& nLocRefConfTypes)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size type = typeMap.toType(typeCentral[atom]);
        for (Size mIndx = 0; mIndx < (Size)nLocRefConfTypes[type]; mIndx++) { counter++; }
    }

    return counter;
}

void KernelPolynomialLookUp::update(const Size&      nAtoms,
                                    const Vec1Int&   typeCentral,
                                    const TypeMap&   typeMap,
                                    const Vec1Int&   nLocRefConfTypes,
                                    const Vec1Int&   centralAtomIndexPerType,
                                    ArrayResizing1D& kernelPolynomialLookUpSize)
{
    bool resize = kernelPolynomialLookUpSize.checkResize(nAtoms);
    if (resize)
    {
        Size numberElements = countElements(nAtoms, typeCentral, typeMap, nLocRefConfTypes);

        mainIndex.resize(numberElements);
        centralAtom.resize(numberElements);
        typeStruc.resize(numberElements);
        typeFF.resize(numberElements);
        kernelEntry.resize(numberElements);
        refill(nAtoms, typeCentral, typeMap, nLocRefConfTypes, centralAtomIndexPerType);
    }
}

void KernelPolynomialLookUp::refill(const Size&    nAtoms,
                                    const Vec1Int& typeCentral,
                                    const TypeMap& typeMap,
                                    const Vec1Int& nLocRefConfTypes,
                                    const Vec1Int& centralAtomIndexPerType)
{
    Size counter = 0;
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        Size type = typeCentral[atom];
        Size typeFFLoc = typeMap.toType(type);
        Size numberLocalRefConfsType = nLocRefConfTypes[typeFFLoc];
        Size increaseCounter = centralAtomIndexPerType[atom] * numberLocalRefConfsType;
        for (Size mIndx = 0; mIndx < (Size)nLocRefConfTypes[typeFFLoc]; mIndx++)
        {
            mainIndex[counter] = counter;
            centralAtom[counter] = atom;
            typeStruc[counter] = type;
            typeFF[counter] = typeFFLoc;
            kernelEntry[counter] = increaseCounter + mIndx;
            counter++;
        }
    }
}

/*================================================================================================+
 |
 | SECTION Actual KernelPolynomial stuff starts here.
 |
 +================================================================================================*/

void KernelPolynomial::allocateDerivativeArraysGPU()
{
    const Size numberAtoms =
        descriptorCollection->getDescriptor("SHS2-2-body").get_neighborList().get_numberAtoms();

    bool resize = centralAtomIndexSize.checkResize(numberAtoms);

    resize = polynomialDerivativeSize.checkResize1Dim(nAtomsPerType.size());
    if (resize)
    {
        polynomialDerivativeSize.resizeArray1Dim(polynomialDerivative, kernelMatrixDim);
        polynomialDerivativeSize.resizeArray1Dim(polynomialKernelNorm, kernelMatrixDim);
        polynomialKernel = kernelMatrix;
    }
    else
    {
        resize = polynomialDerivativeSize.checkResize2Dim(kernelMatrixDim);
        if (resize)
        {
            polynomialDerivativeSize.resizeArray2Dim(polynomialDerivative, kernelMatrixDim);
            polynomialDerivativeSize.resizeArray2Dim(polynomialKernelNorm, kernelMatrixDim);
            polynomialKernel = kernelMatrix;
        }
    }
}

void KernelPolynomial::allocateDerivativeArraysRefitGPU()
{
    const Size numberAtoms =
        descriptorCollection->getDescriptor("SHS2-2-body").get_neighborList().get_numberAtoms();

    bool resize = centralAtomIndexSize.checkResize(numberAtoms);

    resize = polynomialDerivativeSize.checkResize1Dim(nAtomsPerType.size());
    if (resize)
    {
        polynomialDerivativeSize.resizeArray1Dim(polynomialDerivative, kernelMatrixDim);
        polynomialKernel = kernelMatrix;
    }
    else
    {
        resize = polynomialDerivativeSize.checkResize2Dim(kernelMatrixDim);
        if (resize)
        {
            polynomialDerivativeSize.resizeArray2Dim(polynomialDerivative, kernelMatrixDim);
            polynomialKernel = kernelMatrix;
        }
    }
}

void KernelPolynomial::makeLookUpTables(const TypeMap& typeMap)
{
    kernelPolynomialLookUp.update(
        descriptorCollection->getDescriptor("SHS2-2-body").get_neighborList().get_numberAtoms(),
        descriptorCollection->getDescriptor("SHS2-2-body").get_typeIndexCentral(),
        typeMap,
        (*this->numberLocalRefConfs),
        descriptorCollection->getDescriptor("SHS2-2-body")
            .get_neighborList()
            .get_centralAtomIndexPerType(),
        kernelPolynomialLookUpSize);
}

void KernelPolynomial::computePolynomialKernelTerms_stdlib()
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();

    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                kernelPolynomialLookUp.mainIndex.cbegin(),
                kernelPolynomialLookUp.mainIndex.cend(),
                [&](const Size& indx)
                {
                    Real norm = 0;
                    if (normAtoms[kernelPolynomialLookUp.centralAtom[indx]] > constants::EPS_TOL)
                        norm =
                            (Real)kernelPower / normAtoms[kernelPolynomialLookUp.centralAtom[indx]];

                    Real factor = math::intPowMixAlgo(
                        polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                        [kernelPolynomialLookUp.kernelEntry[indx]],
                        kernelPower - 1);
                    polynomialDerivative[kernelPolynomialLookUp.typeStruc[indx]]
                                        [kernelPolynomialLookUp.kernelEntry[indx]] = factor * norm;
                    polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                    [kernelPolynomialLookUp.kernelEntry[indx]] =
                                        factor
                                        * polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                                          [kernelPolynomialLookUp
                                                               .kernelEntry[indx]];
                    polynomialKernelNorm[kernelPolynomialLookUp.typeStruc
                                             [indx]][kernelPolynomialLookUp.kernelEntry[indx]] =
                        norm
                        * polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                          [kernelPolynomialLookUp.kernelEntry[indx]];
                });
        },
        __func__);
}

void KernelPolynomial::computePolynomialKernelTermsRefit_stdlib()
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();

    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                kernelPolynomialLookUp.mainIndex.cbegin(),
                kernelPolynomialLookUp.mainIndex.cend(),
                [&](const Size& indx)
                {
                    Real factor1 = (Real)kernelPower;
                    if (normAtoms[kernelPolynomialLookUp.centralAtom[indx]] > constants::EPS_TOL)
                        factor1 = factor1 / normAtoms[kernelPolynomialLookUp.centralAtom[indx]];

                    Real factor = math::intPowMixAlgo(
                        polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                        [kernelPolynomialLookUp.kernelEntry[indx]],
                        kernelPower - 1);
                    polynomialDerivative[kernelPolynomialLookUp.typeStruc[indx]]
                                        [kernelPolynomialLookUp.kernelEntry[indx]] =
                                            factor * factor1;
                    polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                    [kernelPolynomialLookUp.kernelEntry[indx]] =
                                        factor
                                        * polynomialKernel[kernelPolynomialLookUp.typeStruc[indx]]
                                                          [kernelPolynomialLookUp
                                                               .kernelEntry[indx]];
                });
        },
        __func__);
}
