#include "Predictor.hpp" // IWYU pragma: associated

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "ParallelEnvironment.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"

#include <algorithm> // for fill, for_each
#ifndef VASPML_NO_STDEXECUTION
#include <numeric> // for transform_reduce
#endif

using namespace vaspml;

void StressTensorLookUpTable::update(const Size& nAtoms)
{
    bool resize = (nAtomsOld != nAtoms);
    if (resize)
    {
        nAtomsOld = nAtoms;
        const Size nElements = nAtoms * constants::descriptorKeyList.size();
        mainIndex.resize(nElements);
        centralAtomIndex.resize(nElements);
        descriptorIndex.resize(nElements);
        refill(nAtoms);
    }

    return;
}

void StressTensorLookUpTable::refill(const Size& nAtoms)
{
    Size counter = 0;
    for (Size descInd = 0; descInd < constants::descriptorKeyList.size(); descInd++)
    {
        for (Size centralAtom = 0; centralAtom < nAtoms; centralAtom++)
        {
            mainIndex[counter] = counter;
            centralAtomIndex[counter] = centralAtom;
            descriptorIndex[counter] = constants::descriptorKeyList[descInd];
            counter++;
        }
    }

    return;
}

void AtomicForcesLookUpTable::update(const Int& nAtoms)
{
    bool resize = (nAtomsOld != (Size)nAtoms);
    if (resize)
    {
        nAtomsOld = nAtoms;
        centralAtomIndex.resize(nAtoms);
        mainIndex.resize(nAtoms);
        this->globalIndex.resize(nAtoms);
        refill(nAtoms);
    }

    return;
}

void AtomicForcesLookUpTable::refill(const Int& nAtoms)
{
    for (Size centralAtom = 0; centralAtom < (Size)nAtoms; centralAtom++)
    {
        mainIndex[centralAtom] = centralAtom;
        centralAtomIndex[centralAtom] = centralAtom;
        this->globalIndex[centralAtom] = globalIndex[centralAtom];
    }

    return;
}

void Predictor::allocate_energyArrayGPU(const Size numberTypesStruc)
{
    bool resize = energyArraySize.checkResize1Dim(numberTypesStruc);
    if (resize)
    {
        energyArraySize.resizeArray1Dim(energyArray, energyArrayDim);
    } // second dimension can not change in execution mode. numberLocalRefConfs conserved
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(VASPML_POLICY
                          energyArray.begin(),
                          energyArray.begin() + energyArraySize.act1Dim,
                          [&](Vec1Real& slice) { std::fill(slice.begin(), slice.end(), (Real)0); });
        },
        __func__);

    return;
}

void Predictor::allocate_derivativeMatrixGPU(const Vec1Int& numberAtomType)
{
    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        resize = derivativeMatrixSize[key].checkResize1Dim(numberAtomType.size());
        if (resize)
        {
            derivativeMatrixSize[key].resizeArray1Dim((*derivativeMatrix[key]),
                                                      derivativeMatrixDim[key]);
        }
        else
        {
            resize = derivativeMatrixSize[key].checkResize2Dim(derivativeMatrixDim[key]);
            if (resize)
            {
                derivativeMatrixSize[key].resizeArray2Dim((*derivativeMatrix[key]),
                                                          derivativeMatrixDim[key]);
            }
        }
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            {
                std::for_each(VASPML_POLICY
                              derivativeMatrix[key]->begin(),
                              derivativeMatrix[key]->begin() + derivativeMatrixSize[key].act1Dim,
                              [&](Vec1Real& slice)
                              { std::fill(slice.begin(), slice.end(), (Real)0); });
            },
            __func__);
    }

    return;
}

void Predictor::allocate_forcePreContractGPU(const Size numberAtoms)
{
    allocate_centralAtomIndex(numberAtoms);
    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        // maps are not allowed in parallel region which is later needed in fill
        Vec2Real&        forcePreContract = *this->forcePreContract[key];
        ArrayResizing2D& forcePreContractSize = this->forcePreContractSize[key];
        resize = forcePreContractSize.checkResize1Dim(numberAtoms);
        if (resize)
        {
            forcePreContractSize.resizeArray1Dim(forcePreContract, forcePreContractDim[key]);
        }
        else
        {
            resize = forcePreContractSize.checkResize2Dim(forcePreContractDim[key]);
            if (resize)
            {
                forcePreContractSize.resizeArray2Dim(forcePreContract, forcePreContractDim[key]);
            }
        }
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            {
                std::for_each(VASPML_POLICY
                              centralAtomIndex.cbegin(),
                              centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                              [&](const Size& nAtom) mutable
                              {
                                  std::fill(forcePreContract[nAtom].begin(),
                                            forcePreContract[nAtom].begin()
                                                + forcePreContractSize.actSize[nAtom],
                                            (Real)0);
                              });
            },
            __func__);
    }

    return;
}

void Predictor::allocate_pairForcesGPU(const Size& numberAtoms)
{
    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec2Real&        pairForces = *this->pairForces[key];
        ArrayResizing2D& pairForcesSize = this->pairForcesSize[key];
        resize = pairForcesSize.checkResize1Dim(numberAtoms);
        if (resize) { pairForcesSize.resizeArray1Dim(pairForces, pairForcesDim[key]); }
        else
        {
            resize = pairForcesSize.checkResize2Dim(pairForcesDim[key]);
            if (resize) { pairForcesSize.resizeArray2Dim(pairForces, pairForcesDim[key]); }
        }
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            {
                std::for_each(VASPML_POLICY
                              centralAtomIndex.cbegin(),
                              centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                              [&](const Size& indx)
                              {
                                  std::fill(pairForces[indx].begin(),
                                            pairForces[indx].begin() + pairForcesSize.actSize[indx],
                                            (Real)0);
                              });
            },
            __func__);
    }

    return;
}

void Predictor::allocate_stressTensorGPU()
{
    bool resize = totalStressTensorSize.checkResize(9);
    if (resize)
    {
        totalStressTensor.resize(9);
        for (const String& key : constants::descriptorKeyList)
        {
            if (weights[key] <= 0) continue;
            stressTensor[key]->resize(9);
        }
    }
    // parallel to move/ keep the data on the device
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        { std::fill(VASPML_POLICY
                    totalStressTensor.begin(), totalStressTensor.end(), (Real)0); },
        __func__);
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec1Real& stressTensor = *this->stressTensor[key];
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            { std::fill(VASPML_POLICY
                        stressTensor.begin(), stressTensor.end(), (Real)0); },
            __func__);
    }

    return;
}

void Predictor::allocate_atomicForcesGPU()
{
    // all descriptors have to have the same number of central atoms
    // choose one which is initialized
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        bool resize = atomicForcesSize.checkResize(centralForcesSize[key].actDim);
        if (resize) { atomicForces.resize(centralForcesSize[key].actDim); }
        break;
    }
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::fill(VASPML_POLICY
                      atomicForces.begin(),
                      atomicForces.begin() + atomicForcesSize.actDim,
                      (Real)0);
        },
        __func__);

    return;
}

void Predictor::allocate_atomicForcesGPU(const Size& nAtoms)
{
    bool resize = atomicForcesSize.checkResize(3 * nAtoms);
    if (resize)
    {
        //for ( const auto& key : constants::descriptorKeyList )
        //{
        //    atomicForcesDesc[ key ]->resize( 3 * nAtoms );
        //}
        atomicForces.resize(3 * nAtoms);
    }
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::fill(VASPML_POLICY
                      atomicForces.begin(),
                      atomicForces.begin() + atomicForcesSize.actDim,
                      (Real)0);
        },
        __func__);
    //for ( const auto& key : constants::descriptorKeyList )
    //{
    //    global::parallel.run(
    //        [&]([[maybe_unused]] const auto& policy)
    //        {
    //            std::fill(VASPML_POLICY
    //                      atomicForcesDesc[ key ]->begin(),
    //                      atomicForcesDesc[ key ]->begin() + atomicForcesSize.actDim,
    //                      (Real)0);
    //        },
    //        __func__);
    //}

    return;
}

void Predictor::allocate_tempForceVectorGPU()
{
    bool resize = tempForceVectorSize.checkResize1Dim(tempForceVectorDim.size());
    if (resize) { tempForceVectorSize.resizeArray1Dim(tempForceVector, tempForceVectorDim); }
    else
    {
        resize = tempForceVectorSize.checkResize2Dim(tempForceVectorDim);
        if (resize) { tempForceVectorSize.resizeArray2Dim(tempForceVector, tempForceVectorDim); }
    }

    return;
}

void Predictor::compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection)
{
    for (const String& key : constants::descriptorKeyList)
    {
        const NearestNeighborNSquare& neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec2Real&      pairForces = (*this->pairForces[key]);
        Vec1Real&      centralForces = (*this->centralForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();

        // TODO: Introduces a race condition! Temporarily switched to sequential policy.
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            {
                std::for_each(VASPML_POLICY_SEQ
                              centralAtomIndex.cbegin(),
                              centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                              [&](const Size& centralAtom) mutable
                              {
                                  atomicForces[3 * centralAtom] += centralForces[3 * centralAtom];
                                  atomicForces[3 * centralAtom + 1] +=
                                      centralForces[3 * centralAtom + 1];
                                  atomicForces[3 * centralAtom + 2] +=
                                      centralForces[3 * centralAtom + 2];

                                  //SumFunctorForce functor( atomicForces, pairForces[ centralAtom ]
                                  //); std::for_each( neighborIndex[ centralAtom ].cbegin(),
                                  //neighborIndex[ centralAtom ].cend(), functor );
                                  Size counter = 0;
                                  std::for_each(neighborIndex[centralAtom].cbegin(),
                                                neighborIndex[centralAtom].cend(),
                                                [&](const Int& neighborIndex) mutable
                                                {
                                                    atomicForces[3 * neighborIndex] +=
                                                        pairForces[centralAtom][3 * counter];
                                                    atomicForces[3 * neighborIndex + 1] +=
                                                        pairForces[centralAtom][3 * counter + 1];
                                                    atomicForces[3 * neighborIndex + 2] +=
                                                        pairForces[centralAtom][3 * counter + 2];
                                                    counter++;
                                                });
                              });
            },
            __func__);
    }

    return;
}

void Predictor::compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection,
                                        const Vec1Int&             globalIonIndex,
                                        const Size&                nAtomsGlob)
{
    atomicForcesLookUpTable.update(nAtomsGlob);
    //nvtxRangePush( "ForceComputation" );
    global::parallel.run(
        [&]([[maybe_unused]] const auto& policy)
        {
            std::for_each(
                VASPML_POLICY
                atomicForcesLookUpTable.mainIndex.cbegin(),
                atomicForcesLookUpTable.mainIndex.cend(),
                [&](const Size& centralAtom)
                {
                    for (const auto& key : descriptorKeyList)
                    {
                        const NearestNeighborNSquare& neighborList =
                            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
                        Vec2Real& pairForces = (*this->pairForces.at(key));
                        Vec1Real& centralForces = (*this->centralForces.at(key));
                        for (Size centralAtomLoc = 0; centralAtomLoc < neighborList.get_nAtoms();
                             centralAtomLoc++)
                        {
                            const Vec2Int& neighborIndex = neighborList.get_globalIndex();
                            if ((Size)globalIonIndex[centralAtomLoc] == centralAtom)
                            {
                                atomicForces[3 * centralAtom] +=
                                    centralForces[3 * centralAtomLoc] * FUNIT;
                                atomicForces[3 * centralAtom + 1] +=
                                    centralForces[3 * centralAtomLoc + 1] * FUNIT;
                                atomicForces[3 * centralAtom + 2] +=
                                    centralForces[3 * centralAtomLoc + 2] * FUNIT;
                            }
                            for (Size nIndx = 0; nIndx < neighborList.get_size(centralAtomLoc);
                                 nIndx++)
                            {
                                if ((Size)neighborIndex[centralAtomLoc][nIndx] == centralAtom)
                                {
                                    atomicForces[3 * centralAtom] +=
                                        pairForces[centralAtomLoc][3 * nIndx] * FUNIT;
                                    atomicForces[3 * centralAtom + 1] +=
                                        pairForces[centralAtomLoc][3 * nIndx + 1] * FUNIT;
                                    atomicForces[3 * centralAtom + 2] +=
                                        pairForces[centralAtomLoc][3 * nIndx + 2] * FUNIT;
                                }
                            }
                        }
                    }
                });
        },
        __func__);
    //nvtxRangePop();

    return;
}

void Predictor::compute_stressTensorGPU(const DescriptorCollector& descriptorCollection,
                                        const Real&                volume)
{
    inverseVolume = (Real)1.0 / volume;
    const Size& nAtoms =
        descriptorCollection.getDescriptor("SHS2-2-body").get_neighborList_ptr()->get_nAtoms();
    stressTensorLookUpTable.update(nAtoms);
    Vec1Size entries = {0, 1, 2, 4, 5, 8};
    Size     counter = 0;

#ifdef VASPML_NO_STDEXECUTION
    for (Size dir0 = 0; dir0 < 3; dir0++)
    {
        for (Size dir1 = dir0; dir1 < 3; dir1++)
        {
            totalStressTensor[entries[counter]] = 0.0;
            std::for_each(stressTensorLookUpTable.mainIndex.cbegin(),
                          stressTensorLookUpTable.mainIndex.cend(),
                          [&](const Size& indx)
                          {
                              const String& key = stressTensorLookUpTable.descriptorIndex[indx];
                              const Size&   centralAtom =
                                  stressTensorLookUpTable.centralAtomIndex[indx];
                              const NearestNeighborNSquare& neighborList =
                                  *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
                              const Vec2Real& pairForces = (*this->pairForces.at(key));
                              const Vec1Real& connectionVector =
                                  neighborList.get_connectionVector(centralAtom);
                              Real sumLoc = 0.0;
                              for (Size i = 0; i < neighborList.get_size(centralAtom); i++)
                              {
                                  sumLoc += inverseVolume * pairForces[centralAtom][3 * i + dir0]
                                          * connectionVector[3 * i + dir1];
                              }
                              totalStressTensor[entries[counter]] += sumLoc;

                              return;
                          });
            counter++;
        }
    }
#else
    for (Size dir0 = 0; dir0 < 3; dir0++)
    {
        for (Size dir1 = dir0; dir1 < 3; dir1++)
        {
            totalStressTensor[entries[counter]] = global::parallel.run(
                [&]([[maybe_unused]] const auto& policy)
                {
                    return std::transform_reduce(
                        VASPML_POLICY
                        stressTensorLookUpTable.mainIndex.cbegin(),
                        stressTensorLookUpTable.mainIndex.cend(),
                        0.0,
                        [](const Real& a, const Real& b) { return a + b; },
                        [&, dir0 = dir0, dir1 = dir1, counterS = counter](const Size& indx) mutable
                        {
                            const String& key = stressTensorLookUpTable.descriptorIndex[indx];
                            const Size&   centralAtom =
                                stressTensorLookUpTable.centralAtomIndex[indx];
                            const NearestNeighborNSquare& neighborList =
                                *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
                            const Vec2Real& pairForces = (*this->pairForces.at(key));
                            const Vec1Real& connectionVector =
                                neighborList.get_connectionVector(centralAtom);
                            Real& sum = (*this->stressTensor.at(key))[counterS];
                            Real  sumLoc = 0.0;

                            for (Size counter = 0; counter < neighborList.get_size(centralAtom);
                                 counter++)
                            {
                                // NOTE this syntax is needed because there seems to be a compiler
                                // bug compute xx term
                                //stressTensor[0] += inverseVolume * pairForces[centralAtom][3 *
                                //counter + dir0] * connectionVector[3 * counter + dir1];
                                sumLoc += inverseVolume
                                        * pairForces[centralAtom][3 * counter + dir0]
                                        * connectionVector[3 * counter + dir1];
                                sum += inverseVolume * pairForces[centralAtom][3 * counter + dir0]
                                     * connectionVector[3 * counter + dir1];
                            }
                            return sumLoc;
                        });
                },
                __func__);
            counter++;
        }
    }
#endif

    totalStressTensor[3] = totalStressTensor[1];
    totalStressTensor[6] = totalStressTensor[2];
    totalStressTensor[7] = totalStressTensor[5];

    return;
}

void Predictor::allocate_centralForcesGPU(const Int numberIons)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        bool             resize = centralForcesSize[key].checkResize(3 * numberIons);
        Vec1Real&        centralForces = *this->centralForces[key];
        ArrayResizing1D& centralForcesSize = this->centralForcesSize[key];
        if (resize) { centralForces.resize(3 * numberIons); }
        global::parallel.run(
            [&]([[maybe_unused]] const auto& policy)
            {
                std::fill(VASPML_POLICY
                          centralForces.begin(),
                          centralForces.begin() + centralForcesSize.actDim,
                          (Real)0.0);
            },
            __func__);
    }

    return;
}
