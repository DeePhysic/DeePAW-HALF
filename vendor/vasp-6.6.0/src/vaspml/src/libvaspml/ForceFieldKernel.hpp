#ifndef FORCEFIELDKERNEL_HPP
#define FORCEFIELDKERNEL_HPP

#include "MlMPI.hpp"

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "ForceField.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "Predictor.hpp"
#include "Record.hpp"
#include "Timer.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "nearest_neighbor.hpp"
#include "settings.hpp"
#include "setup.hpp"
#include "types.hpp"

#include <algorithm> // for fill

namespace vaspml
{

class BasisFunctions;

class ForceFieldKernel : public ForceField
{
  public:
    ForceFieldKernel(std::shared_ptr<MlMPI> mpiIn) : mpiMain(mpiIn)
    {
        Record& ffd = forceFieldData;
        io::readMlffAndConvertUnits(ffd);

        basisFunctions = setup::makeBasisFunctions(ffd);

        frame.init(ffd, basisFunctions);
        typeMapPtr = std::make_shared<TypeMap>();
        polynomialKernel.init(ffd, mpiMain);
        predictor.init(ffd, mpiMain);
    }

    void updateSetup(Record& setup) final
    {
        settings::updateFromMlff(forceFieldData, setup);
        settings::checkFeatureSupport(setup);

        return;
    }

    Real get_W1() const final { return forceFieldData.cget<Real>("ML_W1"); }
    Real get_W2() const final { return forceFieldData.cget<Real>("ML_W2"); }
    Real get_RCUT1() const final { return forceFieldData.cget<Real>("ML_RCUT1"); }
    Real get_RCUT2() const final { return forceFieldData.cget<Real>("ML_RCUT2"); }

    void set_typeMap(Vec1String types) final
    {
        Vec1String ffTypes = forceFieldData.cget<Vec1String>("types");
        typeMapPtr->update(ffTypes, types);

        return;
    }

    void set_nAtomsType(Int ntypes, const int* nAtomsType, String key) final
    {
        frame.neighborLists[key]->set_nTypes(ntypes);
        frame.neighborLists[key]->set_totalNumberTypes(ntypes);
        Vec1Int& nAtomsPerType =
            frame.data->get<ShRec>(key + "-neighbor-list")->get<Vec1Int>("nAtomsType");
        nAtomsPerType.resize(ntypes);
        for (Size type = 0; type < (Size)ntypes; type++) { nAtomsPerType[type] = nAtomsType[type]; }

        return;
    }

    void resizeNeighborArrays(Int nions, String key) final
    {
        ShRec& dataNList = frame.data->get<ShRec>(key + "-neighbor-list");
        dataNList->get<Vec2Int>("globalIndex").resize(nions);
        dataNList->get<Vec2Int>("typeIndex").resize(nions);
        dataNList->get<Vec1Int>("typeIndexCentral").resize(nions);
        dataNList->get<Vec2Real>("distances").resize(nions);
        dataNList->get<Vec2Real>("connectionVector").resize(nions);
        dataNList->get<Vec2Real>("connectionVectorNormalized").resize(nions);
        dataNList->get<Vec1Int>("numberNeighbors").resize(nions);
        dataNList->get<Vec2Int>("numberNeighborsType").resize(nions);
        dataNList->get<Vec1Int>("centralAtomIndexPerType").resize(nions);
        frame.neighborLists[key]->set_nAtoms(nions);

        return;
    }

    void fillNeighborArraysSingleAtom(Int           numberNeighbors,
                                      Int           atomNumber,
                                      Int           centralType,
                                      const int*    neighborIndex,
                                      const int*    neighborTypes,
                                      const double* neighborDist,
                                      const double* neighborConnect,
                                      String        key) final
    {
        ShRec& dataNList = frame.data->get<ShRec>(key + "-neighbor-list");
        dataNList->get<Vec1Int>("typeIndexCentral")[atomNumber] = centralType;
        dataNList->get<Vec1Int>("numberNeighbors")[atomNumber] = numberNeighbors;

        // unpack arrays
        Vec1Int&  index = dataNList->get<Vec2Int>("globalIndex")[atomNumber];
        Vec1Int&  nTypes = dataNList->get<Vec2Int>("typeIndex")[atomNumber];
        Vec1Real& nDistance = dataNList->get<Vec2Real>("distances")[atomNumber];
        Vec1Real& nVectorNorm = dataNList->get<Vec2Real>("connectionVectorNormalized")[atomNumber];
        Vec1Real& nVector = dataNList->get<Vec2Real>("connectionVector")[atomNumber];
        Vec1Int&  nTypeCounter = dataNList->get<Vec2Int>("numberNeighborsType")[atomNumber];

        //resize arrays
        index.resize(numberNeighbors);
        nTypes.resize(numberNeighbors);
        nDistance.resize(numberNeighbors);
        nVectorNorm.resize(3 * numberNeighbors);
        nVector.resize(3 * numberNeighbors);
        nTypeCounter.resize(dataNList->get<Vec1Int>("nAtomsType").size());
        std::fill(nTypeCounter.begin(), nTypeCounter.end(), (Real)0);

        for (Size nIndex = 0; nIndex < (Size)numberNeighbors; nIndex++)
        {
            index[nIndex] = neighborIndex[nIndex];
            nTypes[nIndex] = neighborTypes[nIndex];
            nDistance[nIndex] = neighborDist[nIndex];
            nVectorNorm[3 * nIndex] = neighborConnect[3 * nIndex];
            nVectorNorm[3 * nIndex + 1] = neighborConnect[3 * nIndex + 1];
            nVectorNorm[3 * nIndex + 2] = neighborConnect[3 * nIndex + 2];
            nVector[3 * nIndex] = neighborConnect[3 * nIndex] * neighborDist[nIndex];
            nVector[3 * nIndex + 1] = neighborConnect[3 * nIndex + 1] * neighborDist[nIndex];
            nVector[3 * nIndex + 2] = neighborConnect[3 * nIndex + 2] * neighborDist[nIndex];
            // count number of neighors per type ( not really needed )
            nTypeCounter[neighborTypes[nIndex]]++;
        }

        return;
    }

    void fillForceSingleAtom(Int /* nions */, Int ionIndex, Int localIndex, double* forces) const
    {
        //nvtxRangePush( "ForceCompute" );
        const auto& descriptorCollection = polynomialKernel.get_descriptorCollection();
        const std::map<String, ShVec1Real>& centralForces = predictor.get_centralForces();
        const std::map<String, ShVec2Real>& pairForces = predictor.get_pairForces();
        for (const String& key : constants::descriptorKeyList)
        {
            std::shared_ptr<const NearestNeighborNSquare> neighborList =
                descriptorCollection.getDescriptor(key).get_neighborList_ptr();
            Vec1Real& centralForcesVec = *(centralForces.at(key));
            Vec2Real& pairForcesVec = *(pairForces.at(key));
            forces[3 * ionIndex] += centralForcesVec[3 * localIndex];
            forces[3 * ionIndex + 1] += centralForcesVec[3 * localIndex + 1];
            forces[3 * ionIndex + 2] += centralForcesVec[3 * localIndex + 2];
            for (Size neighborAtom = 0; neighborAtom < neighborList->get_size(localIndex);
                 neighborAtom++)
            {
                const Int& neighborIndex = neighborList->get_globalIndex(localIndex, neighborAtom);
                forces[3 * neighborIndex] += pairForcesVec[localIndex][3 * neighborAtom];
                forces[3 * neighborIndex + 1] += pairForcesVec[localIndex][3 * neighborAtom + 1];
                forces[3 * neighborIndex + 2] += pairForcesVec[localIndex][3 * neighborAtom + 2];
            }
        }
        //nvtxRangePop();

        return;
    }

    void computeForces(Int nions, double* forces) final
    {
        //VASPML_PROFILING_START("ForceComputeTotal");
        const auto& descriptorCollection = polynomialKernel.get_descriptorCollection();
        //VASPML_PROFILING_START("ForceComputeComputation");
        predictor.compute_atomicForces(descriptorCollection, globalIonIndx, nions);
        //VASPML_PROFILING_STOP("ForceComputeComputation");
        // copy array into vasp array
        const Vec1Real& forcePredicted = predictor.get_atomicForces();
        //VASPML_PROFILING_START("ForceComputeCopy");
        for (Int i = 0; i < 3 * nions; i++) { forces[i] = forcePredicted[i]; }
        //VASPML_PROFILING_STOP("ForceComputeCopy");
        //VASPML_PROFILING_STOP("ForceComputeTotal");

        return;
    }

    void makeGlobalIonIndex(Int nions) final
    {
        Int nCPU;
        Int node;
        if (mpiMain != nullptr)
        {
            nCPU = mpiMain->get_numberRanks();
            node = mpiMain->get_rank();
        }
        else
        {
            nCPU = 1;
            node = 0;
        }
        globalIonIndx.resize(nions, -1);
        for (Int indxTmp = 0; indxTmp < nions; indxTmp += nCPU)
        {
            Int ion = indxTmp + node;
            if (ion < nions) { globalIonIndx[indxTmp] = ion; }
        }

        return;
    }

    void update(Real volume) final
    {
        VASPML_PROFILING_START("UPDATE");
        VASPML_PROFILING_START("Descriptors");
        frame.update(typeMapPtr);
        VASPML_PROFILING_STOP("Descriptors");
        VASPML_PROFILING_START("Kernel");
        polynomialKernel.updatePolynomialKernel(frame);
        VASPML_PROFILING_STOP("Kernel");
        VASPML_PROFILING_START("Predictor");
        predictor.update(polynomialKernel);
        VASPML_PROFILING_STOP("Predictor");
        VASPML_PROFILING_START("Stress");
        predictor.compute_stressTensor(polynomialKernel.get_descriptorCollection(), volume);
        VASPML_PROFILING_STOP("Stress");
        VASPML_PROFILING_STOP("UPDATE");

        return;
    }

    void getPotentialEnergy(double* energy) const final
    {
        (*energy) += predictor.get_totalEnergy();

        return;
    }

    void getStressTensor(double* stress) const final
    {
        stress[0] = predictor.get_totalStressTensor(0, 0);
        stress[1] = predictor.get_totalStressTensor(0, 1);
        stress[2] = predictor.get_totalStressTensor(0, 2);

        stress[3] = predictor.get_totalStressTensor(1, 0);
        stress[4] = predictor.get_totalStressTensor(1, 1);
        stress[5] = predictor.get_totalStressTensor(1, 2);

        stress[6] = predictor.get_totalStressTensor(2, 0);
        stress[7] = predictor.get_totalStressTensor(2, 1);
        stress[8] = predictor.get_totalStressTensor(2, 2);

        return;
    }

  private:
    Frame frame;
    /*******************************************************************************************
     * manages the input of the ML_FF file. Stores and distributes the input data
     * needed for the force field computation
     *******************************************************************************************/
    Record forceFieldData;
    /*******************************************************************************************
     * main mpi class for vaspml
     *
     * all other mpi communicators have to be created from this instance.
     * From this instance we can create shared memory mpi communicators
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mpiMain = nullptr;
    /*******************************************************************************************
     * stores the basis functions as spherical harmonics and modified spherical Bessel functions.
     *
     * This data structure handles the basis set functions for DescriptorSHS2 and DescriptorSHS3
     *******************************************************************************************/
    std::map<String, std::shared_ptr<BasisFunctions>> basisFunctions;
    /*******************************************************************************************
     * class takes the descriptors DescriptorSHS2 and DescriptorSHS3 as input
     * and computes the kernel matrix and derivatives of the specific kernel matrix.
     *
     * The kernel matrix is computed as a similarity meassure between local reference
     * configurations and local configurations of the current structure.
     * A local reference configuration is a local representation ( within some cutoff ) of the
     * 2 and 3-body interaction in terms of the classes DescriptorSHS2 and DescriptorSHS3
     *******************************************************************************************/
    KernelPolynomial polynomialKernel;
    /*******************************************************************************************
     * class takes the polyonimalKernel as an input and computes energy, forces and stress tensor
     *
     * The class uses the polynomialKernel and the fitting weights obtained from the ML_FF
     * to compute the total energy of the system
     * The forces are computed from the derivative of the polyonomialKernel.
     * The stress tensor is obtained from the computed pairForces with the help
     * of the virial theorem
     *******************************************************************************************/
    Predictor predictor;
    /*******************************************************************************************
     * this class manages the correspondence between atom types in the machine learning force field
     * and the structure which is currently analysed
     *
     * The types in the structure always have to be a subset of the types present in the machine
     * learning force field. Otherwise, it is not possible to make any predictions.
     * Because this would indicate there are types in the structure not considered during training
     *******************************************************************************************/
    std::shared_ptr<TypeMap> typeMapPtr;
    /*******************************************************************************************
     * Mapping from locally distributed ion index to the global all atom index.
     *******************************************************************************************/
    Vec1Int globalIonIndx;
};

} // namespace vaspml

#endif
