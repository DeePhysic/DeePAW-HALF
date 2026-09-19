#ifndef FORCEFIELDGRACE_HPP
#define FORCEFIELDGRACE_HPP

#include "ForceField.hpp"
#include "GraceInterface.hpp"
#include "Record.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm> // for fill

namespace vaspml
{

class ForceFieldGrace : public ForceField
{
  public:
    ForceFieldGrace(String modelName) :
        interface(modelName, 1.2),
        neighborListData(std::make_shared<Record>()),
        neighborList(std::make_shared<NearestNeighborNSquare>(interface.getMaximumCutoff()
                                                                  / constants::AUTOA,
                                                              true,
                                                              false,
                                                              neighborListData)) {};

    void updateSetup(Record& setup) final
    {
        setup.get<Real>("ML_RCUT1") = get_RCUT1() * constants::AUTOA;

        return;
    }

    Real get_RCUT1() const final { return interface.getMaximumCutoff() / constants::AUTOA; }

    void set_typeMap(Vec1String types) final
    {
        typeMap = types;
        // Element strings from Fortran may contain spaces!
        for (String& t : typeMap) t = string_tools::trim(t);

        return;
    }

    void set_nAtomsType(Int ntypes, const int* nAtomsType, String /* key */) final
    {
        neighborList->set_nTypes(ntypes);
        Vec1Int& nAtomsPerType = neighborListData->get<Vec1Int>("nAtomsType");
        nAtomsPerType.resize(ntypes);
        for (Size type = 0; type < (Size)ntypes; type++) { nAtomsPerType[type] = nAtomsType[type]; }

        return;
    }

    void resizeNeighborArrays(Int nions, String /* key */) final
    {
        neighborListData->get<Vec2Int>("globalIndex").resize(nions);
        neighborListData->get<Vec2Int>("typeIndex").resize(nions);
        neighborListData->get<Vec1Int>("typeIndexCentral").resize(nions);
        neighborListData->get<Vec2Real>("distances").resize(nions);
        neighborListData->get<Vec2Real>("connectionVector").resize(nions);
        neighborListData->get<Vec2Real>("connectionVectorNormalized").resize(nions);
        neighborListData->get<Vec1Int>("numberNeighbors").resize(nions);
        neighborListData->get<Vec2Int>("numberNeighborsType").resize(nions);
        neighborListData->get<Vec1Int>("centralAtomIndexPerType").resize(nions);
        neighborList->set_nAtoms(nions);

        return;
    }

    void fillNeighborArraysSingleAtom(const Int     numberNeighbors,
                                      const Int     atomNumber,
                                      const Int     centralType,
                                      const int*    neighborIndex,
                                      const int*    neighborTypes,
                                      const double* neighborDist,
                                      const double* neighborConnect,
                                      String /* key */) final
    {
        // TODO: copy from InterfaceVASP::fillNeighborArraysSingleAtomCPU
        // in the future probably better to avoid the whole
        // NearestNeighborNSquare business and fill our list directly from here
        neighborListData->get<Vec1Int>("typeIndexCentral")[atomNumber] = centralType;
        neighborListData->get<Vec1Int>("numberNeighbors")[atomNumber] = numberNeighbors;

        // unpack arrays
        Vec1Int&  index = neighborListData->get<Vec2Int>("globalIndex")[atomNumber];
        Vec1Int&  nTypes = neighborListData->get<Vec2Int>("typeIndex")[atomNumber];
        Vec1Real& nDistance = neighborListData->get<Vec2Real>("distances")[atomNumber];
        Vec1Real& nVectorNorm =
            neighborListData->get<Vec2Real>("connectionVectorNormalized")[atomNumber];
        Vec1Real& nVector = neighborListData->get<Vec2Real>("connectionVector")[atomNumber];
        Vec1Int&  nTypeCounter = neighborListData->get<Vec2Int>("numberNeighborsType")[atomNumber];

        //resize arrays
        index.resize(numberNeighbors);
        nTypes.resize(numberNeighbors);
        nDistance.resize(numberNeighbors);
        nVectorNorm.resize(3 * numberNeighbors);
        nVector.resize(3 * numberNeighbors);
        nTypeCounter.resize(neighborListData->get<Vec1Int>("nAtomsType").size());
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

    void fillForceSingleAtom(Int /* nions */,
                             Int     ionIndex,
                             Int     localIndex,
                             double* forces) const final
    {
        double* AtomForce = forces + 3 * ionIndex;
        for (auto i = 0; i < 3; ++i)
        {
            AtomForce[i] = interface.getForces()[3 * localIndex + i];
            AtomForce[i] /= vaspml::constants::FUNIT;
        }

        return;
    }

    void update(Real volume) final
    {
        interface.updateNeighborList(typeMap, *neighborList);
        interface.update();
        this->volume = volume;

        return;
    }

    void getPotentialEnergy(double* energy) const final
    {
        *energy += interface.getEnergy() / vaspml::constants::EUNIT;

        return;
    }

    void getStressTensor(double* stress) const final
    {
        const Size XX = 0;
        const Size YY = 1;
        const Size ZZ = 2;
        const Size XY = 3;
        const Size YX = 3;
        const Size XZ = 4;
        const Size ZX = 4;
        const Size YZ = 5;
        const Size ZY = 5;

        double factor = 0.5 / (vaspml::constants::RYTOEV * volume);
        stress[0] = factor * interface.getStressTensor()[XX];
        stress[1] = factor * interface.getStressTensor()[XY]; // (0, 1);
        stress[2] = factor * interface.getStressTensor()[XZ]; // (0, 2);

        stress[3] = factor * interface.getStressTensor()[YX]; // (1, 0);
        stress[4] = factor * interface.getStressTensor()[YY]; // (1, 1);
        stress[5] = factor * interface.getStressTensor()[YZ]; // (1, 2);

        stress[6] = factor * interface.getStressTensor()[ZX]; // (2, 0);
        stress[7] = factor * interface.getStressTensor()[ZY]; // (2, 1);
        stress[8] = factor * interface.getStressTensor()[ZZ]; // (2, 2);

        return;
    }

  private:
    GraceInterface                          interface;
    ShRec                                   neighborListData;
    std::shared_ptr<NearestNeighborNSquare> neighborList;
    Vec1String                              typeMap;
    double                                  volume;
};

} // namespace vaspml

#endif
