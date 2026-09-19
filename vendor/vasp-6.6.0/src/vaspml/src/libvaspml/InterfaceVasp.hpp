#ifndef INTERFACEVASP_HPP
#define INTERFACEVASP_HPP

#include "mpi.hpp"

#include "ForceField.hpp"
#include "types.hpp"

namespace vaspml
{

class MlMPI;
struct Record;

class InterfaceVasp
{
  public:
    InterfaceVasp();
    void   prepareSetup(const int*    ntyp,
                        const int*    nsw,
                        const int*    ibrion,
                        const int8_t* lmlabexist,
                        const int8_t* lmlffexist);
    void   readIncar(const char* incarAsChars);
    void   setupFromIncar(const char* typesAsChars);
    void   updateSetup(const int* outblockSync);
    void   setupMpi(MPI_Fint* mpiComm);
    void   setupForceField();
    double get_W1() const;
    double get_W2() const;
    double get_RCUT1() const;
    double get_RCUT2() const;
    void   set_typeMap(const char* typesAsChars);
    void   set_nAtomsType(const int* ntypes, const int* nAtomsType, const char* keyIn) const;
    void   resizeNeighborArrays(const int* nions, const char* keyIn);
    void   fillNeighborArraysSingleAtom(const int*    numberNeighbors,
                                        const int*    atomNumber,
                                        const int*    centralType,
                                        const int*    neighborIndex,
                                        const int*    neighborTypes,
                                        const double* neighborDist,
                                        const double* neighborConnect,
                                        const char*   keyIn);
    void   fillForceSingleAtom(const int* nions,
                               const int* ionIndex,
                               const int* localIndex,
                               double*    forces) const;
    // TODO: This does not seem to be used currently.
    void computeForces(const int* nions, double* forces) const;
    void makeGlobalIonIndex(const int* nions);
    void update(const double* volume);
    void getPotentialEnergy(double* energy) const;
    void getStressTensor(double* stress) const;

    // TODO: Members in this block do not seem to be used currently, keeping them empty for now.
    void resetNeighborArrays(const char* /* keyIn */) const {};
    void copyNeighborArrays() const {};
    void copyForceArrayToCPU() const {};
    void copyForceArrayToGPU() const {};

  private:
    ShRec                       data;
    Record&                     setup;
    Record&                     incar;
    std::unique_ptr<ForceField> forceField;
    std::shared_ptr<MlMPI>      mpi;
};

/**************************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **************************************************************************************************/
MapString dataMapInterfaceVasp();

/**************************************************************************************************
 * Create a Record with "setup" and "incar" entries.
 **************************************************************************************************/
ShRec createSetupRecord();

} // namespace vaspml

#endif
