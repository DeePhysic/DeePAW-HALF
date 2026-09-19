#ifndef FORCEFIELD_HPP
#define FORCEFIELD_HPP

#include "types.hpp"

namespace vaspml
{

struct Record;

class ForceField
{
  public:
    virtual ~ForceField() = default;
    // TODO: This should be const but isn't currently because of unit transformation.
    virtual void updateSetup(Record& /* setup */) {};
    virtual Real get_W1() const { return 1.0; };
    virtual Real get_W2() const { return 0.0; };
    virtual Real get_RCUT1() const = 0;
    virtual Real get_RCUT2() const { return 0.0; };
    virtual void set_typeMap(Vec1String types) = 0;
    virtual void set_nAtomsType(Int ntypes, const int* nAtomsType, String key) = 0;
    virtual void resizeNeighborArrays(Int nions, String key) = 0;
    virtual void fillNeighborArraysSingleAtom(Int           numberNeighbors,
                                              Int           atomNumber,
                                              Int           centralType,
                                              const int*    neighborIndex,
                                              const int*    neighborTypes,
                                              const double* neighborDist,
                                              const double* neighborConnect,
                                              String        key) = 0;
    virtual void fillForceSingleAtom(Int     nions,
                                     Int     ionIndex,
                                     Int     localIndex,
                                     double* forces) const = 0;
    // TODO: This does not seem to be used currently.
    virtual void computeForces(Int /* nions */, double* /* forces */) {};
    virtual void makeGlobalIonIndex(Int /* nions */){};
    virtual void update(Real volume) = 0;
    virtual void getPotentialEnergy(double* energy) const = 0;
    virtual void getStressTensor(double* stress) const = 0;
};

} // namespace vaspml

#endif
