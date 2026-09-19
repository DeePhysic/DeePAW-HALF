#ifndef VASPML_WRAPPERS_HPP
#define VASPML_WRAPPERS_HPP

#include "mpi.hpp"

#include <cstdint>

namespace vaspml
{
class InterfaceVasp;
}

extern "C"
{

vaspml::InterfaceVasp* vaspml_createHandle();
void                   vaspml_destroyHandle(vaspml::InterfaceVasp** handle);
/*================================================================================================+
 | Wrappers for InterfaceVasp member functions.
 +================================================================================================*/
void   vaspml_prepareSetup(vaspml::InterfaceVasp** handle,
                           const int*              ntyp,
                           const int*              nsw,
                           const int*              ibrion,
                           const int8_t*           lmlabexist,
                           const int8_t*           lmlffexist);
void   vaspml_readIncar(vaspml::InterfaceVasp** handle, const char* incarAsChars);
void   vaspml_setupFromIncar(vaspml::InterfaceVasp** handle, const char* typesAsChars);
void   vaspml_updateSetup(vaspml::InterfaceVasp** handle, const int* outblockSync);
void   vaspml_setupMpi(vaspml::InterfaceVasp** handle, MPI_Fint* mpiComm);
void   vaspml_setupForceField(vaspml::InterfaceVasp** handle);
double vaspml_get_W1(vaspml::InterfaceVasp** handle);
double vaspml_get_W2(vaspml::InterfaceVasp** handle);
double vaspml_get_RCUT1(vaspml::InterfaceVasp** handle);
double vaspml_get_RCUT2(vaspml::InterfaceVasp** handle);
void   vaspml_set_typeMap(vaspml::InterfaceVasp** handle, const char* typesAsChars);
void   vaspml_set_nAtomsType(vaspml::InterfaceVasp** handle,
                             const int*              ntypes,
                             const int*              nAtomsType,
                             const char*             keyIn);
void   vaspml_resizeNeighborArrays(vaspml::InterfaceVasp** handle,
                                   const int*              nions,
                                   const char*             keyIn);
void   vaspml_fillNeighborArraysSingleAtom(vaspml::InterfaceVasp** handle,
                                           const int*              numberNeighbors,
                                           const int*              atomNumber,
                                           const int*              centralType,
                                           const int*              neighborIndex,
                                           const int*              neighborTypes,
                                           const double*           neighborDist,
                                           const double*           neighborConnect,
                                           const char*             keyIn);
void   vaspml_fillForceSingleAtom(vaspml::InterfaceVasp** handle,
                                  const int*              nions,
                                  const int*              ionIndex,
                                  const int*              localIndex,
                                  double*                 forces);
void   vaspml_computeForces(vaspml::InterfaceVasp** handle, const int* nions, double* forces);
void   vaspml_makeGlobalIonIndex(vaspml::InterfaceVasp** handle, const int* nions);
void   vaspml_update(vaspml::InterfaceVasp** handle, const double* volume);
void   vaspml_getPotentialEnergy(vaspml::InterfaceVasp** handle, double* energy);
void   vaspml_getStressTensor(vaspml::InterfaceVasp** handle, double* stress);
void   vaspml_resetNeighborArrays(vaspml::InterfaceVasp** handle, const char* keyIn);
void   vaspml_copyNeighborArrays(vaspml::InterfaceVasp** handle);
void   vaspml_copyForceArrayToCPU(vaspml::InterfaceVasp** handle);
void   vaspml_copyForceArrayToGPU(vaspml::InterfaceVasp** handle);

} // extern "C"

#endif
