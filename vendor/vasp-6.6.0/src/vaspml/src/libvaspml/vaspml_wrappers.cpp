#include "vaspml_wrappers.hpp"

#include "InterfaceVasp.hpp"

using namespace vaspml;

extern "C"
{

InterfaceVasp* vaspml_createHandle()
{
    return new InterfaceVasp();
}

void vaspml_destroyHandle(InterfaceVasp** handle)
{
    delete *handle;
    *handle = nullptr;
}

void vaspml_prepareSetup(InterfaceVasp** handle,
                         const int*      ntyp,
                         const int*      nsw,
                         const int*      ibrion,
                         const int8_t*   lmlabexist,
                         const int8_t*   lmlffexist)
{
    (*handle)->prepareSetup(ntyp, nsw, ibrion, lmlabexist, lmlffexist);
}

void vaspml_readIncar(InterfaceVasp** handle, const char* incarAsChars)
{
    (*handle)->readIncar(incarAsChars);
}

void vaspml_setupFromIncar(InterfaceVasp** handle, const char* typesAsChars)
{
    (*handle)->setupFromIncar(typesAsChars);
}

void vaspml_updateSetup(InterfaceVasp** handle, const int* outblockSync)
{
    (*handle)->updateSetup(outblockSync);
}

void vaspml_setupMpi(InterfaceVasp** handle, MPI_Fint* mpiComm)
{
    (*handle)->setupMpi(mpiComm);
}

void vaspml_setupForceField(InterfaceVasp** handle)
{
    (*handle)->setupForceField();
}

double vaspml_get_W1(InterfaceVasp** handle)
{
    return (*handle)->get_W1();
}

double vaspml_get_W2(InterfaceVasp** handle)
{
    return (*handle)->get_W2();
}

double vaspml_get_RCUT1(InterfaceVasp** handle)
{
    return (*handle)->get_RCUT1();
}

double vaspml_get_RCUT2(InterfaceVasp** handle)
{
    return (*handle)->get_RCUT2();
}

void vaspml_set_typeMap(InterfaceVasp** handle, const char* typesAsChars)
{
    (*handle)->set_typeMap(typesAsChars);
}

void vaspml_set_nAtomsType(InterfaceVasp** handle,
                           const int*      ntypes,
                           const int*      nAtomsType,
                           const char*     keyIn)
{
    (*handle)->set_nAtomsType(ntypes, nAtomsType, keyIn);
}

void vaspml_resizeNeighborArrays(InterfaceVasp** handle, const int* nions, const char* keyIn)
{
    (*handle)->resizeNeighborArrays(nions, keyIn);
}

void vaspml_fillNeighborArraysSingleAtom(InterfaceVasp** handle,
                                         const int*      numberNeighbors,
                                         const int*      atomNumber,
                                         const int*      centralType,
                                         const int*      neighborIndex,
                                         const int*      neighborTypes,
                                         const double*   neighborDist,
                                         const double*   neighborConnect,
                                         const char*     keyIn)
{
    (*handle)->fillNeighborArraysSingleAtom(numberNeighbors,
                                            atomNumber,
                                            centralType,
                                            neighborIndex,
                                            neighborTypes,
                                            neighborDist,
                                            neighborConnect,
                                            keyIn);
}

void vaspml_fillForceSingleAtom(InterfaceVasp** handle,
                                const int*      nions,
                                const int*      ionIndex,
                                const int*      localIndex,
                                double*         forces)
{
    (*handle)->fillForceSingleAtom(nions, ionIndex, localIndex, forces);
}

void vaspml_computeForces(InterfaceVasp** handle, const int* nions, double* forces)
{
    (*handle)->computeForces(nions, forces);
}

void vaspml_makeGlobalIonIndex(InterfaceVasp** handle, const int* nions)
{
    (*handle)->makeGlobalIonIndex(nions);
}

void vaspml_update(InterfaceVasp** handle, const double* volume)
{
    (*handle)->update(volume);
}

void vaspml_getPotentialEnergy(InterfaceVasp** handle, double* energy)
{
    (*handle)->getPotentialEnergy(energy);
}

void vaspml_getStressTensor(InterfaceVasp** handle, double* stress)
{
    (*handle)->getStressTensor(stress);
}

void vaspml_resetNeighborArrays(InterfaceVasp** handle, const char* keyIn)
{
    (*handle)->resetNeighborArrays(keyIn);
}

void vaspml_copyNeighborArrays(InterfaceVasp** handle)
{
    (*handle)->copyNeighborArrays();
}

void vaspml_copyForceArrayToCPU(InterfaceVasp** handle)
{
    (*handle)->copyForceArrayToCPU();
}

void vaspml_copyForceArrayToGPU(InterfaceVasp** handle)
{
    (*handle)->copyForceArrayToGPU();
}

} // extern "C"
