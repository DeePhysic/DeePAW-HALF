#ifndef INTERFACELAMMPS_HPP
#define INTERFACELAMMPS_HPP

#include "mpi.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace vaspml
{

class InterfaceLAMMPS
{
  public:
    /// Constructor, creates Impl (Frame) class instance.
    InterfaceLAMMPS();
    /// Default destructor.
    ~InterfaceLAMMPS();
    /// Copy constructor, deleted (we do not allow copying).
    InterfaceLAMMPS(InterfaceLAMMPS const& other) = delete;
    /// Copy assignment operator, deleted (we do not allow copying).
    InterfaceLAMMPS& operator=(InterfaceLAMMPS const& other) = delete;
    /// Move constructor.
    InterfaceLAMMPS(InterfaceLAMMPS&& other);
    /// Move assignment operator.
    InterfaceLAMMPS& operator=(InterfaceLAMMPS&& other);
    /***********************************************************************************************
     * Set up the force field from a given file.
     *
     * @param forceFieldFileName Name of force field file.
     *
     * @return Maximum cutoff radius required for this force field.
     **********************************************************************************************/
    double setupForceField(std::string forceFieldFileName, const MPI_Comm world);
    /***********************************************************************************************
     * Forward LAMMPS types which are mapped to VASP types.
     *
     * @param typesToNames Vector with LAMMPS type numbers map to VASP type strings.
     * @return Multi-line string containing a type mapping report (see below).
     *
     * Consider the following example:
     * * LAMMPS `pair_coeff` invoked with this type string: `NULL Cs NULL Br Pb Br`.
     * * VASP ML_FF file with the following types present in this order: `Ca, Pb, O, Br, Cs`.
     * Then we would like to achieve this mapping (this is also the output of the report):
     *
     *
     *    LAMMPS       pair_coeff      VASP      |             VASP force field
     *     types       names           subtypes  |     types       name         subtypes
     * ----------------------------------------- | -------------------------------------
     *         1 <---> unmapped! <---> unmapped! |         1 <---> Ca     <---> unused!
     *         2 <---> Cs        <---> 3         |         2 <---> Pb     <---> 1
     *         3 <---> unmapped! <---> unmapped! |         3 <---> O      <---> unused!
     *         4 <---> Br        <---> 2         |         4 <---> Br     <---> 2
     *         5 <---> Pb        <---> 1         |         5 <---> Cs     <---> 3
     *         6 <---> Br        <---> 2         |
     *
     **********************************************************************************************/
    std::string setupTypes(std::vector<std::string> typesToNames);
    /***********************************************************************************************
     * Transfer LAMMPS atom tags.
     **********************************************************************************************/
    void setTags(const int* const& tag);
    void setTags(const int64_t* const& tag);
    /***********************************************************************************************
     * Forward LAMMPS neighbor list and compute internal neighbor lists.
     *
     * @param x 2-dimensional array containing the positions of all local + ghost atoms.
     * @param type 1-dimensional array containing the LAMMPS types of all local + ghost atoms.
     * @param inum Number of atoms for which neighbor list entries are present (this may be smaller
     *             than `nlocal` in case there are atoms which are not mapped to types in the force
     *             field, i.e., have `NULL` entry in `pair_coeff` line).
     * @param ilist 1-dimensional array containing the indices of central atoms (in the x and type
     *              arrays) for which neighbor lists are present.
     * @param numneigh 1-dimensional array containing the number of neighbors for each central atom.
     * @param firstneigh 2-dimensional array containing the indices of neighbor atoms (in the x and
     *                   type arrays) for each central atom (index 1: central atom index (result of
     *                   ilist), index 2: neighbor index).
     * @param neighmask Integer mask to remove bits which mark special neighbors in LAMMPS.
     **********************************************************************************************/
    void computeNeighborLists(const double* const* const& x,
                              const int* const&           type,
                              const int&                  inum,
                              const int* const&           ilist,
                              const int* const&           numneigh,
                              const int* const* const&    firstneigh,
                              const int                   neighmask);
    /***********************************************************************************************
     * Perform force field prediction.
     **********************************************************************************************/
    void predict();
    /***********************************************************************************************
     * Get potential energy prediction summed over all local atoms.
     **********************************************************************************************/
    double getLocalEnergy() const;
    /***********************************************************************************************
     * Retrieve forces and store them in LAMMPS arrays.
     **********************************************************************************************/
    void getForces(double* const* const& f) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
};

} // namespace vaspml

#endif
