#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include "Lattice.hpp"
#include "types.hpp"

#include <fstream>
#include <tuple>

namespace vaspml
{

class TypeMap;

class Structure
{
  public:
    /***********************************************************************************************
     * Constructor, optionally pass smart pointers to existing memory.
     **********************************************************************************************/
    Structure(ShRec record = nullptr, const bool& isDirect = false);
    /***********************************************************************************************
     * Read in structure data from POSCAR file.
     *
     * @param fname Name of input file in POSCAR format.
     **********************************************************************************************/
    void readPoscar(const String& fname);
    /***********************************************************************************************
     * Set lattice vectors and compute inverse lattice matrix.
     *
     * @param lattice_in   3x3 matrix in linearized form as xx,xy,xz,yx,yy,yz,zx,zy,zz
     **********************************************************************************************/
    void set_lattice(const Vec1Real& lattice_in);
    /* set positions and shift them to superell
     * @param positions    -> 3*N array containing the atomic positions in either direct or
     * cartesian coordinates
     * @param isDirect_in  -> specify in which type the coordinates are supplied; default is direct
     * coordinates
     *
     * @note positions is supplied as 3*N array in format x1,y1,z1,x2,y2,z2 where the number
     * refers to the atom index
     */
    void set_positions(const Vec1Real& positions, const bool isDirect_in = true);
    /* set atom types and atom numbers; set atom_type_index array
     *
     * @param atom_types_in      Vector containing string determining the atom types
     * @param number_atoms_in    integer array of same length as atom types containing number of
     * atoms per type
     */
    void set_types(const Vec1String& atom_types_in, const Vec1Int& number_atoms_in);

    /* set the whole poscar structure at once
     * @param lattice_in   Lattice to which lattice has to be set, format xx,xy,xz,yx,yy,yz,zx,zy,zz
     * @param positions    -> 3*N array containing the atomic positions in either direct or
     * cartesian coordinates
     * @param isDirect_in  -> specify in which type the coordinates are supplied; default is direct
     * coordinates
     * @param atom_types_in    ->  Vector containing string determining the atom types
     * @param number_atoms_in  ->  integer array of same length as atom types containing number of
     * atoms per type
     *
     * @note positions is supplied as 3*N array in format x1,y1,z1,x2,y2,z2 where the number
     * refers to the atom index
     */
    void importPoscarData(const Vec1Real&   lattice_in,
                          const Vec1String& atom_types_in,
                          const Vec1Int&    number_atoms_in,
                          const Vec1Real&   positions,
                          const bool        isDirect_in = true);

    /// make cartesian positions from direct
    void directToCartesian();
    /// make direct positions from cartesian
    void cartesianToDirect();
    /// check if coordinates are direct or cartesian
    bool isDirect() const;
    /**
     * get xyz coordinates of atom
     *
     * @param atomIndx index of atom in positions array
     *
     * @return xyz as const tuple
     */
    const std::tuple<const Real&, const Real&, const Real&> getAtom(const Size atomIndx) const;
    /*******************************************************************************************
     * get xyz coordinates of atom
     *
     * @param atomIndx index of atom in positions array
     *
     * @return xyz as rvalue refernce tuple
     *******************************************************************************************/
    std::tuple<Real&&, Real&&, Real&&> getAtom(const Size atomIndx);
    // getters
    /*******************************************************************************************
     * get positions as a const reference
     * for max performance call with const Matrix<Real>& positions = get_positions()
     * to not call copy constructor
     *
     * @return positions 3*N positions array
     *
     * @note format is in x1,y1,z1,x2,y2,z2...
     *******************************************************************************************/
    const Vec1Real& get_positions() const;
    /* get lattice as a const reference
     * for max performance call with const Lattice& lattice = get_lattice()
     * to not call copy constructor
     */
    const Lattice& get_lattice() const;
    /*******************************************************************************************
     * get inverse lattice as a const reference
     *******************************************************************************************/
    const Lattice& get_inverseLattice() const;
    /*******************************************************************************************
     * get atom_types as Integer indices
     *******************************************************************************************/
    const Vec1Int& get_types() const;
    /*******************************************************************************************
     * get atom_types as strings
     *******************************************************************************************/
    const Vec1String& get_typeNames() const;
    /*******************************************************************************************
     * @brief Retrieves the list of atom type names.
     *
     * @return A constant reference to a vector of strings containing the atom type names.
     *******************************************************************************************/
    const Vec1String& get_atomTypeNames() const;
    /*******************************************************************************************
     * @brief Retrieves the atom type name for a given atom index.
     *
     * @param atomIdx The index of the atom whose type name is to be retrieved.
     * @return A constant reference to the atom type name corresponding to the given index.
     *
     * @note The caller must ensure that the provided index is within the valid range of
     *`atomTypeNames`.
     *******************************************************************************************/
    const String& get_atomTypeNames(const Size& atomIdx) const;
    /* get total number of atom in current structure
     */
    const Vec1Int& get_numAtomsPerType() const;
    /*******************************************************************************************
     * get the total number of atoms
     *******************************************************************************************/
    Size get_numAtoms() const;
    /*******************************************************************************************
     * Return indexOffset for batched learning index trick.
     *******************************************************************************************/
    Size get_indexOffset() const;
    /*******************************************************************************************
     * Return typeMap
     *******************************************************************************************/
    std::shared_ptr<TypeMap> get_typeMap() const;
    /*******************************************************************************************
     * Write structure information to screen.
     *******************************************************************************************/
    void writeToScreen() const;
    /*******************************************************************************************
     * Setting indexOffest
     *******************************************************************************************/
    void set_indexOffset(const Size& offset);
    /*******************************************************************************************
     * Setting up typeMap for refit mode.
     *
     * @param superTypes atom types present in the super set of types.
     *******************************************************************************************/
    void set_typeMap(const Vec1String superTypes);
    /*******************************************************************************************
     * @brief Converts the units of the structure's positions and lattice by a given factor.
     *
     * This function rescales the positions and lattice of the structure by the specified factor.
     * If the coordinates are in direct mode, no rescaling of positions is performed, as the units
     * are contained within the lattice. The lattice and its inverse are updated accordingly.
     *
     * @param factor The scaling factor to convert the units.
     *******************************************************************************************/
    void convertUnits(const Real& factor);
    /*******************************************************************************************
     * @brief Sets the atom indices for the structure.
     *
     * @param atomIndx A vector of integers representing the atom indices.
     *******************************************************************************************/
    void set_atomIndx(const Vec1Int& atomIndx);
    /*******************************************************************************************
     * @brief Retrieves the atom indices of the structure.
     *
     * @return A constant reference to a vector of integers representing the atom indices.
     *******************************************************************************************/
    const Vec1Int& get_atomIndx() const;
    /*******************************************************************************************
     * @brief Retrieves a filtered list of local type IDs based on atom indices.
     *
     * @return Vec1Int A vector containing the filtered type IDs.
     *******************************************************************************************/
    Vec1Int get_typesLocal() const;
    /*******************************************************************************************
     * @brief Retrieves a filtered list of local type names based on atom indices.
     *
     * @return Vec1String A vector containing the filtered type names.
     *******************************************************************************************/
    //Vec1String get_atomTypeNamesLocal( void ) const;

  private:
    /// read lattice from poscar
    void read_lattice(std::ifstream& infile_stream);
    /// read positions and types from poscar file
    void read_positions(std::ifstream& infile_stream);
    /// read atom types and numbers
    void read_types(std::ifstream& infile_stream);
    /// shift atoms to primitive cell
    void shift_atoms_primitive_cell();
    /// direct to cartesian for supplied array
    void rescale_coordinates(const Lattice& lattice);
    /**
     * set up an integer array for atom types
     * types can be used to point to
     */
    void build_types();
    /*******************************************************************************************
     * Memory which can be supplied from outside to store data arrays
     *******************************************************************************************/
    ShRec data;
    /// check if positions are in cartesian or direct coordinates
    bool direct;
    /// Lattice vectors of given structure
    Lattice lattice;
    /// inverse lattice vectors of given structure
    Lattice inverseLattice;
    /**
     * positions stored in format x1,y1,z1,x2,y2,z2,x3,y3,z3...
     */
    Vec1Real& positions;
    /// stores atom types
    Vec1String& typeNames;
    /// type index starting at 0 and going up to N number atom types
    Vec1Int&    typeIDs;
    Vec1String& atomTypeNames;
    /// number of atoms per type
    Vec1Int& numAtomsPerType;
    /*******************************************************************************************
     * Storing atom indexes which are MPI distributed. These indexes can be used
     * for neighbor list construction.
     *
     * The index array can either be set up in round robin or block distributed form.
     *******************************************************************************************/
    Vec1Int& atomIndx;
    /// number of atoms
    Size numAtoms;
    /// Scaling factor.
    Real scale;
    /*******************************************************************************************
     * Index offset is only used if the structure is part of a batch.
     *
     * This Marks the indexes of a batched neighbor list with the offset supplie by the variable
     * neighborIndex = neighborIndex + indexOffset. Currently onluy used in Refit mode.
     *******************************************************************************************/
    Size indexOffset = 0;
    /*******************************************************************************************
     * Mapping from structure types to some super set of types. For example supertypes could be
     * from ML_AB.
     *
     * This variable is currently only used in the Refit mode and later in the on the fly
     * training mode.
     *******************************************************************************************/
    std::shared_ptr<TypeMap> typeMap = nullptr;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapStructure();

void convertUnits(Structure& structure, const Real& factor);

} // namespace vaspml

#endif
