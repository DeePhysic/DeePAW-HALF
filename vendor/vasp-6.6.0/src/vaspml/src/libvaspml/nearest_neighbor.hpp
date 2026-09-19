#ifndef NEAREST_NEIGHBOR_HPP
#define NEAREST_NEIGHBOR_HPP

#include "BatchMap.hpp"
#include "Structure.hpp"
#include "types.hpp"

#include <algorithm>
#include <numeric>
#include <tuple>

namespace vaspml
{

class Lattice;

/// Nsquare nearest neighbor algorithm
class NearestNeighborNSquare
{
  public:
    /*******************************************************************************************
     * @class NearestNeighborNSquare
     * @brief A class for performing nearest neighbor calculations using an N-square algorithm.
     *
     * This class provides functionality to find the nearest neighbors of a set of points
     * using an N-square approach, which involves comparing all pairs of points.
     *******************************************************************************************/

    /*******************************************************************************************
     * @brief Default constructor for the NearestNeighborNSquare class.
     *
     * This constructor initializes the NearestNeighborNSquare object with default parameters.
     *
     * @param record Optional shared record object (default is nullptr).
     *               This can be used to store or retrieve additional information during
     *               the nearest neighbor calculations.
     *******************************************************************************************/
    NearestNeighborNSquare(ShRec record = nullptr);
    /*******************************************************************************************
     * @brief Parameterized constructor for the NearestNeighborNSquare class.
     *
     * This constructor initializes the NearestNeighborNSquare object with user-defined parameters.
     *
     * @param Cutoff The cutoff distance for nearest neighbor calculations.
     *               Points separated by a distance greater than this value will not be considered
     *neighbors.
     * @param isTypeSorted_in A boolean flag indicating whether the input data is type-sorted.
     *                        If true, the data is assumed to be sorted by type, which may optimize
     *calculations.
     * @param isDistSorted_in A boolean flag indicating whether the input data is distance-sorted.
     *                        If true, the data is assumed to be sorted by distance, which may
     *optimize calculations.
     * @param record Optional shared record object (default is nullptr).
     *               This can be used to store or retrieve additional information during
     *               the nearest neighbor calculations.
     *******************************************************************************************/
    NearestNeighborNSquare(Real  Cutoff,
                           bool  isTypeSorted_in = false,
                           bool  isDistSorted_in = false,
                           ShRec record = nullptr);
    void init(const Real Cutoff,
              const bool isTypeSorted_in = false,
              const bool isDistSorted_in = false);
    /** compute nearest neighbor arrays in terms of direct coordinates
     *  the results, distance vectors, distance, normalized distance
     *  will still be in cartesian coordinates
     *
     * @param positions strucure containing information as positions, lattice, atom types
     */
    void computeNearestNeighborsDirectCoordinates(const Structure& positions);
    /*******************************************************************************************
     * @brief Computes the nearest neighbors for a set of central atoms using direct coordinates.
     *
     * This function calculates the nearest neighbors for the specified central atoms in the given
     *structure. It assumes that the positions are provided in direct coordinates. The function
     *allocates necessary storage, extracts atom types for the central atoms, computes unique types,
     *and processes each central atom to determine its nearest neighbors.
     *
     * @param positions The structure containing atomic positions and types. Must be in direct
     *coordinates.
     * @param centralAtoms A vector of indices specifying the central atoms for which nearest
     *neighbors are computed.
     *
     * @note The function will throw an error if the positions are not in direct coordinates.
     *******************************************************************************************/
    void computeNearestNeighborsDirectCoordinates(const Structure& positions,
                                                  const Vec1Int    centralAtoms);
    /*******************************************************************************************
     * @brief Computes neighbors for a batch of atoms using direct coordinates.
     *
     * This function processes a range of structures and computes the neighbors for each atom
     * in the batch using their direct coordinates. The batch of atoms is specified using the
     * provided `batchMap`.
     *
     * @tparam Iterator The type of the iterator used to traverse the structures. Typically, this
     *         will be an iterator over a container of structures.
     * @param strucBegin Iterator pointing to the beginning of the range of structures to process.
     * @param strucEnd Iterator pointing to the end of the range of structures to process.
     * @param batchMap A mapping of atom batches that specifies which atoms are part of the batch
     *        for which neighbors need to be computed. This is typically a data structure that
     *        associates atoms with their batch identifiers.
     *
     * @note This function assumes that the structures and batch map are properly initialized
     *       and valid. It does not perform any validation of the input data.
     *
     * @warning The behavior of this function is undefined if the iterators do not form a valid
     *range or if the `batchMap` contains invalid data.
     *
     * @pre The range `[strucBegin, strucEnd)` must be valid and contain structures with atoms
     *      that can be processed.
     * @pre The `batchMap` must be properly initialized and contain valid mappings for the atoms.
     *
     * @post Neighbors for the specified batch of atoms are computed and stored in the appropriate
     *       data structures.
     *
     * @example
     * @code
     * std::vector<Structure> structures = loadStructures();
     * AtomBatchMap batchMap = createBatchMap();
     * computeNeighborsDirectCoordinatesBatch(structures.begin(), structures.end(), batchMap);
     * @endcode
     *******************************************************************************************/
    template<typename Iterator>
    void computeNeighborsDirectCoordinatesBatch(const Iterator      strucBegin,
                                                const Iterator      strucEnd,
                                                const AtomBatchMap& batchMap);
    /** compute nearest neighbor arrays in terms of Cartesian coordinates
     *  the results, distance vectors, distance, normalized distance
     *  will be in cartesian coordinates
     *
     * @param positions strucure containing information as positions, lattice, atom types
     */
    void computeNearestNeighborsCartesianCoordinates(const Structure& positions);
    /*******************************************************************************************
     * @brief Computes the nearest neighbors for a set of central atoms in Cartesian coordinates.
     *
     * This function calculates the nearest neighbors for a given set of central atoms based on
     *their Cartesian coordinates. It ensures that the positions are in direct coordinates,
     *allocates necessary storage arrays, determines the unique types of central atoms, and computes
     *the nearest neighbors for each central atom. Additionally, it computes the mapping of central
     *atom indices per type.
     *
     * @param positions The structure containing atomic positions and related information.
     *                  It must be in direct coordinates.
     * @param centralAtoms A vector of indices representing the central atoms for which nearest
     *neighbors are to be computed.
     *
     * @note The function assumes that the input positions are in direct coordinates. If not, an
     *error is raised.
     *******************************************************************************************/
    void computeNearestNeighborsCartesianCoordinates(const Structure& positions,
                                                     const Vec1Int    centralAtoms);
    /*******************************************************************************************
     * @brief Computes the nearest Cartesian coordinates for a batch of structures.
     *
     * This function processes a range of structures and computes the nearest Cartesian coordinates
     * for each structure in the batch. The mapping between atoms and batches is provided by the
     * `batchMap` parameter.
     *
     * @tparam Iterator The type of the iterator used to traverse the range of structures.
     *Typically, this would be an iterator for a container holding the structures.
     *
     * @param strucBegin Iterator pointing to the beginning of the range of structures to process.
     * @param strucEnd Iterator pointing to the end of the range of structures to process.
     * @param batchMap A mapping of atom batches, represented by an `AtomBatchMap` object. This map
     *                 provides the necessary information to compute the nearest Cartesian
     *coordinates for each batch.
     *
     * @note The function assumes that the range defined by `strucBegin` and `strucEnd` is valid and
     *       that `batchMap` contains the required data for processing.
     *
     * @warning Ensure that the iterators and the `batchMap` are properly initialized before calling
     *          this function to avoid undefined behavior.
     *
     * @pre The range `[strucBegin, strucEnd)` must be valid and non-empty.
     * @pre `batchMap` must contain valid mappings for the atom batches.
     *
     * @post The nearest Cartesian coordinates for each structure in the batch are computed and
     *stored as per the implementation details.
     *
     * @throws std::runtime_error If an error occurs during computation or if the input data is
     *invalid.
     *
     * @example
     * @code
     * std::vector<Structure> structures = { ... }; // A container of structures
     * AtomBatchMap batchMap = { ... };            // A mapping of atom batches
     * computeNearestCartesianCoordinatesBatch(structures.begin(), structures.end(), batchMap);
     * @endcode
     *******************************************************************************************/
    template<typename Iterator>
    void computeNeighborsCartesianCoordinatesBatch(const Iterator      strucBegin,
                                                   const Iterator      strucEnd,
                                                   const AtomBatchMap& batchMap);
    /// write nearest neighbor list to screen
    void writeListToScreen() const;
    /*******************************************************************************************
     * get global index of neighbor in position array
     * @param atomIndx    central atom of which a neigbor is requested
     * @param indxNeigh   number of neighor in nearest neighbor array to return
     *******************************************************************************************/
    const Int& get_globalIndex(const Size atomIndx, const Size indxNeigh) const;
    /**  get whole global index array of whole nearest neighbor list
     */
    const Vec2Int& get_globalIndex() const;
    /**  get global index array for certain central atom
     *
     * @param atomIndx index of central atom for which global index array is requested
     */
    const Vec1Int& get_globalIndex(const Size atomIndx) const;
    /** get distance between atoms atomIndx and atom get_index_list( atomIndx, inxdNeigh )
     *  distances are ordered in the same way as index_list
     *@param atomIndx   -> central atom from which a distance to
     *@param indxNeigh  -> points to the list entry of the neighbor
     */
    const Real& get_distances(const Size atomIndx, const Size indxNeigh) const;
    /**  get all distances between atoms in neighbor list
     *
     * @note format of return array is [central_atom][ vector of neighbor distances ]
     */
    const Vec2Real& get_distances() const;
    /**
     * get distances of nearest neighbors for given central atom
     *
     * @param atomIndx
     */
    const Vec1Real& get_distances(Size atomIndx) const;
    /**
     * get x component of distance vector between atom pair
     *
     * distances are ordered in the same way as index_list
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVector_x(const Size atomIndx, Size indxNeigh) const;
    /**  get y component of distance vector between atoms atomIndx and atom get_index_list(
     * atomIndx, inxdNeigh ) distances are ordered in the same way as index_list
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVector_y(const Size atomIndx, Size indxNeigh) const;

    /**  get z componeent of distance vector between atoms atomIndx and atom get_index_list(
     *atomIndx, inxdNeigh ) distances are ordered in the same way as index_list
     *@param atomIndx   -> central atom from which a distance to
     *@param indxNeigh  -> points to the list entry of the neighbor
     */
    const Real& get_connectionVector_z(const Size atomIndx, Size indxNeigh) const;
    /**
     * get x,y,z component of distance vector between atom pair
     *
     * distances are ordered in the same way as index_list for optimal
     * performance call as auto & [ x, y, z ] = get_connectionVector( atomIndx, indxNeigh )
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const std::tuple<const Real&, const Real&, const Real&> get_connectionVector(
        Size atomIndx,
        Size indxNeigh) const;
    /**  get x,y,z component of distance vector between all atoms
     *
     * first index is the central atom
     * format is xyz atom1 neighbor 1, xyz atom1 neighbor 2,...., xyz atom2, neighbor1...
     */
    const Vec2Real& get_connectionVector() const;
    /**
     * get connection vectors for cetrain central atom defined by atomIndx
     *
     * data is of format xyz neighbo1, xyz neighbor2, xyz neighbor3,...
     * @param atomIndx list index of central atom
     */
    const Vec1Real& get_connectionVector(Size atomIndx) const;
    /** get x component of normalized distance vector between atom pair
     *
     *  distances are ordered in the same way as index_list
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_x(const Size atomIndx, Size indxNeigh) const;
    /**  get y component of normalized distance vector between atom pairs
     *
     * distances are ordered in the same way as index_list
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_y(const Size atomIndx, Size indxNeigh) const;
    /**  get z component of normalized distance vector between atom pairs
     *
     *   distances are ordered in the same way as index_list
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_z(const Size atomIndx, Size indxNeigh) const;
    /**  get x,y,z component of distance vector between atom pairs
     *
     * distances are ordered in the same way as index_list for optimal
     * performance call as auto & [ x, y, z ] = get_connectionVector( atomIndx, indxNeigh )
     *
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const std::tuple<const Real&, const Real&, const Real&> get_connectionVectorNormalized(
        Size atomIndx,
        Size indxNeigh) const;
    /**  get x,y,z component of distance vector between all atoms
     *  format is xyz atom1 neighbor 1, xyz atom1 neighbor 2,...., xyz atom2, neighbor1...
     */
    const Vec2Real& get_connectionVectorNormalized() const;
    /**
     * get normalized connection vectors for certain central atom
     *
     * @param atomIndx atom index of desired central atom
     */
    const Vec1Real& get_connectionVectorNormalized(Size atomIndx) const;
    /**
     * get number of nearest neighbors for atom atomIndx
     *@param atomIndx central atom index
     */
    Size get_size(const Size atomIndx) const;
    /**
     * get number of nearest neighbors whole array
     */
    const Vec1Int& get_size() const;
    /*******************************************************************************************
     * get type index of central atom in nearest neighbor list
     *
     * @param atomIndex of central atom from which type has to be retrieved
     *******************************************************************************************/
    const Int& get_typeIndexCentral(const Size atomIndex) const;
    /*******************************************************************************************
     * @brief Retrieves the type key central for a specified atom index.
     *
     * This function returns a reference to the type key central associated with the given atom
     *index. The type key central is typically used to identify or categorize the atom in a specific
     *context.
     *
     * @param atomIndex The index of the atom for which the type key central is requested.
     *                  This should be a valid index within the bounds of the atom collection.
     * @return A constant reference to the type key central associated with the specified atom
     *index. The returned value is guaranteed to remain valid as long as the object owning this
     *method is not modified or destroyed.
     *
     * @note Ensure that the provided `atomIndex` is within the valid range to avoid undefined
     *behavior.
     *******************************************************************************************/
    const String& get_typeKeyCentral(const Size atomIndex) const;
    /*******************************************************************************************
     * get types of all central atoms
     *
     * for performance call with const Vec1Int& typeIndexCentral = get_typeIndexCentral();
     *******************************************************************************************/
    const Vec1Int& get_typeIndexCentral() const;
    /*******************************************************************************************
     * @brief Retrieves the central type key.
     *
     * This function returns a constant reference to the `typeKeyCentral` member,
     * which represents the central type key in the system. The returned value
     * is immutable and cannot be modified.
     *
     * @return A constant reference to a `Vec1String` object representing the central type key.
     *******************************************************************************************/
    const Vec1String& get_typeKeyCentral() const;
    /**
     * get type index of neighbor atom in nearest neighbor list
     *@param atomIndex central atom index from which neighbor the types is retrieved
     *@param indxNeigh number of neighbor atom from which type is wanted
     */
    const Int& get_typeIndex(const Size atomIndex, const Size indxNeigh) const;
    /*******************************************************************************************
     * @brief Retrieves the type key associated with a specific atom and its neighbor.
     *
     * This function returns a reference to the type key for a given atom and its neighbor
     * based on their indices. The type key is typically used to identify or classify
     * the relationship or interaction between the atom and its neighbor.
     *
     * @param atomIndex The index of the atom for which the type key is being retrieved.
     * @param indxNeigh The index of the neighbor associated with the atom.
     * @return A constant reference to the type key as a String.
     *
     * @note The caller must ensure that the indices provided are valid and within bounds.
     *       Behavior is undefined if the indices are out of range.
     *******************************************************************************************/
    const String& get_typeKey(const Size atomIndex, const Size indxNeigh) const;
    /*******************************************************************************************
     * get type index of neighbors for given central atom
     * @param atomIndex central atom index from which neighbor the types is retrieved
     *******************************************************************************************/
    const Vec1Int& get_typeIndex(const Size atomIndex) const;
    /*******************************************************************************************
     * @brief Retrieves the type key associated with a specific atom index.
     *
     * This function returns a reference to a `Vec1String` object that represents
     * the type key for the atom specified by the given index. The type key is
     * typically used to identify or categorize the atom in a molecular or structural
     * context.
     *
     * @param atomIndex The index of the atom for which the type key is requested.
     *                  This should be a valid index within the range of atoms managed
     *                  by the class.
     * @return A constant reference to a `Vec1String` object representing the type key
     *         of the specified atom.
     *
     * @note The caller must ensure that `atomIndex` is within the valid range of indices.
     *       Accessing an invalid index may result in undefined behavior.
     *******************************************************************************************/
    const Vec1String& get_typeKey(const Size atomIndex) const;
    /**
     * get types of all neighbor atoms
     *
     * for performance call with const Vec2Int& typeIndex = get_typeIndex();
     */
    const Vec2Int& get_typeIndex() const;
    /*******************************************************************************************
     * @brief Retrieves the type key associated with the object.
     *
     * This function returns a constant reference to the type key, which is a
     * `Vec2String` object. The type key typically represents a unique identifier
     * or classification for the object.
     *
     * @return A constant reference to the `Vec2String` representing the type key.
     *******************************************************************************************/
    const Vec2String& get_typeKey() const;
    /**
     * get total number of atoms for types
     *
     */
    const Vec1Int& get_nAtomsType() const;
    /**
     * get total number of atoms for certain type
     *@param type index of type for which atom number should be obtained
     */
    const Int& get_nAtomsType(const Size type) const;
    /// get total number of atoms
    const Size& get_nAtoms() const;
    /*******************************************************************************************
     * @brief Retrieves the number of neighbors for all atoms.
     *
     * This function returns a constant reference to a `Vec1Int` object that contains
     * the number of neighbors for all atoms in the system. The `Vec1Int` is assumed
     * to be a container (e.g., vector or array) holding integer values, where each
     * value corresponds to the number of neighbors for a specific atom.
     *
     * @return A constant reference to a `Vec1Int` object representing the number of neighbors for
     *all atoms.
     *******************************************************************************************/
    const Vec1Int& get_numberNeighbors() const;
    /*******************************************************************************************
     * @brief Retrieves the number of neighbors for a specific atom.
     *
     * This function returns the number of neighbors for the atom specified by its index.
     * The `atomIdx` parameter is used to identify the atom within the system, and the
     * function returns a constant reference to an integer representing the number of neighbors.
     *
     * @param atomIdx The index of the atom for which the number of neighbors is requested.
     *                It is of type `Size`, which is typically an unsigned integer or similar type.
     * @return A constant reference to an `Int` object representing the number of neighbors for the
     *specified atom.
     * @throws std::out_of_range If the provided `atomIdx` is invalid or exceeds the bounds of the
     *atom list.
     *******************************************************************************************/
    const Int& get_numberNeighbors(const Size atomIdx) const;
    /// number of nearest neighbors per type for central atom
    const Vec2Int& get_numberNeighborsType() const;
    /*******************************************************************************************
     * @brief Retrieves the number of neighbors of a specific type for a given atom.
     *
     * This function returns a constant reference to a `Vec1Int` object that represents
     * the number of neighbors of a specific type for the atom identified by the given index.
     *
     * @param atomIdx The index of the atom for which the number of neighbors of a specific type is
     *requested. This should be a valid index within the range of atoms in the system.
     * @return A constant reference to a `Vec1Int` object containing the number of neighbors of the
     *specified type.
     *
     * @note The caller must ensure that `atomIdx` is within the valid range of atom indices.
     * @warning Accessing an invalid `atomIdx` may result in undefined behavior.
     *******************************************************************************************/
    const Vec1Int& get_numberNeighborsType(const Size atomIdx) const;
    /**
     * get all nearest neighbor data from one atom pair
     *
     * call function as auto &[ type, r, x, y, z, x_norm, y_norm, z_norm ] = get_neighborData(
     * atomIndx, neighIndx )
     * @param atomIndex integer defining the central atom
     * @param indxNeigh integer defining entry of neighbor atom
     */
    const std::tuple<const Int&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&>
    get_neighborData(const Size atomIndex, const Size indxNeigh) const;
    /**
     * return the number of atoms in current neighbor list
     */
    Size get_numberAtoms() const;
    /**
     *
     * get the number of different types in current structure
     */
    Size get_numberTypes() const;
    /**
     * getter of control variable if neighbor list is type sorted
     */
    bool is_typeSorted() const;
    /*******************************************************************************************
     * get array which stores the index of a central atom per type
     *******************************************************************************************/
    const Vec1Int& get_centralAtomIndexPerType() const;
    /*******************************************************************************************
     * @brief Retrieves the central atom index for a specific atom type.
     *
     * This function returns a reference to the central atom index associated with a given atom
     *type, identified by its index. The central atom index is typically used in molecular modeling
     *or computational chemistry to represent the key atom in a group or structure.
     *
     * @param atomIdx The index of the atom type for which the central atom index is requested.
     *                This should be a valid index within the range of atom types.
     * @return A constant reference to the central atom index corresponding to the specified atom
     *type.
     *
     * @note The caller must ensure that the provided `atomIdx` is valid and within bounds.
     *       Accessing an invalid index may result in undefined behavior.
     *
     * @warning The returned reference is constant, meaning the value cannot be modified directly.
     *******************************************************************************************/
    const Int& get_centralAtomIndexPerType(const Size atomIdx) const;
    /**
     *returns the latticeVolume
     */
    const Real& get_latticeVolume() const;
    /// Get cutoff radius used for this neighbor list.
    Real get_cutOff() const;
    /*******************************************************************************************
     * @brief Checks if the type data is sorted.
     *
     * This function returns a boolean value indicating whether the type data
     * is sorted in the current context. Sorting may depend on specific criteria
     * defined elsewhere in the implementation.
     *
     * @return `true` if the type data is sorted, `false` otherwise.
     *******************************************************************************************/
    bool get_isTypeSorted() const;
    /*******************************************************************************************
     * @brief Checks if the distance data is sorted.
     *
     * This function returns a boolean value indicating whether the distance data
     * is sorted in the current context. Sorting may depend on specific criteria
     * defined elsewhere in the implementation.
     *
     * @return `true` if the distance data is sorted, `false` otherwise.
     *******************************************************************************************/
    bool get_isDistSorted() const;
    /*******************************************************************************************
     * setting the total number of atoms central atoms in the neighbor.
     *
     * @param nAtoms number of central atoms in current neighbor list
     *
     * @note this variable has to be set when the neighbor list is filled via an interface
     * for example to VASP. If the neighbor list is computed with it's own routines
     * the variable will be automatically set.
     *******************************************************************************************/
    void set_nAtoms(const Int nAtoms);
    /*******************************************************************************************
     * setting the total number of atom types in neighbor list
     *
     * @param nTypes number of unique atom types in neighbor list
     *
     * @note this variable has to be set when the neighbor list is filled via an interface
     * for example to VASP. If the neighbor list is computed with it's own routines
     * the variable will be automatically set.
     *******************************************************************************************/
    void set_nTypes(const Int nTypes);
    void set_totalNumberTypes(const Int nTypes);
    /*******************************************************************************************
     * setting up the centralAtomIndexPerType array
     *
     * Array stores for given global atom index the number of the atom within it's type
     *******************************************************************************************/
    void compute_centralAtomIndexPerType();
    /*******************************************************************************************
     * @brief Resets neighbor-related arrays with the provided data.
     *
     * @param globalIndex Global indices of atoms.
     * @param typeIndex Type indices of atoms.
     * @param typeIndexCentral Type index of the central atom.
     * @param distances Distances between atoms.
     * @param connectionVector Connection vectors between atoms.
     * @param connectionVectorNormalized Normalized connection vectors.
     * @param numberNeighbors Number of neighbors for each atom.
     * @param numberNeighborsType Number of neighbors per type for each atom.
     * @param centralAtomIndexPerType Central atom indices per type.
     *******************************************************************************************/
    void resetNeighborArrays(const ShVec2Int&  globalIndex,
                             const ShVec2Int&  typeIndex,
                             const ShVec1Int&  typeIndexCentral,
                             const ShVec2Real& distances,
                             const ShVec2Real& connectionVector,
                             const ShVec2Real& connectionVectorNormalized,
                             const ShVec1Int&  numberNeighbors,
                             const ShVec2Int&  numberNeighborsType,
                             const ShVec1Int&  centralAtomIndexPerType);
    /*******************************************************************************************
     * @brief Detaches neighbor-related arrays, releasing any associated resources.
     *******************************************************************************************/
    void detachNeighborArrays();
    /*******************************************************************************************
     * @brief Allocates arrays for reference configuration based on the provided sizes.
     *
     * @param sizes Sizes of the arrays to be allocated.
     *******************************************************************************************/
    void allocateArraysRefConfWrapper(const Vec1Size& sizes);
    /*******************************************************************************************
     * @brief Reorders the neighbor list using the provided batch map.
     *
     * @param batchMap Mapping of atom batches for reordering.
     *******************************************************************************************/
    void reorderNeighborList(const Vec1Int& newOrder);
    /*******************************************************************************************
     * @brief Sets parameters for neighbor list computation.
     *
     * @param Cutoff Cutoff distance for neighbor detection.
     * @param isTypeSorted_in Flag indicating if neighbors are sorted by type.
     * @param isDistSorted_in Flag indicating if neighbors are sorted by distance.
     *******************************************************************************************/
    void setParameters(Real Cutoff, bool isTypeSorted_in, bool isDistSorted_in);
    /*******************************************************************************************
     * @brief Sets data for neighbor-related arrays.
     *
     * @param globalIndex Global indices of atoms.
     * @param typeIndex Type indices of atoms.
     * @param typeKey Type keys of atoms.
     * @param typeIndexCentral Type index of the central atom.
     * @param typeKeyCentral Type key of the central atom.
     * @param distances Distances between atoms.
     * @param connectionVector Connection vectors between atoms.
     * @param connectionVectorNormalized Normalized connection vectors.
     * @param numberNeighbors Number of neighbors for each atom.
     * @param numberNeighborsType Number of neighbors per type for each atom.
     * @param centralAtomIndexPerType Central atom indices per type.
     *******************************************************************************************/
    void setData(const Vec2Int&    globalIndex,
                 const Vec2Int&    typeIndex,
                 const Vec2String& typeKey,
                 const Vec1Int&    typeIndexCentral,
                 const Vec1String& typeKeyCentral,
                 const Vec2Real&   distances,
                 const Vec2Real&   connectionVector,
                 const Vec2Real&   connectionVectorNormalized,
                 const Vec1Int&    numberNeighbors,
                 const Vec2Int&    numberNeighborsType,
                 const Vec1Int&    centralAtomIndexPerType);
    /*******************************************************************************************
     * @brief Extends existing neighbor-related data with additional entries.
     *
     * @param globalIndex Global indices of atoms.
     * @param typeIndex Type indices of atoms.
     * @param typeKey Type keys of atoms.
     * @param typeIndexCentral Type index of the central atom.
     * @param typeKeyCentral Type key of the central atom.
     * @param distances Distances between atoms.
     * @param connectionVector Connection vectors between atoms.
     * @param connectionVectorNormalized Normalized connection vectors.
     * @param numberNeighbors Number of neighbors for each atom.
     * @param numberNeighborsType Number of neighbors per type for each atom.
     * @param centralAtomIndexPerType Central atom indices per type.
     *******************************************************************************************/
    void extendData(const Vec2Int&    globalIndex,
                    const Vec2Int&    typeIndex,
                    const Vec2String& typeKey,
                    const Vec1Int&    typeIndexCentral,
                    const Vec1String& typeKeyCentral,
                    const Vec2Real&   distances,
                    const Vec2Real&   connectionVector,
                    const Vec2Real&   connectionVectorNormalized,
                    const Vec1Int&    numberNeighbors,
                    const Vec2Int&    numberNeighborsType,
                    const Vec1Int&    centralAtomIndexPerType);
    /*******************************************************************************************
     * @brief Adds a single element to the neighbor-related data.
     *
     * @param globalIndex Global index of the atom.
     * @param typeIndex Type index of the atom.
     * @param typeKey Type key of the atom.
     * @param typeIndexCentral Type index of the central atom.
     * @param typeKeyCentral Type key of the central atom.
     * @param distances Distance to the neighbor.
     * @param connectionVector Connection vector to the neighbor.
     * @param connectionVectorNormalized Normalized connection vector to the neighbor.
     * @param numberNeighbors Number of neighbors for the atom.
     * @param numberNeighborsType Number of neighbors per type for the atom.
     * @param centralAtomIndexPerType Central atom index per type.
     *******************************************************************************************/
    void addElement(const Vec1Int&    globalIndex,
                    const Vec1Int&    typeIndex,
                    const Vec1String& typeKey,
                    const Int&        typeIndexCentral,
                    const String&     typeKeyCentral,
                    const Vec1Real&   distances,
                    const Vec1Real&   connectionVector,
                    const Vec1Real&   connectionVectorNormalized,
                    const Int&        numberNeighbors,
                    const Vec1Int&    numberNeighborsType,
                    const Int&        centralAtomIndexPerType);
    /*******************************************************************************************
     * @brief Computes the number of atoms per type based on unique type keys.
     *
     * @param unique List of unique type keys.
     *******************************************************************************************/
    void computeNumberAtomsPerType(const Vec1String& unique);
    /*******************************************************************************************
     * @brief Overloads the less-than operator for comparing two NearestNeighborNSquare objects.
     *
     * @param other The other NearestNeighborNSquare object to compare with.
     * @return true if the cutOff value of this object is less than that of the other object.
     * @return false otherwise.
     *******************************************************************************************/
    bool operator<(const NearestNeighborNSquare& other) const;
    /*******************************************************************************************
     * @brief Overloads the greater-than operator for comparing two NearestNeighborNSquare objects.
     *
     * @param other The other NearestNeighborNSquare object to compare with.
     * @return true if the cutOff value of this object is greater than that of the other object.
     * @return false otherwise.
     *******************************************************************************************/
    bool operator>(const NearestNeighborNSquare& other) const;
    /*******************************************************************************************
     * @brief Overloads the less-than-or-equal-to operator for comparing two NearestNeighborNSquare
     *objects.
     *
     * @param other The other NearestNeighborNSquare object to compare with.
     * @return true if the cutOff value of this object is less than or equal to that of the other
     *object.
     * @return false otherwise.
     *******************************************************************************************/
    bool operator<=(const NearestNeighborNSquare& other) const;
    /*******************************************************************************************
     * @brief Overloads the equality operator for comparing two NearestNeighborNSquare objects.
     *
     * @param other The other NearestNeighborNSquare object to compare with.
     * @return true if the cutOff value of this object is equal to that of the other object.
     * @return false otherwise.
     *******************************************************************************************/
    bool operator==(const NearestNeighborNSquare& other) const;
    /*******************************************************************************************
     * @brief Overloads the inequality operator for comparing two NearestNeighborNSquare objects.
     *
     * @param other The other NearestNeighborNSquare object to compare with.
     * @return true if the cutOff value of this object is not equal to that of the other object.
     * @return false otherwise.
     *******************************************************************************************/
    bool operator!=(const NearestNeighborNSquare& other) const;
    /*******************************************************************************************
     * @brief Computes the number of neighbors per type for each atom.
     *
     * This function calculates a 2D vector where each row corresponds to an atom,
     * and each column represents the count of neighbors of a specific type for that atom.
     *
     * @return Vec2Int A 2D vector where `numberNeighborsPerType[atom][type]`
     *         contains the count of neighbors of `type` for the given `atom`.
     *******************************************************************************************/
    Vec2Int computeNumberNeighborsPerType() const;
    /*******************************************************************************************
     * @brief Computes the maximum number of neighbors per type across all atoms.
     *
     * This function calculates a 1D vector where each element corresponds to the
     * maximum number of neighbors of a specific type across all atoms.
     *
     * @return Vec1Int A 1D vector where `maxNeighborsPerType[type]` contains the
     *         maximum count of neighbors of `type` across all atoms.
     *******************************************************************************************/
    Vec1Int      computeMaxNeighborsPerType() const;
    const ShRec& getNeighborRecord() const;
    ShRec&       getNeighborRecord();
    Vec1String   getUniqueCentralTypes() const;

  private:
    /// check for unique atom types
    void compute_unique_types(const Vec1Int& types);
    /** compute nearest neighbor data for a single atom starting from
     *  direct coordinates, results will be given in Cartesian coordinates
     *
     * @param indx  unsigned integer denoting the actual central atom
     * @param entryIndex position to which data will be stored
     * @param positions structure containing containing information about atomic structure
     */
    void computeNearestNeighborsSingleAtomDirect(const Size       indx,
                                                 const Size       entryIndex,
                                                 const Structure& positions);
    /** compute nearest neighbor data for a single atom starting from
     *  Cartesian coordinates, results will be given in Cartesian coordinates
     *
     * @param indx  unsigned integer denoting the actual central atom
     * @param entryIndex position to which data will be stored
     * @param positions structure containing containing information about atomic structure
     */
    void computeNearestNeighborsSingleAtomCartesian(const Size       indx,
                                                    const Size       entryIndex,
                                                    const Structure& positions);
    /** allocate work arrays for computation of nearest neighbors single atom
     *
     * @param nn_list      temporary nearest (N) neighbor list for single atom
     * @param types_index  temporary integer (N) array containing atom types associated to nn_list
     * @param dist         temporary real (N) array containing distances in same order as nn_list
     * @param dist_vecs    temporary real (3*N) array containing distance vectors in same order as
     * nn_list first dimension stores x,y,z
     * @param normed_dist  temporary real (3*N) array containing nromalized distance vectors
     *                     in same order as nn_list first dimension stores x,y,z
     */
    void allocateWorkArrays(Vec1Int&    nn_list,
                            Vec1Int&    types_index,
                            Vec1String& types_key,
                            Vec1Real&   dist,
                            Vec1Real&   dist_vecs,
                            Vec1Real&   normed_dist,
                            Size        N);
    /** fill the nearest neighbor arrays as defined at VARIABLES of this object
     *  for every atom; index defines the central atom
     *
     * @param index        Size storing the index of the actual central atom
     * @param entryIndex   determines at which position data is stored in neighbor arrays
     * @param nelements    Size denoting how many nearest neighbors have been found
     * @param nn_list      temporary nearest (N) neighbor list for single atom
     * @param types_index  temporary integer (N) array containing atom types associated to nn_list
     * @param dist         temporary real (N) array containing distances in same order as nn_list
     * @param dist_vecs    temporary real (3*N) array containing distance vectors in same order as
     * nn_list first dimension stores x,y,z
     * @param normed_dist  temporary real (3*N) array containing nromalized distance vectors
     *                     in same order as nn_list first dimension stores x,y,z
     */
    void fillNeighborArrays(const Size        index,
                            const Size        entryIndex,
                            const Size        nelements,
                            const Vec1Int&    nn_list,
                            const Vec1Int&    types_index,
                            const Vec1String& types_key,
                            const Vec1Real&   dist,
                            const Vec1Real&   dist_vecs,
                            const Vec1Real&   normed_dist);
    /**
     * compute how many periodic boxes have to be searched for
     * nearest neighbors of the given central atom
     *
     * @param pos 3d-vector containing Cartesian coordinates of actual central atom
     * @param inverse_lattice inverse of the supplied bravais lattice
     * @param xmin boxes in negative x,y,z direction that have to be chcked for nearest neighbors
     * @param xmax boxes in positive x,y,z direction that have to be chcked for nearest neighbors
     */
    void set_periodic_size_images(const Vec1Real& pos_central,
                                  const Lattice&  inverse_lattice,
                                  Vec1Int&        xmin,
                                  Vec1Int&        xmax);
    /**
     * allocate storage arrays
     *@param size primary size of nearest neighbor array
     * number of central atoms in nearest neighbor list
     */
    void allocateStorageArrays(Size size);

    /*******************************************************************************************
     * Record object managing memory of class members
     *******************************************************************************************/
    ShRec data;
    /// nearest neighbor array containing index to position array
    Vec2Int& globalIndex;
    /// Type index of neighbor atom
    Vec2Int& typeIndex;
    /*******************************************************************************************
     * String array determining the type of the central atom
     *******************************************************************************************/
    Vec2String typeKey;
    /// Type index of central atom
    Vec1Int& typeIndexCentral;
    /*******************************************************************************************
     * String array determining the type of the central atom
     *******************************************************************************************/
    Vec1String& typeKeyCentral;
    /// nearest neighbor distances
    Vec2Real& distances;
    /// nearest neighbor connection vectors
    Vec2Real& connectionVector;
    /// nearest neighbor normalized connection vectors
    Vec2Real& connectionVectorNormalized;
    /// number of nearest neighbors per central atom
    Vec1Int& numberNeighbors;
    /// number of nearest neighbors per type for central atom
    Vec2Int& numberNeighborsType;
    /// total number of atoms per type
    Vec1Int& nAtomsType;
    /*******************************************************************************************
     * gives for for global central atom index the atom number within a type
     *******************************************************************************************/
    Vec1Int& centralAtomIndexPerType;
    /*******************************************************************************************
     * key for neighbor atom type
     *******************************************************************************************/
    Vec2String neighborAtomTypeKey;
    // parameters needed for computation
    /// cutoff radius in units of input coordinates
    Real cutOff = (Real)0.0;
    /// squared cutoff radius
    Real cutOffSquared = (Real)0.0;
    /// minimum distance atoms have to be apart
    /// to avoid self neighbors
    Real eps;
    /// point to array start indices for every atom type / both start and end point to zero
    /// if no atoms of type are present
    ShVec2Int typeStart;
    /// point to array end indices for every atom type
    ShVec2Int typeEnd;
    /// point to array start indices for every atom type
    ShVec2Int neighborType;
    /// unique types in supplied structure
    std::map<Int, Int> unique_types;
    /// unique types in supplied structure
    /// number of different types in supplied atoms structure
    Int nTypes = 0;
    Int totalNumberTypes = 0;
    /// total number of atoms
    Size nAtoms = 0;
    /// decides the storage order of the nearest neighbor array
    /// sorts nearest neighbors according to types
    bool isTypeSorted = false;
    /// decides the storage order of the nearest neighbor array
    /// sorts nearest neighbors according to distances
    bool isDistSorted = false;
    /**
     * volume of box of atoms from which neighbor list is computed
     *
     *the volume is needed for the stress tensor compuation in
     *DescriptorSHS2
     */
    Real latticeVolume = (Real)0;
};

template<typename Iterator>
void NearestNeighborNSquare::computeNeighborsDirectCoordinatesBatch(const Iterator      strucBegin,
                                                                    const Iterator      strucEnd,
                                                                    const AtomBatchMap& batchMap)
{
    nAtoms =
        std::accumulate(strucBegin,
                        strucEnd,
                        0,
                        [](Int sum, const Structure& struc) { return sum + struc.get_numAtoms(); });
    allocateStorageArrays(nAtoms);
    const Vec1Int&    typeIdx = batchMap.get_typesIdx();
    const Vec1String& type = batchMap.get_types();

    std::copy(typeIdx.cbegin(), typeIdx.cend(), typeIndexCentral.begin());
    std::copy(type.cbegin(), type.cend(), typeKeyCentral.begin());
    computeNumberAtomsPerType(batchMap.get_typeOrder());

    std::for_each(strucBegin,
                  strucEnd,
                  [&](const Structure& struc)
                  {
                      for (Size i = 0; i < struc.get_numAtoms(); i++)
                      {
                          computeNearestNeighborsSingleAtomDirect(i, i, struc);
                      }
                  });
    compute_centralAtomIndexPerType();
}

template<typename Iterator>
void NearestNeighborNSquare::computeNeighborsCartesianCoordinatesBatch(const Iterator strucBegin,
                                                                       const Iterator strucEnd,
                                                                       const AtomBatchMap& batchMap)
{
    nAtoms =
        std::accumulate(strucBegin,
                        strucEnd,
                        0,
                        [](Int sum, const Structure& struc) { return sum + struc.get_numAtoms(); });
    allocateStorageArrays(nAtoms);
    const Vec1Int&    typeIdx = batchMap.get_typesIdx();
    const Vec1String& type = batchMap.get_types();
    std::copy(typeIdx.cbegin(), typeIdx.cend(), typeIndexCentral.begin());
    std::copy(type.cbegin(), type.cend(), typeKeyCentral.begin());
    computeNumberAtomsPerType(batchMap.get_typeOrder());
    std::for_each(strucBegin,
                  strucEnd,
                  [&](const Structure& struc)
                  {
                      for (Size i = 0; i < struc.get_numAtoms(); i++)
                      {
                          computeNearestNeighborsSingleAtomCartesian(i, i, struc);
                      }
                  });
    compute_centralAtomIndexPerType();
}

using NeighborListMap = std::map<String, std::shared_ptr<NearestNeighborNSquare>>;

namespace neighbor_list
{
/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapNearestNeighborNSquare();
/*******************************************************************************************
 * @brief Copies the parameters from one NearestNeighborNSquare object to another.
 *
 * @param listA The source NearestNeighborNSquare object from which parameters are copied.
 * @param listB The target NearestNeighborNSquare object to which parameters are copied.
 *******************************************************************************************/
void copyParameters(const NearestNeighborNSquare& listA, NearestNeighborNSquare& listB);
/*******************************************************************************************
 * @brief Creates a new NearestNeighborNSquare object by copying parameters from an existing object.
 *
 * @param listA Pointer to the source NearestNeighborNSquare object from which parameters are
 *copied.
 * @return A pointer to the newly created NearestNeighborNSquare object with copied parameters.
 *******************************************************************************************/
NearestNeighborNSquare* copyParameters(const NearestNeighborNSquare& listA);
/*******************************************************************************************
 * @brief Copies the data from one NearestNeighborNSquare object to another.
 *
 * @param listA The source NearestNeighborNSquare object from which data is copied.
 * @param listB The target NearestNeighborNSquare object to which data is copied.
 *******************************************************************************************/
void copyData(const NearestNeighborNSquare& listA, NearestNeighborNSquare& listB);
/*******************************************************************************************
 * @brief Extends the data of one NearestNeighborNSquare object by appending data from another
 *object.
 *
 * @param listA The source NearestNeighborNSquare object from which data is appended.
 * @param listB The target NearestNeighborNSquare object to which data is extended.
 *******************************************************************************************/
void extendData(const NearestNeighborNSquare& listA, NearestNeighborNSquare& listB);
/*******************************************************************************************
 * @brief Adds data for a specific element from one NearestNeighborNSquare object to another.
 *
 * @param listA The source NearestNeighborNSquare object from which element data is copied.
 * @param listB The target NearestNeighborNSquare object to which element data is added.
 * @param atomIndx The index of the atom whose data is being added.
 *******************************************************************************************/
void addElementData(const NearestNeighborNSquare& listA,
                    NearestNeighborNSquare&       listB,
                    const Size                    atomIndx);
/*******************************************************************************************
 * @brief Retrieves the shared pointer to the maximum NearestNeighborNSquare object from a map.
 *
 * This function finds the NearestNeighborNSquare object with the maximum value in the given map
 * and returns a shared pointer to it. The comparison is performed using the `<` operator on the
 * objects pointed to by the shared pointers.
 *
 * @param neighborList A map where the key is a String and the value is a shared pointer to a
 *                     constant NearestNeighborNSquare object.
 * @return A shared pointer to the NearestNeighborNSquare object with the maximum value in the map.
 *         If the map is empty, the behavior is undefined.
 *******************************************************************************************/
std::shared_ptr<const NearestNeighborNSquare> getMaxNeighborListPtr(
    const std::map<String, std::shared_ptr<const NearestNeighborNSquare>>& neighborLists);

ShRec extractNeighborListFromBatchList(const NearestNeighborNSquare& neighborList,
                                       const AtomBatchMap&           batchMap,
                                       const Int                     nStruc);
/*******************************************************************************************
 * @brief Extracts and filters neighbor elements from a given neighbor list based on a specified set
 *of atoms.
 *
 * This function creates a new neighbor list by filtering the data from the input `neighborList`
 * according to the indices specified in the `atoms` vector. The filtered data includes various
 * properties such as global indices, type indices, distances, connection vectors, and neighbor
 *counts.
 *
 * @param neighborList The input neighbor list of type `NearestNeighborNSquare` containing the full
 *set of neighbors.
 * @param atoms A vector of atom indices (`Vec1Int`) used to filter the neighbor list.
 * @return A new `NearestNeighborNSquare` object containing the filtered neighbor data.
 *
 * The following fields are filtered and included in the resulting neighbor list:
 * - `globalIndex`: Global indices of neighbors.
 * - `typeIndex`: Type indices of neighbors.
 * - `typeIndexCentral`: Type indices of central atoms.
 * - `typeKeyCentral`: Type keys of central atoms.
 * - `distances`: Distances between atoms and their neighbors.
 * - `connectionVector`: Connection vectors between atoms and their neighbors.
 * - `connectionVectorNormalized`: Normalized connection vectors.
 * - `numberNeighbors`: Number of neighbors for each atom.
 * - `numberNeighborsType`: Number of neighbors of each type for each atom.
 *
 * The resulting neighbor list retains the cutoff, type sorting, and distance sorting properties
 * of the input `neighborList`.
 *******************************************************************************************/
NearestNeighborNSquare getNeighborElements(const NearestNeighborNSquare& neighborList,
                                           const Vec1Int&                atoms);
} // namespace neighbor_list
} //namespace vaspml

#endif
