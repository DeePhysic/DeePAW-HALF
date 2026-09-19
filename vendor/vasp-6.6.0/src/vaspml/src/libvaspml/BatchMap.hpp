#ifndef BATCHMAP_HPP
#define BATCHMAP_HPP

#include "types.hpp"

#include <algorithm>      // for for_each
#include <iterator>       // for distance
#include <unordered_map>  // for unordered_map

namespace vaspml
{

class Structure;

class AtomBatchMap
{
  public:
    AtomBatchMap() = default;
    /*******************************************************************************************
     * Create all necessary maps for needed for batching.
     *
     * @param structures vector of set up Structure types
     * @param typeOrder optional vector to supply a desired type order according to which is sorted
     *
     * Creates a nstruc,atom -> batch index map
     * Creates a nstruc,atom -> batch index map in type ordered form
     * Creates a batch indx -> structure map
     * Creates a batch indx type ordered -> structure map
     * Creates a batch indx -> atom map
     * Creates a batch indx type ordered -> atom map
     *******************************************************************************************/
    void makeMap(const std::vector<Structure>& structures,
                 const Vec1String&             typeOrder = Vec1String{});
    /*******************************************************************************************
     * Create all necessary maps for needed for batching.
     *
     * @param strucStart start iterator of vector of set up Structure types
     * @param strucEnd end iterator of vector of set up Structure types
     * @param typeOrder optional vector to supply a desired type order according to which is sorted
     *
     * Creates a nstruc,atom -> batch index map
     * Creates a nstruc,atom -> batch index map in type ordered form
     * Creates a batch indx -> structure map
     * Creates a batch indx type ordered -> structure map
     * Creates a batch indx -> atom map
     * Creates a batch indx type ordered -> atom map
     *******************************************************************************************/
    template<typename Iterator>
    void makeMap(const Iterator    strucStart,
                 const Iterator    strucEnd,
                 const Vec1String& typeOrder = Vec1String{});
    /*******************************************************************************************
     * Get the index of the atom in the batch before type ordering.
     *
     * @param strucIdx supply structure indx
     * @param atomIdx supply atom indx within the structure
     *******************************************************************************************/
    Int get_mapStrucAtom_Batch(const Int& strucIdx, const Int& atomIdx) const;
    /*******************************************************************************************
     * Get the index of the atom back in the batch after type ordering was done
     *
     * Completely broken and luckily not used
     *
     * @param strucIdx supply structure indx as in input too batch
     * @param atomIdx supply atom indx within the structure
     *******************************************************************************************/
    Int get_mapStrucAtom_BatchTO(const Int& strucIdx, const Int& atomIdx) const;
    /*******************************************************************************************
     * Get the structure index of some atom in the batch before type ordering.
     *
     * @param flatBatchIdx index of atom of interest in super structure batch which is just
     * the atoms of all structures in batchh concatenated.
     *******************************************************************************************/
    Int get_strucIndex(const Int& flatBatchIdx) const;
    /*******************************************************************************************
     * Get the structure index of some atom in the batch after type ordering
     *
     * @param flatBatchIdx index of atom of interest in super structure batch.
     *******************************************************************************************/
    Int get_strucIndexTO(const Int& flatBatchIdx) const;
    /*******************************************************************************************
     * Get the original atom index within the input structure back before type ordering
     *
     * @param flatBatchIdx index of atom of interest in super structure batch.
     *******************************************************************************************/
    Int get_atomIndex(const Int& flatBatchIdx) const;
    /*******************************************************************************************
     * Get the original atom index within the input structure back before type ordering
     *******************************************************************************************/
    const Vec1Int& get_atomIndex() const;
    /*******************************************************************************************
     * Get the original atom index within the input structure back after type ordering
     *
     * @param flatBatchIdx index of atom of interest in super structure batch.
     *******************************************************************************************/
    Int get_atomIndexTO(const Int& flatBatchIdx) const;
    /*******************************************************************************************
     * Get the original atom index within the input structure back after type ordering
     *******************************************************************************************/
    const Vec1Int& get_atomIndexTO() const;
    /*******************************************************************************************
     * Get the original atom flattened atom index back by supplying type ordered index
     *
     * @param flatBatchIdx index of atom of interest in super structure batch.
     *******************************************************************************************/
    Int get_mapTO_OrigOrder(const Int& flatBatchIdx) const;
    Int get_mapOrigOrder_TO(const Int& flatBatchIdx) const;
    /*******************************************************************************************
     * Supply a flattened batch index after type ordering was done and get original index back
     *
     * @param inputIdx index of atom of interest in super structure batch.
     *******************************************************************************************/
    Int get_mapTO(const Int& inputIdx) const;
    /*******************************************************************************************
     * Get type index of a atom by suppling flattened batch index before type reordering
     *
     * @param atomIndx flattened batch index.
     *******************************************************************************************/
    String get_types(const Int& atomIndx) const;
    /*******************************************************************************************
     * Get type index of a atom before type reordering
     *******************************************************************************************/
    const Vec1String& get_types() const;
    /*******************************************************************************************
     * Get type index of a atom by suppling flattened batch index after type reordering
     *
     * @param atomIndx flattened batch index.
     *******************************************************************************************/
    Int get_typesIdxTO(const Int& atomIndx) const;
    /*******************************************************************************************
     * Get type index of a atom after type reordering
     *******************************************************************************************/
    String get_typesTO(const Int& atomIndx) const;
    /*******************************************************************************************
     * Get type indexes of all atoms in batch before type reordering
     *
     * @param atomIndx flattened batch index.
     *******************************************************************************************/
    Int get_typesIdx(const Int& atomIndx) const;
    /*******************************************************************************************
     * Get type indexes of all atoms in batch before type reordering
     *******************************************************************************************/
    const Vec1Int& get_typesIdx() const;
    /*******************************************************************************************
     * Get type indexes of all atoms in batch after type reordering
     *
     * @param atomIndx flattened batch index.
     *******************************************************************************************/
    const Vec1Int& get_typesIdxTO() const;
    /*******************************************************************************************
     * Get atom types after types reordering was done.
     *******************************************************************************************/
    const Vec1String& get_typesTO() const;
    /*******************************************************************************************
     * Get the unique types in the supplied structures.
     *
     * Defines the type order in which atoms are sorted after type sorting was done.
     *******************************************************************************************/
    const Vec1String& get_typeOrder() const;

    Int            get_numberAtomsPerStruc(const Int nStruc) const;
    const Vec1Int& get_numberAtomsPerStruc() const;

  private:
    /*******************************************************************************************
     * @brief Sets the number of atoms per structure in the batch.
     *
     * This function resizes the `numberAtomsPerStruc` vector to match the number
     * of structures in the input and populates it with the number of atoms in
     * each structure.
     *
     * @param structures A vector of `Structure` objects, each representing a structure
     *                   whose atom count will be stored.
     *******************************************************************************************/
    void set_numberAtomsPerStruc(const std::vector<Structure>& structures);
    /*******************************************************************************************
     * Creating flattening maps before type ordering was done.
     *
     * @param structures vector of structure classes which habve to be batched.
     *
     * Creates a nstruc,atom -> batch index map
     * Creates a batch indx -> structure map
     * Creates a batch indx -> atom map
     *******************************************************************************************/
    void makeFlatMaps(const std::vector<Structure>& structures);
    /*******************************************************************************************
     * Creating flattening maps before type ordering was done.
     *
     * @param strucBegin start iteratoor of vector of Structure classes which have to be batched.
     * @param strucBegin end iteratoor of vector of Structure classes which have to be batched.
     *
     * Creates a nstruc,atom -> batch index map
     * Creates a batch indx -> structure map
     * Creates a batch indx -> atom map
     *******************************************************************************************/
    template<typename Iterator>
    void makeFlatMaps(const Iterator strucBegin, const Iterator strucEnd);
    /*******************************************************************************************
     * Creating flattening maps before type ordering was done.
     *
     * Creates a nstruc,atom -> batch index map in type ordered form
     * Creates a batch indx type ordered -> structure map
     * Creates a batch indx type ordered -> atom map
     *******************************************************************************************/
    void makeTypeOrderMaps();
    /*******************************************************************************************
     * Find unique types in combined types array.
     *******************************************************************************************/
    void makeUniqueTypeOrder(const Vec1String& typeOrder);
    /*******************************************************************************************
     * Assign every atom type coded as a string a unique identifier coded as int
     *******************************************************************************************/
    void makeIntTypes();
    /*******************************************************************************************
     * Mapping from structure index and atom index to entry in flattened batch
     *******************************************************************************************/
    Vec2Int mapStrucAtom_Batch;
    /*******************************************************************************************
     * Mapping from structure index and atom index to entry in flattened batch after type reordering
     *******************************************************************************************/
    Vec2Int mapStrucAtom_BatchTO;
    /*******************************************************************************************
     * Map between flattaened arrays and flattened arrays after type reordering
     *******************************************************************************************/
    Vec1Int mapTO;
    /*******************************************************************************************
     * Storage of type strings in a flattened form. First atoms within structures change then strucs
     *******************************************************************************************/
    Vec1String types;
    /*******************************************************************************************
     * Storage of type strings in a flattened form. First atoms within structures change then strucs
     *
     * This is the order of types after atoms are sorted according to their types.
     *******************************************************************************************/
    Vec1String typesTO;
    /*******************************************************************************************
     * Maps from flattened structure atom index back to the structure entry of the atom before flat
     *******************************************************************************************/
    Vec1Int strucIndex;
    /*******************************************************************************************
     * Maps from flattened structure atom index back to the structure entry of the atom before flat
     *
     * Contains the order after type reordering was done.
     *******************************************************************************************/
    Vec1Int strucIndexTO;
    /*******************************************************************************************
     * Maps from flattened structure atom index back to the atom entry within structure before flat
     *******************************************************************************************/
    Vec1Int atomIndex;
    /*******************************************************************************************
     * Maps from flattened structure atom index back to the atom entry within structure before flat
     * in type order format
     *******************************************************************************************/
    Vec1Int atomIndexTO;
    /*******************************************************************************************
     * Maps from flattened type ordered index to original index.
     *
     * If a index in type order is supplied to the array it will return the atom index before
     * type ordering.
     *******************************************************************************************/
    Vec1Int mapTO_OrigOrder;
    /*******************************************************************************************
     * Maps from flattened original index to type order index.
     *******************************************************************************************/
    Vec1Int mapOrigOrder_TO;
    /*******************************************************************************************
     * Unique elements in combined element lists types of combined structures.
     *******************************************************************************************/
    Vec1String typeOrder;
    /*******************************************************************************************
     * Number of structures
     *******************************************************************************************/
    Size nStruc;
    /*******************************************************************************************
     * Stores the number of atom per structure.
     *******************************************************************************************/
    Vec1Int nAtomsInStruc;
    /*******************************************************************************************
     * Int descriptions for types converted to integer indexes starting at 0
     *******************************************************************************************/
    Vec1Int typesIdx;
    /*******************************************************************************************
     * Int descriptions for types converted to integer indexes starting at 0 after type reorder
     *******************************************************************************************/
    Vec1Int typesIdxTO;
    /*******************************************************************************************
     * Map converting a string to an integer index before type sorting.
     *
     * The mapping string -> int is different before and after type sorting. First string
     * occuring in Vec1String will get index 0.
     *******************************************************************************************/
    std::unordered_map<String, Int> strTypeToIntType;
    /*******************************************************************************************
     * Number of atoms in a certain structure. Order of vector is same as input structures
     *******************************************************************************************/
    Vec1Int numberAtomsPerStruc;
};

template<typename Iterator>
void AtomBatchMap::makeMap(const Iterator    strucBegin,
                           const Iterator    strucEnd,
                           const Vec1String& typeOrder)
{
    nStruc = std::distance(strucBegin, strucEnd);
    makeFlatMaps(strucBegin, strucEnd);
    makeUniqueTypeOrder(typeOrder);
    makeIntTypes();
    makeTypeOrderMaps();
}

template<typename Iterator>
void AtomBatchMap::makeFlatMaps(const Iterator strucBegin, const Iterator strucEnd)
{
    Size totalAtoms = 0;
    std::for_each(strucBegin,
                  strucEnd,
                  [&totalAtoms](const auto& struc) { totalAtoms += struc.get_numAtoms(); });

    types.resize(totalAtoms);
    strucIndex.resize(totalAtoms);
    atomIndex.resize(totalAtoms);
    atomIndexTO.resize(totalAtoms);
    mapStrucAtom_Batch.resize(nStruc);
    mapStrucAtom_BatchTO.resize(nStruc);
    nAtomsInStruc.resize(nStruc);
    Size counter = 0;
    auto struc = strucBegin;
    for (Size strucIdx = 0; strucIdx < nStruc; strucIdx++)
    {
        const auto& keys = struc->get_atomTypeNames();
        mapStrucAtom_Batch[strucIdx].resize(keys.size());
        mapStrucAtom_BatchTO[strucIdx].resize(keys.size());
        nAtomsInStruc[strucIdx] = keys.size();
        for (Size atom = 0; atom < keys.size(); atom++)
        {
            types[counter] = keys[atom];
            strucIndex[counter] = strucIdx;
            mapStrucAtom_Batch[strucIdx][atom] = counter;
            atomIndex[counter] = atom;
            counter++;
        }
        struc++;
    }
}

/*******************************************************************************************
 * @brief Manages the distribution of structures into batches of a fixed size.
 *
 * This class provides methods to determine the batch index for a given
 * structure, retrieve the indices of structures within a specific batch,
 * calculate the total number of batches, and find the index of a structure
 * within its batch.
 *******************************************************************************************/
class BatchManager
{

  public:
    /*******************************************************************************************
     * @brief Constructs a BatchManager with a specified batch size.
     *
     * @param batchSize The desired size of each batch. Must be greater than 0.
     * @throws bug message If the provided batchSize is 0.
     *******************************************************************************************/
    BatchManager(const Size& batchSize);
    /*******************************************************************************************
     * @brief Gets the index of the batch that a given structure belongs to.
     *
     * @param structureIndex The index of the structure.
     * @return The index of the batch containing the structure.
     * @throws std::out_of_range If the structureIndex is out of bounds.
     *******************************************************************************************/
    Size getBatchIndexForStructure(const Size& structureIndex) const;
    /*******************************************************************************************
     * @brief Gets the indices of all structures within a specific batch.
     *
     * @param batchIndex The index of the batch.
     * @return A vector containing the indices of the structures in the batch.
     *******************************************************************************************/
    std::vector<Size> getStructureIndicesInBatch(const Size& batchIndex) const;
    /*******************************************************************************************
     * @brief Gets the index of a structure within its assigned batch.
     *
     * @param structureIndex The index of the structure.
     * @return The index of the structure within its batch (0-based).
     * @throws bug message If the structureIndex is out of bound
     *******************************************************************************************/
    Size getStructureIndexWithinBatch(const Size& structureIndex) const;
    /*******************************************************************************************
     * @brief Gets the total number of batches required to hold all structures.
     *
     * @return The total number of batches.
     *******************************************************************************************/
    Size getNumberOfBatches() const;
    /*******************************************************************************************
     * @brief Sets the total number of structures managed by the BatchManager.
     *
     * This is necessary for bounds checking in other methods.
     *
     * @param numStructures The total number of structures.
     *******************************************************************************************/
    void setNumberOfStructures(const Size& numStructures);

  private:
    /*******************************************************************************************
     * @brief Sets the number of atoms per structure in the batch.
     *
     * This function resizes the `numberAtomsPerStruc` vector to match the number
     * of structures in the input and populates it with the number of atoms in
     * each structure.
     *
     * @param structures A vector of `Structure` objects, each representing a structure
     *                   whose atom count will be stored.
     *******************************************************************************************/
    void set_numberAtomsPerStruc(const std::vector<Structure>& structures);
    /*******************************************************************************************
     * The size of the batch.
     *******************************************************************************************/
    Size batchSize;
    /*******************************************************************************************
     * The total number of structures.
     *******************************************************************************************/
    Size numStructures = 0;
};

} //namespace vaspml

#endif
