#ifndef DESCRIPTORCOLLECTOR_HPP
#define DESCRIPTORCOLLECTOR_HPP

#include "ArrayResizing.hpp"
#include "types.hpp"

namespace vaspml
{

class Descriptor;
class NearestNeighborNSquare;
class TypeMap;

enum class DescriptorStorage
{
    Type,
    CentralAtom,
};

class DescriptorCollector
{
  public:
    /**
     *@param storageOrder determines the storage order of the arrays.
     *Setting 0 makes type order and setting it 1
     *makes central atom ordering
     *@param normalizedDescriptors here a map array can be supplied
     *if the memory should not be managed by the class
     *itself
     @param normsAtom supply a storage space where the descriptor norms are stored per central atom
     */
    DescriptorCollector(
        Int                                 storageOrder = 0,
        const std::map<String, ShVec2Real>& normalizedDescriptors = std::map<String, ShVec2Real>(),
        ShRec                               record = nullptr);
    /**
     * updating the descriptor collector class
     *
     *@param descriptors map which contains the descriptors that should be collected and normalized.
     *The objects in the map have to be inherited from the Descriptor class
     */
    void updateCollector();
    /**
     *get the normalized descriptors as array by const reference
     */
    const std::map<String, ShVec2Real>& get_descriptorsNormalized() const;
    /*******************************************************************************************
     * @brief Retrieves the collection of descriptors.
     *
     * This method provides access to the internal map of descriptors, where each descriptor
     * is identified by a unique string key and stored as a shared pointer.
     *
     * @return A constant reference to the map of descriptors.
     *******************************************************************************************/
    const std::map<String, std::shared_ptr<Descriptor>>& get_descriptors() const;
    /*******************************************************************************************
     * get norm over descriptor list for single atom
     *
     * @param atomIndx index of central atom
     *******************************************************************************************/
    Real get_normAtom(Size atomIndx) const;
    /*******************************************************************************************
     * get normalization factors of central atoms
     *
     * normalization factors are computed over all descriptors of the specific atom
     *
     * @note call with const Vec1Real& x = get_normAtom(); to avoid making copy
     *******************************************************************************************/
    const Vec1Real& get_normAtom() const;
    /*******************************************************************************************
     * rearrange normalized descriptors in case that ML_FF file contains more types than
     * the actual structure.
     *
     * this routine selects CPU or GPU version
     *******************************************************************************************/
    void rearrangeSHS2Body(const TypeMap& typeMap);
    /*******************************************************************************************
     * rearrange normalized descriptors - cpu version using std::lib
     *
     * no parallel execution policy is used in routine
     *******************************************************************************************/
    void rearrangeSHS2BodyTypeOrderCPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * @brief Rearranges the SHS2 2-body descriptors for atoms on the CPU based on the atom type
     *mapping.
     *
     * This function modifies the normalized SHS2 2-body descriptors by rearranging their order
     * according to the force field types and structure types defined in the provided type map.
     *
     * @param typeMap A mapping between structure types and force field types used for rearranging
     *the descriptors.
     *
     * The function performs the following steps:
     * - Retrieves the normalized SHS2 2-body descriptors and the number of descriptors per type.
     * - Iterates over each descriptor slice and resizes it to accommodate the rearranged data.
     * - Rearranges the descriptor values based on the mapping between structure types and force
     *field types.
     *******************************************************************************************/
    void rearrangeSHS2BodyAtomOrderCPU(const TypeMap& typeMap);

    /*******************************************************************************************
     * rearrange normalized descriptors - gpu version using std::lib
     *******************************************************************************************/
    void rearrangeSHS2BodyTypeOrderGPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * returns one of the collected Descriptor in non normalized form.
     *
     * Can for example be used to retrieve neighbor lists
     *******************************************************************************************/
    const Descriptor& getDescriptor(const String& key) const;
    /*******************************************************************************************
     * Set local shared pointer to the  input descriptor map
     *
     * @param descriptors descriptor map to which the local shared escriptor pointer will be set
     *******************************************************************************************/
    void setDescriptorMap(const std::map<String, std::shared_ptr<Descriptor>>& descriptors);
    /*******************************************************************************************
     * Write the normalized descriptors to a files in single column format
     *******************************************************************************************/
    void writeDescriptorCollector() const;
    /*******************************************************************************************
     * @brief Retrieves the nearest neighbor list associated with a given key.
     *
     * This function accesses the descriptor map using the provided key and
     * returns the nearest neighbor list associated with the corresponding descriptor.
     *
     * @param key The key used to identify the desired descriptor in the map.
     * @return A constant reference to the NearestNeighborNSquare object representing the neighbor
     *list.
     * @throws std::out_of_range If the key is not found in the descriptor map.
     *******************************************************************************************/
    const NearestNeighborNSquare& get_neighborList(const String& key) const;
    /*******************************************************************************************
     * @brief Retrieves the neighbor list with the maximum weight from the descriptors.
     *
     * This function iterates through a predefined list of descriptor keys and collects
     * the neighbor lists associated with descriptors that have a positive weight. It then
     * determines and returns the neighbor list with the maximum weight.
     *
     * @return A shared pointer to the neighbor list with the maximum weight.
     *******************************************************************************************/
    std::shared_ptr<const NearestNeighborNSquare> getMaxNeighborList() const;

    const Real* get_descriptorNormalizedPtrTO(const Int     type,
                                              const Int     atom,
                                              const String& key) const;
    void        computeBracketTerms(const Int                                            kAtom,
                                    const std::map<String, std::shared_ptr<Descriptor>>& descriptors,
                                    const NearestNeighborNSquare&                        maxNeighborList);
    void        compute_innerBracketSlowDerivative(const Int                     kAtom,
                                                   const NearestNeighborNSquare& maxNeighborList);
    void        compute_innerBracketSlowDerivativeHead(const Int                     kAtom,
                                                       const NearestNeighborNSquare& maxNeighborList);
    void        compute_outerBracketSlowDerivative();

  private:
    void allocate_innerBracketSlowDerivative(const Int numberNeighbors);
    void allocate_innerBracketSlowDerivativeHead(const Int numberTypes);
    void allocate_outerBracketSlowDerivative();
    /**
     *Normalizing the supplied descriptors and storing them in the chosen format.
     *Choses whether CPU or GPU version is used.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptors();
    /**
     *Normalizing the supplied descriptors and storing them in the chosen format
     *
     *CPU implementation using for loops and standard library algorithms.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptorsCPU();
    /**
     *normalizing the supplied descriptors and storing them in the chosen format
     *GPU implementation using standard library algorithms.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptorsGPU();
    /*******************************************************************************************
     * allocate arrays descriptorsNormalized in type ordered or central atom format
     *******************************************************************************************/
    void allocateArrays();
    /*******************************************************************************************
     * Allocate arrays for normalized descriptors for the CPU. Resize is done by a std vector resize
     *******************************************************************************************/
    void allocateArraysCPU();
    /*******************************************************************************************
     * Allocate arrays for normalized descriptors for the GPU. Resize is done to avoid
     * device host copies because the is currently no resize on the device possible.
     *******************************************************************************************/
    void allocateArraysGPU();
    /*******************************************************************************************
     * compute norm of a single atom over the set of all descriptors
     *
     * @param atomIndex index of central atom for which the descriptor norm is computed
     *******************************************************************************************/
    Real computeNormSingleAtom(Size atomIndex);
    /*******************************************************************************************
     * Compute normalized descriptors on the GPU
     *
     * @param atomIndex index of central atom to choose descriptor array to normalize.
     *******************************************************************************************/
    Real computeNormSingleAtomGPU(const Size& atomIndex);
    /**
     *rescale descriptors by the computed norms and multiply times the
     *descriptor weight store in type order format
     *
     *@param atomIndex index of central atom
     *@param norm computed norm of descriptor for certain central atom and summed over descriptor
     *types
     *@param descriptors map of classes to obtain the non-normalized descriptor elements
     */
    void normalizeDescriptorsSingleAtomTypeOrder(const Size& atomIndex, const Real norm);
    /**
     *rescale descriptors by the computed norms and multiply times the
     *descriptor weight; store in Central atom format
     *
     *@param atomIndex index of central atom
     *@param norm computed norm of descriptor for certain central atom and summed over descriptor
     *types
     *@param descriptors map of classes to obtain the non-normalized descriptor elements
     */
    void normalizeDescriptorsSingleAtomCentralOrder(Size atomIndex, const Real norm);
    /*******************************************************************************************
     * Resize descriptorsNormalized and fill with 0.
     *
     * This function will only be called if the number of force field types does not match
     * the number of types in the structure.
     *******************************************************************************************/
    void allocateDescriptorSHS2BodyNormalized(const TypeMap& typeMap);

    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /**
     * vector of the descriptors needed by the descriptor collector class
     */
    std::map<String, std::shared_ptr<Descriptor>> descriptors;
    /**
     * std::map to store the normalized for the supplied descriptors.
     * @NOTE
     * the array can either be stored as [type][ atoms of type x number descriptors ] or
     * [ central atom ][ number descriptors ]
     */
    std::map<String, ShVec2Real> descriptorsNormalized;
    /*******************************************************************************************
     * Object to controll if the arrays for the normalized descriptors have to be resized.
     *******************************************************************************************/
    std::map<String, ArrayResizing2D> descriptorsNormalizedSize;
    /**
     *smart enum to determine the storage order
     *
     *supported order types are type ordered. Type index determines first index
     *of storage arrays; Second index is a combination index of atoms x number of descriptors.
     *Central atom index ordering is second possibility. Then the first index is the central atom
     *index and the second index of the storage arrays is the number of descriptors
     */
    DescriptorStorage storageOrder;
    /**
     * storing the length of the array descriptorsNormalized
     *@NOTE
     * can be stored in a type or central atom fashion
     */
    std::map<String, ShVec1Int> length_descriptorsNormalized;
    /*******************************************************************************************
     * Object to controll if the lengths of the normalized descriptors have to be resized.
     *
     * Needed for the GPU version.
     *******************************************************************************************/
    std::map<String, ArrayResizing1D> length_descriptorsNormalizedSize;
    /**
     * norm computed over descriptor list
     *
     *@Note order is always over central atom. This order will match those of the descriptors
     *so when doing
     *loop type:
     *    loop atom per type:
     *the order matches normsAtom when incrementing by one
     */
    Vec1Real&       normsAtom;
    Vec1Size        centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
    // std::ofstream fileSHS2;
    // std::ofstream fileSHS3;
    /***********************************************************************
     * Number of atoms in the system
     ************************************************************************/
    Int        numberAtoms;
    Vec1String descriptorKeyList;
    // fit variables
    Vec1Real                   innerBracketSlowDerivative;
    Vec2Real                   innerBracketSlowDerivativeHead;
    std::map<String, Vec2Real> outerBracketSlowDerivative;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapDescriptorCollector();

} // namespace vaspml

#endif
