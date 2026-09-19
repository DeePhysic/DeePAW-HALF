#ifndef FRAME_HPP
#define FRAME_HPP

#include "BasisFunctions.hpp"
#include "BatchMap.hpp"
#include "Descriptor.hpp"
#include "Structure.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

namespace vaspml
{

class TypeMap;
enum class Descriptor3BodyType;
struct Record;

/***************************************************************************************************
 * This class contains one structure, the neighbor lists and its descriptors.
 **************************************************************************************************/
class Frame
{
  public:
    /***********************************************************************************************
     * Set up empty frame class instance.
     *
     * @warning Before use, the init() member function needs to be called.
     **********************************************************************************************/
    Frame(ShRec record = nullptr);
    /*******************************************************************************************
     * @brief Initializes the Frame object with the provided parameters and basis functions.
     *
     * This method sets up the Frame object using the given input parameters, a specific
     * 3-body basis function, and the execution policy. It also configures the basis functions
     * for SHS2 processing.
     *
     * @param inputParameters The input parameters for initializing the Frame.
     * @param basisFunctions  A map of basis functions, where the "3-body" key is used.
     *******************************************************************************************/
    void init(const Record& inputParameters, const BasisFunctionMap& basisFunctions);
    /*******************************************************************************************
     * Initialze default constructed Frame class for refitting.
     *
     * Calls init() and set_basisFunctionsSHS2
     *
     * @param incar Record instance conatining incar with force field parameters.
     * @param types vector containing the types of the structure as strings
     * @param BasisFunctionMap map containing the basis functions for the 2 and 3-body descriptors
     *******************************************************************************************/
    void initRefit(const Record&           incar,
                   const Vec1String&       types,
                   const BasisFunctionMap& basisFunctions);
    /*******************************************************************************************
     * @brief Updates neighbor lists and structure descriptors of the Frame object based on the
     *provided structure.
     *
     * This function performs two main tasks:
     * 1. Updates the neighbor lists of the Frame using the given `Structure` object.
     * 2. Updates the 2body and 3-body descriptors.
     *
     * @param structure Reference to a `Structure` object used to update the Frame's state.
     *******************************************************************************************/
    void update(Structure& structure, const Vec1Int& distributed = Vec1Int());
    /***********************************************************************************************
     * Compute all descriptors, neighbor list is assumed to be present and correctly filled.
     *
     * @note BasisFunctions for Descriptors has to be set manually. And neighbor lists
     * have to be filled previously.
     **********************************************************************************************/
    void update();
    /*******************************************************************************************
     * updating reading structure from POSCAR neighbor lists for the SHS2 and SHS3 descriptors.
     *
     * @param fileName name of file from which the strcuture is taken.
     * @param unitFactor unit conversion factor for length dimensions.
     *
     * @note the file format has to be a VASP POSCAR file.
     *******************************************************************************************/
    void updateNeighborLists(const String& fileName, const Real unitFactor = 1.0);
    /*******************************************************************************************
     * @brief Updates the neighbor lists for the current frame based on the given structure.
     *
     * This function performs the following steps:
     * - Converts the structure's coordinates from Cartesian to direct.
     * - Retrieves the structure's type names and initializes a type mapping.
     * - Updates all neighbor lists by computing nearest neighbors in direct coordinates.
     * - Associates updated neighbor lists with corresponding descriptors.
     *
     * @param structure Reference to the Structure object containing atomic data.
     *******************************************************************************************/
    void updateNeighborLists(Structure& structure, const Vec1Int& distributed = Vec1Int());
    /*******************************************************************************************
     * updating neighbor lists with supplied Record containing positions and lattice.
     *
     * @param data contains positions and lattice
     *******************************************************************************************/
    void updateNeighborLists(const ShRec& data);
    /*******************************************************************************************
     * assign used basis functions for the DescriptorSHS2 class
     *
     * @param basisFunctions dictionary containing the basis functions to be used
     *******************************************************************************************/
    void set_basisFunctionsSHS2(const BasisFunctionMap& basisFunctions);
    /*******************************************************************************************
     * update the array sizes of the working arrays for the DescriptorSHS2 classes
     *
     * @note before calling update_arraySizesSHS2 the basis functions have to be set by
     * calling set_basisFunctionsSHS2 and the neighbor lists have to be set by updateNeighborLists
     *******************************************************************************************/
    void update_arraySizesSHS2();
    /*******************************************************************************************
     * @brief Prepares the SHS3 descriptor by setting its parameters and generating lookup tables if
     *required.
     *
     * This function initializes the SHS3 descriptor by retrieving the SHS2 descriptor, casting it
     *appropriately, and setting the necessary parameters. If the execution policy is set to GPU
     *standard library, it also generates lookup tables for the SHS3 descriptor.
     *******************************************************************************************/
    void prepare_DescriptorSHS3();
    /*******************************************************************************************
     * Read structure from file in POSCAR format and compute neighbor lists and descriptors.
     *
     * @param fileName Name of file containing structure in POSCAR format.
     * @param factor unit conversion factor for length dimensions.
     *******************************************************************************************/
    void update(const String& fileName, const Real unitFactor = 1.0);
    /*******************************************************************************************
     * update descriptors from previosuly filled neighbor list array
     *
     * @note SHS2-2-body and SHS2-3-body
     * use BasisFunctionsRadialSpline. SHS3-3-body uses BasisFunctionsAngular.
     * All three kinds are stored in the map basisFunctions.
     *******************************************************************************************/
    void update(const std::shared_ptr<TypeMap>& typeMap);
    //**********************************************************************************************
    // GETTERS
    //**********************************************************************************************
    /*******************************************************************************************
     * Get descriptors as constant reference.
     *
     * Call with const DescriptorMap& descriptor = get_descriptors() to avoid copy constructors
     *******************************************************************************************/
    const DescriptorMap& get_descriptors() const;
    /*******************************************************************************************
     * Getting the type map via shared ptr. Only shared pointer creator will be called.
     *******************************************************************************************/
    std::shared_ptr<const TypeMap> get_typeMap() const;
    /*******************************************************************************************
     * @brief Retrieves the neighbor lists associated with the frame.
     *
     * This function provides access to the map of neighbor lists, which
     * represents the relationships or connections between elements in the frame.
     *
     * @return A constant reference to the NeighborListMap containing the neighbor lists.
     *******************************************************************************************/
    const NeighborListMap& get_neighborLists() const;
    /*******************************************************************************************
     * @brief Retrieves the associated structure of the frame.
     *
     * @return A shared pointer to the associated Structure object.
     *******************************************************************************************/
    std::shared_ptr<Structure> get_structure() const;

    /***********************************************************************************************
     * Record collecting all memory required for this Frame and its member instances of
     * NeighborList and Descriptor.
     *
     * This will be set up in the init() function by createFrameMemory().
     **********************************************************************************************/
    ShRec data;
    /*******************************************************************************************
     * A single structure for which neighbor lists and descriptors are calculated in this Frame.
     *******************************************************************************************/
    std::shared_ptr<Structure> structure;
    /***********************************************************************************************
     * Map of neighbor lists corresponding to different descriptors, cutoffs, etc.
     *
     * This will be set up in the init() function by createNeighborLists().
     **********************************************************************************************/
    NeighborListMap neighborLists;
    /***********************************************************************************************
     * Map of descriptors required for force field evaluation.
     *
     * This will be set up in the init() function by createDescriptors().
     **********************************************************************************************/
    DescriptorMap descriptors;
    /*******************************************************************************************
     * Map from force field types to types in structure.
     *******************************************************************************************/
    std::shared_ptr<TypeMap> typeMap;
    /*******************************************************************************************
     * List of types present in force field.
     *******************************************************************************************/
    ShCVec1String forceFieldTypes;
    /*******************************************************************************************
     * List of types present in structure.
     *******************************************************************************************/
    Vec1String structureTypes;
    /*******************************************************************************************
     * Variable selecting the descriptor type for the 3-body descriptor
     *******************************************************************************************/
    Descriptor3BodyType descriptor3BodyType;
    AtomBatchMap        batchMap;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapFrame();

} // namespace vaspml

#endif
