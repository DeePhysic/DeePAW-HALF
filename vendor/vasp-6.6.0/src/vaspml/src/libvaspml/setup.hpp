#ifndef SETUP_HPP
#define SETUP_HPP

#include "BasisFunctions.hpp"
#include "Descriptor.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

#include <variant>

namespace vaspml
{

class MlMPI;
struct Record;
class Selector;
class Structure;
template <class T> class ShmemArray2DVariableLen;

enum class Descriptor3BodyType
{
    StandardDescriptor,
    LinearScalingDescriptor
};

} //namespace vaspml

namespace vaspml::setup
{

/*******************************************************************************************
 * @brief Constructs a map of basis functions for 2-body and 3-body interactions.
 *
 * This function initializes and returns a `BasisFunctionMap` containing
 * radial spline basis functions for 2-body interactions and angular basis
 * functions for 3-body interactions. The parameters for these basis functions
 * are extracted from the provided `Record` object.
 *
 * @param record A `Record` object containing the configuration parameters
 *               for the basis functions.
 *               - "ML_ICUT1": Cutoff type for 2-body interactions (integer).
 *               - "ML_RCUT1": Radial cutoff for 2-body interactions (real).
 *               - "ML_SION1": Smoothing parameter for 2-body interactions (real).
 *               - "ML_NR1": Number of radial points for 2-body interactions (integer).
 *               - "ML_MRB1": Maximum radial basis for 2-body interactions (integer).
 *               - "ML_ICUT2": Cutoff type for 3-body interactions (integer).
 *               - "ML_RCUT2": Radial cutoff for 3-body interactions (real).
 *               - "ML_SION2": Smoothing parameter for 3-body interactions (real).
 *               - "ML_NR2": Number of radial points for 3-body interactions (integer).
 *               - "ML_LMAX2": Maximum angular momentum for 3-body interactions (integer).
 *               - "ML_MRB2": Maximum radial basis for 3-body interactions (integer).
 *
 * @return A `BasisFunctionMap` containing:
 *         - "2-body": A shared pointer to a `BasisFunctionsRadialSpline` object.
 *         - "3-body": A shared pointer to a `BasisFunctionsAngular` object.
 *******************************************************************************************/
BasisFunctionMap makeBasisFunctions(const Record& record);
/*******************************************************************************************
 * @brief Creates and initializes a storage container for structures.
 *
 * This function generates a vector of `Structure` objects based on the input
 * data provided by the `mlab` object. It retrieves the number of structures
 * and their corresponding data, then constructs each `Structure` and stores
 * it in the resulting vector.
 *
 * @param mlab A shared reference to the input data container (`ShRec`) that
 *             provides the number of structures and their associated data.
 * @return A vector of `Structure` objects initialized with the input data.
 *******************************************************************************************/
std::vector<Structure> makeStructureStorage(const ShRec& ML_ABdata);
/***************************************************************************************************
 * Generate dictionary with all neighbor lists required by a given force field setup.
 *
 * @param cutoff Map with cutoffs for different neighbor lists.
 * @param data Frame memory (Record).
 *
 * @return Map with all required instances of neighbor lists.
 *
 * @pre Currently createFrameMemory() needs to be called in advance to create the desired multi-type
 * map entries.
 **************************************************************************************************/
NeighborListMap createNeighborLists(const std::map<String, Real>& cutoff, ShRec data);
/*******************************************************************************************
 * @brief Creates a map of 3-body descriptors based on the input parameters.
 *
 * This function generates descriptors for a given record using specified basis functions,
 * descriptor type, and execution policy. The descriptors are used for representing
 * 3-body interactions in a computationally efficient manner.
 *
 * @param incar              Input record containing the necessary data for descriptor generation.
 * @param basis3Body         Basis functions used for constructing the 3-body descriptors.
 * @param descriptor3BodyType Type of 3-body descriptor to be created (e.g., specific algorithm or
 *format).
 * @param data               Shared data or resources required for descriptor computation.
 *
 * @return DescriptorMap     A map containing the generated 3-body descriptors.
 *******************************************************************************************/
DescriptorMap createDescriptors(const Record&              incar,
                                const Size&                nTypes,
                                const BasisFunctions&      basis3Body,
                                const Descriptor3BodyType& descriptor3BodyType,
                                ShRec                      data);
/*******************************************************************************************
 * @brief Creates a map of descriptors for 2-body and 3-body interactions based on the provided
 *setup and parameters.
 *
 * This function initializes and configures descriptors for 2-body and 3-body interactions using the
 *provided input parameters, basis functions, and execution policy. The descriptors are stored in a
 *map, where the keys represent the descriptor names and the values are shared pointers to the
 *corresponding descriptor objects.
 *
 * @param setup                The input record containing configuration parameters.
 * @param basis3Body           Basis functions for 3-body interactions.
 * @param descriptor3BodyType  Type of 3-body descriptor to create (e.g., standard or linear
 *scaling).
 * @param data                 Shared record containing additional data required for descriptor
 *initialization.
 *
 * @return DescriptorMap       A map containing the initialized descriptors for 2-body and 3-body
 *interactions.
 *
 * @note The function supports two types of 3-body descriptors:
 *       - StandardDescriptor: A standard 3-body descriptor.
 *       - LinearScalingDescriptor: A linear scaling 3-body descriptor (currently commented out).
 *
 * @details
 * - The "SHS2-2-body" and "SHS2-3-body" descriptors are always created using the `DescriptorSHS2`
 *class.
 * - If the `descriptor3BodyType` is `StandardDescriptor`, an additional "SHS3-3-body" descriptor is
 *created using the `DescriptorSHS3` class.
 * - The `LinearScalingDescriptor` type is currently commented out and not implemented in this
 *function.
 *******************************************************************************************/
DescriptorMap createDescriptors(const Record&              incar,
                                const BasisFunctions&      basis3Body,
                                const Descriptor3BodyType& descriptor3BodyType,
                                ShRec                      data);
/*******************************************************************************************
 * @brief Creates a map of number descriptors for SHS2 (2-body) and SHS3 (3-body) interactions.
 *
 * This function generates a map where the keys are descriptor names ("SHS3-3-body" and
 *"SHS2-2-body"), and the values are shared pointers to vectors of integers representing the number
 *of descriptors per type for each interaction type.
 *
 * @param setup A reference to a Record object containing the configuration data.
 *              It is expected to provide the following keys:
 *              - "SHS3-3-body-number-descriptors-per-type": A vector of integers for 3-body
 *descriptors.
 *              - "SHS2-2-body-number-descriptors-per-type": A single integer for 2-body
 *descriptors.
 *
 * @return A map where:
 *         - The key "SHS3-3-body" maps to a shared pointer to a vector of integers for 3-body
 *descriptors.
 *         - The key "SHS2-2-body" maps to a shared pointer to a vector of integers for 2-body
 *descriptors, with each element initialized to the same value.
 *******************************************************************************************/
std::map<String, ShVec1Int> make_numberDescriptorsMap(const Record& setup);
/*******************************************************************************************
 * @brief Generates a map representing the feature space configuration for a given setup.
 *
 * This function creates a mapping between feature space keys (e.g., descriptor types)
 * and their corresponding sizes, which are computed based on the provided setup record.
 * The sizes are determined using scalar multiplication or element-wise product operations
 * with a reference vector (`numLrc`).
 *
 * @param setup The input record containing configuration data, including the reference vector
 *              (`numLrc`) and descriptor-related parameters.
 * @return A map where the keys are feature space identifiers (e.g., "SHS2-2-body") and the
 *         values are shared pointers to vectors of integers representing the computed sizes.
 *
 * @note The function currently handles specific descriptor keys ("SHS2-2-body" and "SHS3-3-body")
 *       explicitly. The commented-out section suggests a more generic implementation using a
 *       descriptor key list and conditional checks, which is yet to be implemented.
 *
 * @todo Refactor to use a more generic approach (e.g., `getif`) for handling descriptor keys
 *       dynamically, if possible.
 *******************************************************************************************/
std::map<String, ShVec1Int> generate_featureSpaceMap(const Record& setup);
/*******************************************************************************************
 * @brief Creates and initializes a shared memory 2D array for regression coefficients.
 *
 * This function constructs a shared memory 2D array with variable-length rows to store
 * regression coefficients. The array is initialized with values provided in the input
 * parameters and is distributed across MPI processes.
 *
 * @param inputParameters A record containing the input parameters. It must include:
 *                        - "regression-coeff" (ShVec2Real): The regression coefficients.
 *                        - "numLrc" (Vec1Int): The sizes of the rows in the 2D array.
 * @param mpiIn A shared pointer to an MlMPI object for MPI communication.
 * @return A shared pointer to the initialized ShmemArray2DVariableLen<Real> object.
 *******************************************************************************************/
std::shared_ptr<ShmemArray2DVariableLen<Real>> make_regressionCoefficients(
    const Record&                 inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn);
/*******************************************************************************************
 * @brief Creates a map of weights based on input parameters.
 *
 * This function generates a `std::map` where the keys are descriptive strings
 * representing different weight categories, and the values are the corresponding
 * weights retrieved from the provided input parameters.
 *
 * @param inputParameters A constant reference to a `constRecord` object containing
 *                        the input parameters from which the weights are extracted.
 * @return A `std::map<String, Real>` where:
 *         - "SHS2-2-body" maps to the value of "ML_W1" from `inputParameters`.
 *         - "SHS3-3-body" maps to the value of "ML_W2" from `inputParameters`.
 *******************************************************************************************/
std::map<String, Real> make_weightsMap(const Record& inputParameters);

/*******************************************************************************************
 * @brief Creates and initializes a Selector object.
 *
 * This function constructs a Selector object and initializes it using the provided input records
 * and MPI shared pointer.
 *
 * @param incar Input record for initialization.
 * @param mlab Label record (not used in this function, but part of the signature).
 * @param data Data record for initialization.
 * @param mlmpi Shared pointer to an MlMPI object for initialization.
 * @return A fully initialized Selector object.
 *******************************************************************************************/
Selector makeSelector(const ShRec&                  incar,
                      const ShRec&                  mlab,
                      const ShRec&                  data = nullptr,
                      const std::shared_ptr<MlMPI>& mlmpi = nullptr);
/*******************************************************************************************
 * @brief Generates a vector of sample counts for the farthest point selection process.
 *
 * This function retrieves the number of local reference configuration selections (`ML_SNCONF`)
 * from the input `incar` and ensures that it matches the number of atom types (`maxTypes`)
 * specified in the `mlab` object. If the sizes do not match, an error message is logged.
 *
 * @param incar The input configuration object containing the `ML_SNCONF` parameter.
 * @param mlab The machine learning configuration object containing the `maxTypes` parameter.
 * @return A vector of integers (`Vec1Int`) representing the number of samples for each atom type.
 *
 * @throws Logs an error if the size of `ML_SNCONF` does not match the number of atom types in
 *`mlab`.
 *******************************************************************************************/
Vec1Int makeFarthestPointNumberSamples(const ShRec& incar, const ShRec& mlab);
/*******************************************************************************************
 * @brief Creates a vector of farthest point thresholds based on input configurations.
 *
 * This function retrieves the farthest point thresholds from the input configuration (`incar`)
 * and validates that the number of thresholds matches the number of atom types in the machine
 * learning atomic basis (`mlab`). If the sizes do not match, an error is logged.
 *
 * @param incar Input configuration containing the farthest point thresholds.
 * @param mlab Machine learning atomic basis configuration.
 * @return A vector of farthest point thresholds.
 *
 * @throws Logs an error if the number of thresholds does not match the number of atom types.
 *******************************************************************************************/
Vec1Real makeFarthestPointTresholds(const ShRec& incar, const ShRec& mlab);
/*******************************************************************************************
 * @brief Creates a farthest point criterion based on the specified algorithm.
 *
 * This function determines the farthest point criterion to use based on the
 * "ML_SALGO" parameter in the input configuration. If the parameter is set to
 * "fpsn", it invokes `makeFarthestPointNumberSamples`. Otherwise, it invokes
 * `makeFarthestPointTresholds`.
 *
 * @param incar The input configuration record.
 * @param mlab The machine learning auxiliary data.
 * @return A `std::variant` containing either `Vec1Int` or `Vec1Real`,
 *         representing the farthest point criterion.
 *******************************************************************************************/
std::variant<Vec1Int, Vec1Real> makeFarhtestPointCriterion(const ShRec& incar, const ShRec& mlab);
/*******************************************************************************************
 * @brief Rescales specific INCAR parameters by a given distance factor.
 *
 * This function modifies the values of the INCAR parameters "ML_RCUT1",
 * "ML_RCUT2", "ML_SION1", and "ML_SION2" by multiplying them with the
 * provided distance factor.
 *
 * @param incar Reference to the INCAR record containing the parameters to be rescaled.
 * @param distFactor The factor by which the parameters are scaled.
 *******************************************************************************************/
void rescaleIncarUnits(Record& incar, const Real distFactor = 1.0);
/*******************************************************************************************
 * @brief Rescales various units in the given `mlab` record by applying the specified scaling
 *factors.
 *
 * This function adjusts the units of energy, distance, force, stress, and lattice parameters
 * in the `mlab` record and its associated structures. The scaling factors are applied to the
 * respective fields to ensure consistency with the desired unit system.
 *
 * @param mlab         The record containing the data to be rescaled.
 * @param energyFactor The scaling factor for energy-related quantities.
 * @param distFactor   The scaling factor for distance-related quantities.
 * @param forceFactor  The scaling factor for force-related quantities.
 * @param stressFactor The scaling factor for stress-related quantities.
 *
 * The following fields are rescaled:
 * - `atomicRefEnergy` is scaled by `energyFactor`.
 * - Each structure in the `structures` field:
 *   - `positions` is scaled by `distFactor`.
 *   - `energy` is scaled by `energyFactor`.
 *   - `forces` is scaled by `forceFactor`.
 *   - `stress` is scaled by `stressFactor`.
 *   - `lattice` is scaled by `distFactor`.
 *******************************************************************************************/
void rescaleMlabUnits(Record&    mlab,
                      const Real energyFactor = 1.0,
                      const Real distFactor = 1.0,
                      const Real forceFactor = 1.0,
                      const Real stressFactor = 1.0);

void rescaleMlffUnits(Record&    mlff,
                      const Real energyFactor = 1.0,
                      const Real distFactor = 1.0,
                      const Real forceFactor = 1.0,
                      const Real stressFactor = 1.0);

} //namespace vaspml::setup

#endif
