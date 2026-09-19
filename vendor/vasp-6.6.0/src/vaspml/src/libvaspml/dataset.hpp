#ifndef DATASET_HPP
#define DATASET_HPP

#include "types.hpp"

#include <tuple>  // for tuple

namespace vaspml
{
class MlMPI;
class Structure;
class TypeMap;
struct Record;
}

namespace vaspml::dataset
{

void getPerStructureLrc(const Record& set, Vec1Int& lrcStructure, Vec1Int& lrcAtom);
void setPerStructureLrc(const Vec2Int& lrcStructure, const Vec2Int& lrcAtom, Record& set);
void concat(Record& set1, const Record& set2);

/*******************************************************************************************
 * Make a type map which allows the reorder structures to the form written in the
 * "The atom types in the data file" in the ML_AB file which was read.
 *
 * @param subTypes are the types of the current structures
 * @param types are the types which have to contain all types of subTypes and maybe more.
 *******************************************************************************************/
Vec1Int makeTypeMap(const Vec1String& subTypes, const Vec1String& types);
/*******************************************************************************************
 * Reordering the position, forces, types and number of types to the from
 * suggested by atomMap and typeMap.
 *
 * The typeMap is created by makeTypeMap and the atomMap is created by createReorderMap
 *******************************************************************************************/
void reorderSingleStructure(Record& strucSet, const Vec1Int& atomMap, const Vec1Int& typeMap);
/*******************************************************************************************
 * Make the type map for a whole set which contains several structures.
 *
 * @param Record containing an ML_AB file.
 * Calls makeTypeMap for every structure found in set
 *******************************************************************************************/
Vec2Int makeTypeMapSet(Record& set);
/*******************************************************************************************
 * Creating atom reorder maps for a whole ML_AB files.
 *
 * @param set containing ML_AB file
 * @param typeMaps contains type maps created by makeTypeMapSet to reorder atoms.
 *******************************************************************************************/
Vec2Int createReorderMap(const Record& set, const Vec2Int& typeMaps);
/*******************************************************************************************
 * Reordering atom, forces, positions, atom types and number atom per types for a whole set.
 *
 * @param structures is the structure entry from the ML_AB Record.
 * @param atomMap reordering map created by createReorderMap for this ML_AB file.
 * @param typeMap reordering map created by makeTypeMapSet for this ML_AB file.
 *******************************************************************************************/
void reorderStructures(Vec1ShRec& structures, const Vec2Int& atomMap, const Vec2Int& typeMap);
/*******************************************************************************************
 * Reorders the local reference configurations list written after the header in ML_AB.
 *
 * @param set contains a Record which has to have fields lrcAtom and lrcStructure.
 *******************************************************************************************/
void reorderLocalReferenceConfigurations(Record& set, const Vec2Int& atomMap);
/*******************************************************************************************
 * Reorder a whole record.
 *
 * @param atomMap reordering map for atoms created by createReorderMap
 * @param typeMap reordering map for types created by makeTypeMapSet
 *
 * Which means reordering positions, local reference configurations, reference forces,
 * type names and number atoms per type.
 *
 *******************************************************************************************/
void reorderRecord(Record& set, const Vec2Int& atomMap, const Vec2Int& typeMap);
/*******************************************************************************************
 * Reorder structures contained in ML_AB such that the order of types in all structures is
 * same as in type line of ML_AB file
 *
 * @param set set containing a set of structures and a line types according to which the
 * set will be reordered.
 *******************************************************************************************/
void reorderSet(Record& set);
/*******************************************************************************************
 * Creating a vector os Structure classes from an ML_AB Record.
 *
 * @param mlab Record containing ML_AB input file.
 *
 * Before Structures are created the Record is resorted such that the type order in
 * every Structure is the same as in the top type line in the ML_AB file.
 *******************************************************************************************/
std::vector<Structure> prepareStructures(Record& mlab);
/*******************************************************************************************
 * @brief Prepares a vector of Structure objects with atom indices distributed in a round-robin
 *fashion.
 *
 * This function first prepares a vector of Structure objects using the provided Record object.
 * It then assigns atom indices to each Structure object based on a round-robin distribution
 * across MPI ranks using the provided MlMPI instance.
 *
 * @param mlab Reference to a Record object used to prepare the initial structures.
 * @param mlmpi Shared pointer to an MlMPI instance for MPI rank and distribution information.
 * @return A vector of Structure objects with atom indices set for parallel processing.
 *******************************************************************************************/
std::vector<Structure> prepareStructures(Record& mlab, const std::shared_ptr<MlMPI>& mlmpi);
/*******************************************************************************************
 * @brief Extracts energy values from a given dataset record.
 *
 * This function retrieves a list of structures from the provided `Record` object,
 * extracts the energy values associated with each structure, and returns them
 * as a vector of real numbers.
 *
 * @param set The `Record` object containing the dataset from which energy values
 *            are to be extracted. It is expected to have a key "structures" that
 *            maps to a vector of shared pointers to records (`Vec1ShRec`).
 *
 * @return A vector of real numbers (`Vec1Real`) containing the energy values
 *         extracted from the structures in the dataset.
 *
 * @note The function assumes that each structure in the "structures" vector
 *       has a key "energy" that maps to a real number (`Real`).
 *
 * @throws std::runtime_error If the "structures" key is missing or if any structure
 *         does not contain the "energy" key.
 *
 * @code
 * Record datasetRecord = ...; // Initialize the dataset record
 * Vec1Real energies = dataset::extractEnergies(datasetRecord);
 * for (const auto& energy : energies) {
 *     std::cout << "Energy: " << energy << std::endl;
 * }
 * @endcode
 *******************************************************************************************/
Vec1Real extractEnergies(const Record& set);
/*******************************************************************************************
 * @brief Extracts energy values from a collection of structures.
 *
 * This function iterates over a collection of shared pointers to structures (`MLABstrucs`),
 * retrieves the "energy" property from each structure, and returns a vector of these energy values.
 *
 * @param MLABstrucs A vector of shared pointers to structures, each containing an "energy"
 *property.
 * @return Vec1Real A vector of energy values extracted from the input structures.
 *******************************************************************************************/
Vec1Real extractEnergies(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts the number of atoms from a dataset record.
 *
 * This function retrieves the "structures" field from the given dataset record,
 * iterates over the structures, and extracts the "numAtoms" field from each structure.
 *
 * @param set The dataset record containing the "structures" field.
 * @return A vector of integers, where each element represents the number of atoms
 *         in a corresponding structure within the dataset.
 *******************************************************************************************/
Vec1Int extractNumberAtoms(const Record& set);
/*******************************************************************************************
 * @brief Extracts the number of atoms from a collection of molecular structures.
 *
 * This function iterates over a collection of shared pointers to molecular structure records
 * and retrieves the "numAtoms" property from each record. The extracted values are returned
 * as a vector of integers.
 *
 * @param MLABstrucs A vector of shared pointers to constant `Record` objects, representing
 *molecular structures.
 * @return A vector of integers, where each element corresponds to the number of atoms in the
 *respective structure.
 *******************************************************************************************/
Vec1Int extractNumberAtoms(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts the number of atoms per type for all structures in the given dataset record.
 *
 * @param set The dataset record containing structural information.
 * @return A 2D vector (Vec2Int) where each row corresponds to the number of atoms per type for a
 *structure.
 *******************************************************************************************/
Vec2Int extractNumberAtomsPerType(const Record& set);
/*******************************************************************************************
 * @brief Extracts the number of atoms per type for a collection of structures.
 *
 * @param MLABstrucs A vector of shared pointers to constant records representing the structures.
 * @return A 2D vector (Vec2Int) where each row corresponds to the number of atoms per type for a
 *structure.
 *******************************************************************************************/
Vec2Int extractNumberAtomsPerType(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Counts the total number of atoms per type across all structures in the given dataset
 *record.
 *
 * @param set The dataset record containing structural information.
 * @return A 1D vector (Vec1Int) where each element corresponds to the total number of atoms for a
 *specific type.
 *******************************************************************************************/
Vec1Int countTotalAtomsPerType(const Record& set);
/*******************************************************************************************
 * @brief Extracts a specific component (x, y, or z) of atomic forces from a dataset record.
 *
 * This function retrieves the force data for all atoms in a dataset record and extracts the
 * specified component (x, y, or z) for each atom. The extracted force components are returned
 * as a 2D vector, where each inner vector corresponds to the forces of atoms in a structure.
 *
 * @param set The dataset record containing the structures and their associated force data.
 * @param xyzIndex The index of the force component to extract (0 for x, 1 for y, 2 for z).
 *                 Must be in the range [0, 2].
 *
 * @return A 2D vector (Vec2Real) containing the extracted force components for all atoms
 *         in all structures within the dataset record.
 *
 * @throws std::runtime_error If the xyzIndex is not in the range [0, 2].
 *
 * @note This function assumes that the dataset record contains a "structures" field, where
 *       each structure has a "forces" field (a 1D vector of size 3*numAtoms) and a "numAtoms"
 *       field indicating the number of atoms in the structure.
 *
 * @warning Ensure that the dataset record is properly formatted and contains the required
 *          fields; otherwise, undefined behavior may occur.
 *******************************************************************************************/
Vec2Real extractForceComponent(const Record& set, const Int xyzIndex);
/*******************************************************************************************
 * @brief Extracts a specific Cartesian component (x, y, or z) of atomic forces from a dataset.
 *
 * This function processes a collection of molecular structures and extracts the specified
 * Cartesian component (x, y, or z) of the atomic forces for each structure. The extracted
 * forces are returned as a 2D vector, where each inner vector corresponds to the forces
 * of a single structure.
 *
 * @param MLABstrucs A collection of molecular structures, represented as a vector of
 *                   shared pointers to constant `Record` objects.
 * @param xyzIndex   The index of the Cartesian component to extract (0 for x, 1 for y, 2 for z).
 *                   Must be in the range [0, 2].
 *
 * @return A 2D vector (`Vec2Real`) containing the extracted force components. Each inner
 *         vector corresponds to the forces of a single structure.
 *
 * @throws std::runtime_error If `xyzIndex` is not in the range [0, 2].
 *
 * @note The function assumes that each `Record` contains the following:
 *       - A vector of forces (`Vec1Real`) stored under the key "forces".
 *       - The number of atoms (`Int`) stored under the key "numAtoms".
 *
 * @warning This function performs a debug-level check to ensure `xyzIndex` is valid. If
 *          the check fails, an error message is logged, but the program may continue execution.
 *******************************************************************************************/
Vec2Real extractForceComponent(const Vec1ShCRec& MLABstrucs, const Int xyzIndex);
/*******************************************************************************************
 * @brief Extracts the X-component of forces from a given dataset record.
 *
 * This function retrieves the X-component (index 0) of the forces stored in the provided dataset
 *record.
 *
 * @param set The dataset record from which the X-component of forces will be extracted.
 * @return Vec2Real A vector containing the extracted X-component of forces.
 *******************************************************************************************/
Vec2Real extractForcesX(const Record& set);
/*******************************************************************************************
 * @brief Extracts the X-component of forces from the given MLAB structures.
 *
 * This function retrieves the X-component (index 0) of the forces from the provided
 * MLAB structures and returns them as a 2D vector of real numbers.
 *
 * @param MLABstrucs A collection of MLAB structures (Vec1ShCRec) to extract forces from.
 * @return Vec2Real A 2D vector containing the X-components of the forces.
 *******************************************************************************************/
Vec2Real extractForcesX(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts the Y-component of forces from a given dataset record.
 *
 * This function retrieves the Y-component (index 0) of the forces stored in the provided dataset
 *record.
 *
 * @param set The dataset record from which the Y-component of forces will be extracted.
 * @return Vec2Real A vector containing the extracted Y-component of forces.
 *******************************************************************************************/
Vec2Real extractForcesY(const Record& set);
/*******************************************************************************************
 * @brief Extracts the Y-component of forces from the given MLAB structures.
 *
 * This function retrieves the Y-component (index 1) of the forces from the provided
 * MLAB structures and returns them as a vector of real values.
 *
 * @param MLABstrucs A vector of shared pointers to MLAB structures (Vec1ShCRec).
 * @return Vec2Real A vector containing the Y-components of the forces.
 *******************************************************************************************/
Vec2Real extractForcesY(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts the Z-component of forces from a given dataset record.
 *
 * This function retrieves the Z-component (index 0) of the forces stored in the provided dataset
 *record.
 *
 * @param set The dataset record from which the Z-component of forces will be extracted.
 * @return Vec2Real A vector containing the extracted Z-component of forces.
 *******************************************************************************************/
Vec2Real extractForcesZ(const Record& set);
/*******************************************************************************************
 * @brief Extracts the Z-component of forces from the given dataset.
 *
 * This function retrieves the Z-component (index 2) of the forces from the provided
 * MLAB structures and returns them as a 2D vector of real numbers.
 *
 * @param MLABstrucs A vector of shared pointers to MLAB structures containing force data.
 * @return Vec2Real A 2D vector of real numbers representing the Z-component of the forces.
 *******************************************************************************************/
Vec2Real extractForcesZ(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts force vectors from a dataset record.
 *
 * This function retrieves the "structures" field from the given dataset record,
 * iterates over each structure, and extracts the "forces" field from each structure.
 * The extracted forces are returned as a 2D vector of real numbers.
 *
 * @param set The dataset record containing the "structures" field.
 * @return Vec2Real A 2D vector where each inner vector represents the forces
 *         associated with a structure in the dataset.
 *******************************************************************************************/
Vec2Real extractForces(const Record& set);
/*******************************************************************************************
 * @brief Extracts force vectors from a collection of MLAB structures.
 *
 * This function processes a collection of shared pointers to `Record` objects
 * and extracts the "forces" field from each record. The extracted forces are
 * returned as a 2D vector of real numbers.
 *
 * @param MLABstrucs A vector of shared pointers to constant `Record` objects,
 *                   representing the MLAB structures to process.
 * @return A 2D vector (`Vec2Real`) containing the extracted force vectors
 *         from the input MLAB structures.
 *******************************************************************************************/
Vec2Real extractForces(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts a specific stress component from a dataset.
 *
 * This function retrieves the stress values for a given component (specified by `xyzIndex`)
 * from the "structures" field of the provided dataset record. The stress values are extracted
 * from each structure in the dataset and returned as a vector.
 *
 * @param set The dataset record containing the "structures" field.
 * @param xyzIndex The index of the stress component to extract
 * (e.g., 0 for xx, 1 for xy, 2 for xz, 3 for yy, 4 for yz, 5 for zz).
 * @return A vector of stress values corresponding to the specified component for all structures.
 *******************************************************************************************/
Vec2Real extractStress(const Record& set);
/*******************************************************************************************
 * @brief Extracts the "stress" data from a collection of records.
 *
 * This function processes a vector of shared pointers to constant `Record` objects
 * and extracts the "stress" field from each record. The extracted data is returned
 * as a 2D vector of real numbers (`Vec2Real`).
 *
 * @param MLABstrucs A vector of shared pointers to constant `Record` objects
 *                   containing the data to extract.
 * @return A 2D vector (`Vec2Real`) containing the "stress" data extracted
 *         from the input records.
 *******************************************************************************************/
Vec2Real extractStress(const Vec1ShCRec& MLABstrucs);
/*******************************************************************************************
 * @brief Extracts a specific stress component from a dataset.
 *
 * This function retrieves the stress component specified by `xyzIndex` from the
 * "structures" field of the given dataset record. The stress components are
 * extracted for all structures in the dataset and returned as a vector.
 *
 * @param set The dataset record containing the "structures" field.
 * @param xyzIndex The index of the stress component to extract (0 to 5).
 *                 Valid indices correspond to the components of the stress tensor:
 *                 0 - σ_xx, 1 - σ_xy, 2 - σ_xz, 3 - σ_yy, 4 - σ_yz, 5 - σ_zz.
 * @return A vector of the extracted stress component values for all structures.
 *
 * @throws std::runtime_error If `xyzIndex` is not in the range [0, 5].
 *
 * @note This function assumes that each structure in the dataset contains a
 *       "stress" field, which is a vector of stress components.
 *******************************************************************************************/
Vec1Real extractStressComponent(const Record& set, const Int xyzIndex);
/*******************************************************************************************
 * @brief Extracts a specific stress component from a dataset of structures.
 *
 * This function retrieves the stress component specified by `xyzIndex` from each structure
 * in the provided dataset `MLABstrucs`. The stress components are stored in a vector and returned.
 *
 * @param MLABstrucs A vector of shared pointers to `Record` objects, each containing stress data.
 * @param xyzIndex The index of the stress component to extract (must be 0, 1, 2, 3, 4, or 5).
 * @return A vector (`Vec1Real`) containing the extracted stress component for each structure.
 *
 * @throws If `xyzIndex` is out of the valid range [0, 5], an error is logged via
 *`global::tutor.error`.
 *
 * @note The function assumes that each `Record` in `MLABstrucs` contains a stress vector
 *       accessible via the key `"stress"`.
 *******************************************************************************************/
Vec1Real extractStressComponent(const Vec1ShCRec& MLABstrucs, const Int xyzIndex);
/*******************************************************************************************
 * Creates a ML_FF Record from the setup (INCAR) Record.
 *
 * @param setup Record containing the setup (INCAR) information.
 *
 * Copies keys from the setup Record to the ML_FF Record. Only the tags specified in the
 * INCAR file are copied. If a tag is not found in the setup Record, an error is raised.
 *******************************************************************************************/
ShRec createMlffFromSetup(const Record& setup);
/*******************************************************************************************
 * Creates a ML_FF Record from the setup (INCAR) Record.
 *
 * @param setup Record containing the setup (INCAR) information.
 * @param mlff Record to which the needed INCAR input is copied to.
 *
 * Copies keys from the setup Record to the ML_FF Record. Only the tags specified in the
 * INCAR file are copied. If a tag is not found in the setup Record, an error is raised.
 *******************************************************************************************/
void createMlffFromSetup(const Record& setup, Record& mlff);
/*******************************************************************************************
 * Enters unused tags which are not used in the ML_AB file.
 *
 * These tags are needed to write ML_FF files but are not used during refit or
 * production. Routine will create a Record with these tags and return it.
 *******************************************************************************************/
ShRec createUnusedMlffTags();
/*******************************************************************************************
 * @brief Transfers specific data fields from an `mlab` record to an `mlff` record.
 *
 * This function extracts key information from the `mlab` record and populates
 * the corresponding fields in the `mlff` record. The transferred fields include
 * metadata about types, training structures, local reference configurations,
 * atomic masses, and other related properties.
 *
 * @param[in] mlab The source record containing the data to be transferred.
 * @param[out] mlff The destination record where the data will be stored.
 *
 * The following fields are transferred:
 * - "numTypes" from "maxTypes" in `mlab`
 * - "types" from "types" in `mlab`
 * - "training_structures" from "numStructures" in `mlab`
 * - "local_reference_cfgs" from "numLrc" in `mlab`
 * - "atomicMass" from "atomicMass" in `mlab`
 * - "numLrc" from "numLrc" in `mlab`
 *******************************************************************************************/
void addMlabToMlff(const Record& mlab, Record& mlff);
/*******************************************************************************************
 * @brief Computes statistical properties of a dataset and stores them in a record.
 *
 * This function calculates various statistical metrics (mean and variance) for
 * energies, forces, and stress components extracted from the input dataset (`mlab`)
 * and stores the results in the output record (`mlabStat`).
 *
 * @param[out] mlabStat A record where the computed statistical properties will be stored.
 *                      The following keys are populated:
 *                      - "average-energy-per-atom": Mean energy per atom.
 *                      - "variance-training-energies": Variance of the training energies.
 *                      - "variance-training-force": Variance of the training forces
 *                        along X, Y, and Z directions.
 *                      - "variance-training-stress": Variance of the training stress
 *                        components (6 components).
 * @param[in] mlab A record containing the input dataset from which energies, forces,
 *                 and stress components are extracted for statistical analysis.
 *******************************************************************************************/
void mlabStatistics(Record& mlabStat, const Record& mlab);
/*******************************************************************************************
 * @brief Rescales the energy values of all structures in the given dataset.
 *
 * This function adjusts the energy values of all structures within the provided
 * `Record` by multiplying them with the specified energy scaling factor.
 *
 * @param mlab The dataset record containing the structures to be rescaled.
 *             It must have a key "structures" that maps to a collection of `ShRec` objects.
 * @param energyFactor The scaling factor by which the energy values will be multiplied.
 *******************************************************************************************/
void rescaleEnergy(Record& mlab, const Real& energyFactor = 1.0);
/*******************************************************************************************
 * @brief Rescales the forces in a dataset by a given factor.
 *
 * This function iterates over all structures in the provided dataset (`mlab`)
 * and rescales their forces by multiplying each force component by the specified `forceFactor`.
 *
 * @param mlab        Reference to the dataset containing structures and their associated forces.
 * @param forceFactor Scaling factor to be applied to the forces.
 *
 * The function assumes that:
 * - The dataset (`mlab`) contains a key "structures" that maps to a collection of structures
 *(`Vec1ShRec`).
 * - Each structure contains a key "forces" that maps to a collection of force components
 *(`Vec1Real`).
 *******************************************************************************************/
void rescaleForces(Record& mlab, const Real& forceFactor = 1.0);
/**
 * @brief Rescales the stress values of all structures in a dataset by a given factor.
 *
 * This function iterates over all structures in the provided dataset record (`mlab`)
 * and rescales their stress values by multiplying each stress component by the specified
 * `stressFactor`.
 *
 * @param mlab Reference to the dataset record containing the structures and their stress data.
 * @param stressFactor The factor by which to scale the stress values.
 *
 * The function assumes that:
 * - The dataset record (`mlab`) contains a key "structures" that maps to a collection of
 * structures.
 * - Each structure contains a key "stress" that maps to a collection of stress values.
 */
void rescaleStress(Record& mlab, const Real& forceFactor = 1.0);
/*******************************************************************************************
 * @brief Computes the mean energy from a dataset.
 *
 * This function extracts energy values from the provided dataset and calculates their mean.
 *
 * @param data A dataset containing energy-related records.
 * @return The mean energy value as a Real.
 *******************************************************************************************/
Real meanEnergy(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the mean force from the given dataset.
 *
 * This function extracts force values from the provided dataset, converts them
 * into a 1D vector, and calculates their mean.
 *
 * @param data A constant reference to the dataset containing force records.
 * @return The mean value of the forces as a Real.
 *******************************************************************************************/
Real meanForce(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the mean stress from the given dataset.
 *
 * This function extracts stress values from the input dataset, converts them
 * to a 1D vector, and calculates the mean stress using statistical tools.
 *
 * @param data A constant reference to the input dataset of type Vec1ShCRec.
 * @return The mean stress as a Real value.
 *******************************************************************************************/
Real meanStress(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the mean values of energy, force, and stress from the given dataset.
 *
 * This function calculates the mean energy, mean force, and mean stress
 * for the provided dataset and returns them as a tuple.
 *
 * @param data A constant reference to the dataset (Vec1ShCRec) to process.
 * @return A tuple containing the mean energy, mean force, and mean stress
 *         in the order (meanEnergy, meanForce, meanStress).
 *******************************************************************************************/
std::tuple<Real, Real, Real> meanDatsets(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the variance of energy values in the given dataset.
 *
 * This function extracts energy values from the provided dataset, calculates
 * their mean, and then computes the variance of the energy values.
 *
 * @param data A dataset containing energy-related records.
 * @return The variance of the energy values in the dataset.
 *******************************************************************************************/
Real varianceEnergy(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the variance of forces extracted from the given dataset.
 *
 * This function extracts force values from the input dataset, calculates their mean,
 * and then computes the variance of the forces.
 *
 * @param data A dataset containing force-related records.
 * @return The variance of the extracted forces.
 *******************************************************************************************/
Real varianceForce(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the variance of stress values from the given dataset.
 *
 * This function extracts stress values from the input dataset, calculates their mean,
 * and then computes the variance of the stress values.
 *
 * @param data A constant reference to the input dataset of type Vec1ShCRec.
 * @return The variance of the stress values as a Real.
 *******************************************************************************************/
Real varianceStress(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the variances of energy, force, and stress from the given dataset.
 *
 * This function calculates the variance for three key properties (energy, force, and stress)
 * from the provided dataset and returns them as a tuple.
 *
 * @param data A constant reference to the dataset (`Vec1ShCRec`) from which variances are computed.
 * @return A tuple containing the variances of energy, force, and stress, in that order.
 *******************************************************************************************/
std::tuple<Real, Real, Real> varianceDatsets(const Vec1ShCRec& data);
/*******************************************************************************************
 * @brief Computes the Root Mean Square Error (RMSE) of energy values between predicted and
 *reference datasets.
 *
 * This function extracts energy values from the given predicted and reference datasets,
 * and calculates the RMSE to quantify the difference between the two sets of energy values.
 *
 * @param predicted The predicted dataset containing energy values.
 * @param reference The reference dataset containing energy values.
 * @return The RMSE of the energy values as a Real number.
 *******************************************************************************************/
Real rmseEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Root Mean Square Error (RMSE) of forces between predicted and reference
 *datasets.
 *
 * This function extracts force vectors from the predicted and reference datasets, converts them to
 *1D vectors, and calculates the RMSE to quantify the difference between the two sets of forces.
 *
 * @param predicted The predicted dataset containing force information.
 * @param reference The reference dataset containing force information.
 * @return The RMSE value as a measure of the difference between predicted and reference forces.
 *******************************************************************************************/
Real rmseForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Root Mean Square Error (RMSE) between predicted and reference stress values.
 *
 * This function extracts stress values from the given predicted and reference datasets,
 * converts them to 1D vectors, and calculates the RMSE to quantify the difference.
 *
 * @param predicted The predicted dataset containing stress values (Vec1ShCRec format).
 * @param reference The reference dataset containing stress values (Vec1ShCRec format).
 * @return The RMSE value as a Real type, representing the error between predicted and reference
 *stress values.
 *******************************************************************************************/
Real rmseStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Root Mean Square Error (RMSE) for energy, force, and stress between predicted
 *and reference datasets.
 *
 * @param predicted A vector of predicted dataset records.
 * @param reference A vector of reference dataset records.
 * @return A tuple containing the RMSE values for energy, force, and stress, respectively.
 *******************************************************************************************/
std::tuple<Real, Real, Real> rmseDatsets(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Mean Absolute Error (MAE) of energy values between predicted and reference
 *datasets.
 *
 * This function extracts energy values from the given predicted and reference datasets,
 * and calculates the Mean Absolute Error (MAE) between them.
 *
 * @param predicted A vector of predicted records containing energy values.
 * @param reference A vector of reference records containing energy values.
 * @return The computed MAE of the energy values as a Real number.
 *******************************************************************************************/
Real maeEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Mean Absolute Error (MAE) of forces between predicted and reference datasets.
 *
 * This function extracts force vectors from the predicted and reference datasets, converts them
 * to 1D vectors, and calculates the Mean Absolute Error (MAE) between the two.
 *
 * @param predicted The predicted dataset containing force vectors.
 * @param reference The reference dataset containing force vectors.
 * @return The computed MAE as a Real value.
 *******************************************************************************************/
Real maeForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Mean Absolute Error (MAE) of stress values between predicted and reference
 *datasets.
 *
 * This function extracts stress values from the given predicted and reference datasets, converts
 *them to 1D vectors, and calculates the MAE between the two.
 *
 * @param predicted The predicted dataset containing stress values (Vec1ShCRec format).
 * @param reference The reference dataset containing stress values (Vec1ShCRec format).
 * @return The computed Mean Absolute Error (MAE) as a Real value.
 *******************************************************************************************/
Real maeStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Mean Absolute Error (MAE) for energy, force, and stress between predicted and
 *reference datasets.
 *
 * This function calculates the MAE for three physical quantities: energy, force, and stress,
 * by comparing the predicted dataset with the reference dataset. The results are returned as a
 *tuple.
 *
 * @param predicted The predicted dataset (Vec1ShCRec).
 * @param reference The reference dataset (Vec1ShCRec).
 * @return A tuple containing the MAE for energy, force, and stress, in that order.
 *******************************************************************************************/
std::tuple<Real, Real, Real> maeDatsets(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Pearson correlation coefficient between predicted and reference energy
 *datasets.
 *
 * This function extracts energy values from the given predicted and reference datasets
 * and calculates the Pearson correlation coefficient to measure the linear relationship
 * between the two datasets.
 *
 * @param predicted A vector of predicted records containing energy data.
 * @param reference A vector of reference records containing energy data.
 * @return The Pearson correlation coefficient as a Real value.
 *******************************************************************************************/
Real pearsonCorrEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Pearson correlation coefficient between predicted and reference force
 *datasets.
 *
 * This function extracts force components from the given predicted and reference datasets,
 * converts them into 1D vectors, and calculates the Pearson correlation coefficient to
 * measure the linear relationship between the two datasets.
 *
 * @param predicted The predicted dataset containing force components (Vec1ShCRec format).
 * @param reference The reference dataset containing force components (Vec1ShCRec format).
 * @return The Pearson correlation coefficient (Real) between the predicted and reference forces.
 *******************************************************************************************/
Real pearsonCorrForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Pearson correlation coefficient between predicted and reference stress
 *values.
 *
 * This function extracts stress values from the given predicted and reference datasets,
 * converts them to 1D vectors, and calculates the Pearson correlation coefficient to
 * measure the linear relationship between the two sets of stress values.
 *
 * @param predicted A vector of predicted stress data (Vec1ShCRec format).
 * @param reference A vector of reference stress data (Vec1ShCRec format).
 * @return The Pearson correlation coefficient (Real) between the predicted and reference stress
 *values.
 *******************************************************************************************/
Real pearsonCorrStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the Pearson correlation coefficients for energy, force, and stress between two
 *datasets.
 *
 * This function calculates the Pearson correlation coefficients for three different properties:
 * energy, force, and stress, by comparing the predicted dataset with the reference dataset.
 *
 * @param predicted The predicted dataset (Vec1ShCRec).
 * @param reference The reference dataset (Vec1ShCRec).
 * @return A tuple containing the Pearson correlation coefficients for energy, force, and stress,
 *respectively.
 *******************************************************************************************/
std::tuple<Real, Real, Real> pearsonCorrDatsets(const Vec1ShCRec& predicted,
                                                const Vec1ShCRec& reference);
/*******************************************************************************************
 * @brief Computes the energy shift for an isolated atom based on reference energies and ab initio
 *calculations.
 *
 * This function calculates the energy shift for an isolated atom by comparing the ab initio
 *energies of molecular structures with the reference energies derived from the number of atoms of
 *each type. The energy shift is normalized by the total number of atoms in each structure and
 *averaged across all structures.
 *
 * @param mlab The record containing molecular structure data and ab initio energies.
 * @param atomRefEnergy A vector of reference energies for each atom type.
 * @return The computed energy shift for the isolated atom.
 *******************************************************************************************/
Real computeEnergyShiftIsolatedAtom(const Record& mlab, const Vec1Real& atomRefEnergy);
/*******************************************************************************************
 * @brief Computes the average energy shift between ab initio energies and reference energies for
 *training configurations.
 *
 * This function calculates the average energy shift by comparing the ab initio energies of
 *structures in the dataset with their corresponding reference energies derived from atomic
 *reference values. The energy shift is computed for each configuration and averaged over all
 *configurations in the dataset.
 *
 * @param mlab The dataset record containing structural information and ab initio energies.
 * @param atomRefEnergy A vector of reference energies for each atom type.
 * @return The average energy shift across all training configurations.
 *
 * @details
 * - The function extracts the list of structures and their ab initio energies from the dataset
 *record.
 * - For each structure, it computes the reference energy based on the number of
 *   atoms per type and their corresponding reference energies.
 * - The energy shift is calculated as the difference between the ab initio energy and the reference
 *energy for each configuration.
 * - The average energy shift is returned by dividing the total energy shift by the number of
 *configurations.
 *******************************************************************************************/
Real computeEnergyShiftAverageTrainEnergy(const Record& mlab, const Vec1Real& atomRefEnergy);
/*******************************************************************************************
 * @brief Prepares total energy for refit with average isolated atom reference.
 *
 * This function computes the energy refit for a set of structures by subtracting the reference
 *energy contribution of isolated atoms and applying an energy shift. The adjusted energy values are
 *normalized by the number of atoms in each structure.
 *
 * @param mlab The record containing crystal structure data and associated properties.
 * @param atomRefEnergy A vector of reference energies for isolated atoms, indexed by atom type.
 * @return A vector of refitted energy values for the molecular structures.
 *
 * @details
 * - Extracts molecular structure data and energy values from the input record.
 * - Computes the energy shift for isolated atoms using the provided reference energies.
 * - Iterates over each molecular structure to:
 *   - Calculate the reference energy contribution based on the number of atoms per type.
 *   - Adjust the energy value by subtracting the reference energy and normalizing by the total
 *number of atoms.
 *   - Applies the computed energy shift.
 *******************************************************************************************/
Vec1Real prepareEnergyRefitIsolatedAtom(const Record& mlab, const Vec1Real& atomRefEnergy);
/*******************************************************************************************
 * @brief Prepares a vector of TypeMap objects based on the provided Record.
 *
 * This function extracts structure and type information from the given Record
 * and uses it to populate a vector of TypeMap objects. Each TypeMap is updated
 * with the global type names and the specific type names associated with each structure.
 *
 * @param mlab The Record object containing the data for structures and type names.
 *             It is expected to have the following keys:
 *             - "structures": A vector of shared pointers to records (Vec1ShRec).
 *             - "typeNames": A vector of strings (Vec1String) representing all type names.
 *
 * @return A vector of TypeMap objects, where each TypeMap corresponds to a structure
 *         in the input Record and is updated with the relevant type information.
 *******************************************************************************************/
std::vector<TypeMap> prepareTypeMaps(const Record& mlab);

} // namespace vaspml::dataset

#endif
