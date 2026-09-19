#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include "BasisFunctions.hpp"
#include "Random.hpp"
#include "Structure.hpp"
#include "types.hpp"

namespace vaspml
{

class MlMPI;
struct Record;
template<class T>
class ShmemArray2D;

enum class SelectionMethod
{
    FarthestPointSamplingN,
    RandomSelection,
    FarthestPointSamplingThreshold,
    SequentialLikelihoodCover,
    None
};

enum class ParallelMode
{
    OpenMP,
    MPI,
    None
};

enum class SelectionMetric
{
    L1Norm,
    L2Norm,
    LpNorm,
    PolyKernelNorm,
    antiPearsonNorm2,
    None
};

enum class DataDist
{
    roundRobin,
    block,
    None
};

/*******************************************************************************************
 * @brief Chooses a selection method based on the input configuration.
 *
 * This function determines the selection method to be used by interpreting the
 * "ML_SALGO" parameter from the provided `incar` object. The supported selection
 * methods are:
 * - "FPSN": Farthest Point Sampling (N variant).
 * - "RANN": Random Selection.
 * - "FPST": Farthest Point Sampling (Threshold variant).
 *
 * If the "ML_SALGO" parameter is not recognized, an error message is logged.
 *
 * @param incar The input configuration object containing the "ML_SALGO" parameter.
 * @return The selected `SelectionMethod`. Returns `SelectionMethod::None` if the
 *         method is not recognized.
 *******************************************************************************************/
SelectionMethod chooseSelectionMethod(const ShRec& incar);
/*******************************************************************************************
 * @brief Selects a metric based on the input configuration.
 *
 * This function determines the appropriate selection metric based on the
 * value of the "ML_SMETRIC" key in the provided `ShRec` object. If the
 * specified metric is not recognized, an error message is logged.
 *
 * @param incar The input configuration object (`ShRec`) containing the
 *              "ML_SMETRIC" key.
 * @return The selected metric as a `SelectionMetric` enum value. Returns
 *         `SelectionMetric::None` if the metric is not recognized.
 *
 * @note Recognized values for "ML_SMETRIC" and their corresponding metrics:
 *       - "l1norm" -> `SelectionMetric::L1Norm`
 *       - "l2norm" -> `SelectionMetric::L2Norm`
 *       - "lpnorm" -> `SelectionMetric::LpNorm`
 *       - "pkernel4" -> `SelectionMetric::PolyKernelNorm`
 *       - "pcorr" -> `SelectionMetric::antiPearsonNorm2`
 *
 * @warning If the "ML_SMETRIC" value is not recognized, an error message
 *          is logged, and the function returns `SelectionMetric::None`.
 *******************************************************************************************/
SelectionMetric chooseMetric(const ShRec& incar);

class Selector
{
  public:
    /*******************************************************************************************
     * @brief Constructs a Selector object and initializes its members.
     *
     * This constructor initializes the Selector object with the provided input, mlab, and data
     *records. It sets up basis functions and prepares structures based on the input and mlab
     *records. Additionally, it configures the local configuration parameters and adds a
     *"descriptors" entry to the data record.
     *
     * @param input The input record (ShRec) containing configuration and parameters.
     * @param mlab The mlab record (ShRec) used for preparing structures and other settings.
     * @param data The data record (ShRec) to be used or initialized if null.
     *******************************************************************************************/
    Selector(const ShRec&                  input,
             const ShRec&                  mlab,
             const ShRec&                  data = nullptr,
             const std::shared_ptr<MlMPI>& mlmpi = nullptr);
    /*******************************************************************************************
     * @brief Constructs a Selector object with the given data.
     *
     * This constructor initializes the `incar` and `mlab` members as shared pointers
     * to new `Record` objects. The `data` member is initialized to the provided `data`
     * parameter if it is not null; otherwise, it is initialized to a new `Record` object.
     * The `selectionMethod` is set to `SelectionMethod::None`.
     *
     * @param data A shared pointer to a `Record` object. If `nullptr`, a new `Record` is created.
     *******************************************************************************************/
    Selector(const ShRec& data = nullptr);
    /*******************************************************************************************
     * @brief Initializes the Selector object with the provided input data and configurations.
     *
     * This method sets up the Selector object by initializing its internal data members,
     * preparing basis functions, structures, and configuring the selection method.
     *
     * @param input The input configuration record (ShRec) containing initialization parameters.
     * @param mlab The machine learning auxiliary record (ShRec) used for preparing structures.
     * @param data The data record (ShRec) to be used. If nullptr, a new shared Record is created.
     * @param mlmpi A shared pointer to the MlMPI object for parallel processing support.
     *******************************************************************************************/
    void init(const ShRec&                  input,
              const ShRec&                  mlab,
              const ShRec&                  data = nullptr,
              const std::shared_ptr<MlMPI>& mlmpi = nullptr);
    /*******************************************************************************************
     * @brief Executes the selection process based on the specified selection method.
     *
     * This function performs the selection of atoms or configurations using the
     * specified selection method. It supports multiple selection methods, including
     * random selection and farthest point sampling. The function also computes
     * descriptors and manages atom lists as needed.
     *
     * @details
     * - If no selection method is chosen (`SelectionMethod::None`), an error is logged.
     * - For non-random selection methods, descriptors are computed before proceeding.
     * - Atom lists are created and processed for each type up to the maximum number of types.
     * - Depending on the selection method:
     *   - `SelectionMethod::RandomSelection`: Randomly selects atoms and adds them to the list.
     *   - `SelectionMethod::FarthestPointSampling`: Reorders descriptors, performs farthest point
     *sampling, and adds the results to the list.
     *
     * @note The function logs the start and end of the selection process for debugging purposes.
     *
     * @throws Error if no selection method is chosen.
     *******************************************************************************************/
    void select();
    /*******************************************************************************************
     * @brief Writes the current MLab data to a file.
     *
     * This function saves the MLab data associated with the Selector object
     * to the specified file in MLab format.
     *
     * @param filename The name of the file to write the MLab data to.
     *******************************************************************************************/
    void writeMlabFile(const String& filename) const;
    /*******************************************************************************************
     * @brief Sets up the data distribution for the selector.
     *
     * This function initializes the `centralIndexes` vector to match the size of the `structures`
     *vector. For each structure, it determines the data distribution based on the number of atoms
     *in the structure and stores the result in the corresponding index of `centralIndexes`.
     *******************************************************************************************/
    void setUpDataDistribution();
    /*******************************************************************************************
     * @brief Retrieves the current value of the `mlab` member.
     *
     * @return ShCRec The value of the `mlab` member.
     *******************************************************************************************/
    const ShCRec get_mlab() const;

  private:
    /*******************************************************************************************
     * @brief Collects distributed data based on the current parallel mode.
     *
     * This function determines the parallel mode and invokes the appropriate
     * method to collect distributed data. It supports OpenMP and MPI modes.
     *
     * - If the parallel mode is `ParallelMode::OpenMP`, it calls `collectDistributedDataOMP()`.
     * - If the parallel mode is `ParallelMode::MPI`, it calls `collectDistributedDataShmemArray()`.
     *******************************************************************************************/
    void collectDistributedData();
    /*******************************************************************************************
     * @brief Collects distributed data across MPI processes using OpenMP.
     *
     * This function gathers distributed data from multiple MPI processes and consolidates it
     * into a single structure. It ensures that the data is reordered to match the original
     * indexing. The function uses OpenMP for parallel processing and relies on shared memory
     * for efficient data handling.
     *
     * @details
     * - The function starts by initializing a timer for performance measurement.
     * - If the MPI environment (`mlmpi`) is available, it iterates over all structures:
     *   - Retrieves the shared record (`ShRec`) for descriptors.
     *   - Gathers the central indexes and descriptors from all MPI processes using
     *`rec::mpi::gatherV`.
     *   - If the current process is the root (rank 0), it:
     *     - Updates the shared record with the gathered data.
     *     - Computes the new order of indices to undo the round-robin distribution.
     *     - Reorders the descriptor data (`Vec2Real`) to match the original indexing.
     * - Stops the timer after completing the data collection.
     *
     * @note This function assumes that the data being gathered is of type `Vec2Real` and uses
     *       specific tools like `vector_tools::invertIndex` and
     *`vector_tools::inplaceReorderVector` for reordering.
     *******************************************************************************************/
    void collectDistributedDataOMP();
    /*******************************************************************************************
     * @brief Collects distributed data into a shared memory array.
     *
     * This function gathers distributed data from multiple processes using MPI,
     * reorders the data to undo round-robin distribution, and prepares it for
     * shared memory storage. The collected data is stored in the `descriptors`
     * shared memory array.
     *
     * - Starts and stops a timer for performance measurement.
     * - Uses MPI to gather distributed data.
     * - Reorders the gathered data based on the central indexes.
     * - Allocates and fills shared memory for descriptors.
     *
     * @note This function is only executed if `mlmpi` is not null.
     *******************************************************************************************/
    void collectDistributedDataShmemArray();
    /*******************************************************************************************
     * @brief Computes descriptors for a set of structures.
     *
     * This function is responsible for computing descriptors for a collection of structures.
     * It initializes the necessary data structures, processes each structure, and collects
     * the computed descriptors. The descriptors are stored in a shared data object for further use.
     *
     * The function performs the following steps:
     * - Initializes the `Frame` object with the required parameters.
     * - Prepares the data storage for descriptors.
     * - Iterates over the structures and computes descriptors for each structure.
     * - Collects and rearranges the computed descriptors using a `DescriptorCollector`.
     *
     * @note This function assumes that the `mlab`, `incar`, `data`, and `structures` members
     *       of the `Selector` class are properly initialized before calling this function.
     *
     * @warning Ensure that the `structures` vector is not empty before calling this function,
     *          as it iterates over the structures to compute descriptors.
     *******************************************************************************************/
    void computeDescriptors();
    /*******************************************************************************************
     * @brief Reorders descriptors for a specific type index in the Selector class.
     *
     * This function processes and reorganizes descriptor data for a given type index.
     * It iterates over the structures and their associated atoms, extracting and
     * reordering descriptor values based on the provided type index. The reordered
     * descriptors are stored in the `typeOrderDesc` member variable.
     *
     * @param typeIndex The index of the type for which descriptors are to be reordered.
     *
     * @details
     * - The function initializes descriptor offsets for each key in `constants::descriptorKeyList`.
     * - It retrieves the number of atoms per type for each structure and processes the descriptors
     *   for the specified type index.
     * - Descriptor values are extracted and stored in the `typeOrderDesc` matrix, which is indexed
     *   by the total atom count and descriptor index.
     * - The function assumes that `resizeTypeOrderDesc`, `structures`, `data`, `descSizes`,
     *   `typeOrderDesc`, and `descOffset` are properly initialized and accessible within the class.
     *
     * @note The function contains commented-out code for debugging and file output, which is
     *currently not active.
     *
     * @todo Implement the conversion of `typeIndex` to a subtype index.
     *******************************************************************************************/
    void reorderDescriptors(const Size typeIndex);
    /*******************************************************************************************
     * @brief Resizes the `typeOrderDesc` container based on the number of atoms of a specific type
     *and the total number of descriptors.
     *
     * This function calculates the total number of atoms for a given type index across all
     *structures and resizes the `typeOrderDesc` container accordingly. It also determines the total
     *number of descriptors and resizes each element of `typeOrderDesc` to accommodate the
     *descriptors.
     *
     * @param typeIndex The index of the type for which the `typeOrderDesc` container is resized.
     *
     * The function performs the following steps:
     * 1. Computes the total number of atoms for the specified type index across all structures.
     * 2. Resizes the `typeOrderDesc` container to match the total number of atoms.
     * 3. Retrieves the total number of descriptors from the descriptor key list.
     * 4. Resizes each element of `typeOrderDesc` to accommodate the computed number of descriptors.
     *******************************************************************************************/
    void resizeTypeOrderDesc(const Size typeIndex);
    /*******************************************************************************************
     * @brief Resizes the index arrays based on the maximum number of types.
     *
     * This function retrieves the maximum number of types (`maxTypes`) from the `mlab` object
     * and resizes the `strucIndex` and `atomIndexInStruc` arrays to match this value.
     *******************************************************************************************/
    void resizeIndexArrays();
    /*******************************************************************************************
     * @brief Allocates shared memory for storing descriptors.
     *
     * This function calculates the total number of atoms per type and the total
     * number of descriptors. It then allocates shared memory arrays for each
     * atom type to store the descriptors.
     *
     * - Logs the start and end of the allocation process if the rank is 0.
     * - Uses `ShmemArray2D` for shared memory allocation.
     *
     * @note This function is only executed if `mlmpi` is not null.
     *******************************************************************************************/
    void allocate_descriptorsShmem();
    /*******************************************************************************************
     * @brief Fills the shared memory arrays with descriptor data.
     *
     * This function populates the shared memory arrays with descriptor data
     * for each atom type. It iterates over the structures, retrieves the
     * descriptors, and writes them into the shared memory arrays.
     *
     * - Handles data for multiple atom types and structures.
     * - Uses offsets and descriptor sizes to correctly place data in shared memory.
     * - Resets the descriptor records after processing.
     * - Synchronizes processes using an MPI barrier.
     *
     * @note This function is only executed if `mlmpi` is not null.
     *******************************************************************************************/
    void fill_descriptorsShmem();
    /*******************************************************************************************
     * @brief Performs a random selection and updates the corresponding list.
     *
     * This function generates a random selection based on the specified type
     * and updates the associated list in the Mlab system.
     *
     * @param type The type identifier used to determine the selection criteria.
     *******************************************************************************************/
    void randomSelection(const Size type);
    /*******************************************************************************************
     * @brief Randomly selects a specified number of samples of a given type.
     *
     * This function selects a random subset of samples from a collection based on the specified
     *type. It uses a random number generator to draw indices within the range of available samples
     *for the given type.
     *
     * @param type The type of samples to select from. This corresponds to an index in the
     *`nLocConf` and `strucIndex` arrays.
     *
     * The function performs the following steps:
     * - Logs the start of the selection process and the number of samples to be selected.
     * - Reinitializes the random number generator to the range of available indices for the
     *specified type.
     * - Draws a random subset of indices based on the number of samples specified in
     *`nLocConf[type]`.
     * - Logs the completion of the selection process.
     *******************************************************************************************/
    void makeRandomSelection(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling for a given number of points.
     *
     * This function executes the farthest point sampling algorithm by reordering descriptors,
     * performing the sampling, and adding the results to MATLAB for further processing.
     *
     * @param type The number of points to sample or a parameter specifying the sampling type.
     *******************************************************************************************/
    void farthestPointSamplingNPoints(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling to select a specified number of samples.
     *
     * This function selects `nLocRef` samples using the farthest point sampling (FPS) algorithm.
     * The selected samples and their corresponding indices are stored in `descSelection` and
     * `indexSelection`, respectively. The function also outputs the progress to the console.
     *
     * @param nLocRef The number of samples to select.
     *
     * @note The FPS algorithm is invoked through the `fps::farthestPointSampling` function.
     *       The detailed indices of the selected configurations are not printed but can be
     *       uncommented for debugging purposes.
     *******************************************************************************************/
    void makeFarthestPointSamplingNPointsOMP(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling to select a specified number of points using MPI.
     *
     * This function selects a specified number of points from a dataset using the farthest point
     *sampling (FPS) algorithm. The selection is based on the specified metric (e.g., L1 norm, L2
     *norm, polynomial kernel norm, or anti-Pearson norm). The function supports parallel execution
     *using MPI.
     *
     * @param type The type of configuration to process, which determines the dataset and number of
     *points to sample.
     *
     * The function performs the following steps:
     * - Logs the start of the sampling process (only on rank 0 if MPI is enabled).
     * - Initializes a random seed for the sampling process.
     * - Selects the appropriate metric function based on the `selectionMetric` member variable.
     * - Calls the `farthestPointSamplingKPointsMPI` function to perform the sampling.
     * - Logs an error if an unsupported configuration is encountered (e.g., unsupported polynomial
     *kernel power).
     * - Stops the timer and logs the completion of the sampling process (only on rank 0 if MPI is
     *enabled).
     *
     * @note This function assumes that the `mlmpi` object is properly initialized for MPI
     *communication.
     * @note The `descriptorsShmem` and `nLocConf` member variables must be properly initialized
     *before calling this function.
     * @note The `selectionMetric` member variable determines the metric used for sampling.
     *
     * @throws std::runtime_error If an unsupported configuration is encountered (e.g., unsupported
     *polynomial kernel power).
     *******************************************************************************************/
    void makeFarthestPointSamplingNPointsMPI(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling with a threshold for the given type.
     *
     * This function reorders descriptors, applies farthest point sampling with a threshold,
     * and adds the resulting list to Mlab for the specified type.
     *
     * @param type The type of sampling to be performed.
     *******************************************************************************************/
    void farthestPointSamplingThreshold(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling based on a specified threshold and selection metric.
     *
     * This function selects a subset of points from the input data using the farthest point
     *sampling (FPS) algorithm. The selection is based on a predefined threshold and the specified
     *distance metric. The function supports multiple distance metrics, including L1 norm, L2 norm,
     *polynomial kernel norm, and anti-Pearson norm.
     *
     * @param type The size/type parameter used to determine the sampling configuration.
     *
     * The function operates as follows:
     * - A threshold value (`tresh`) is set to 0.999.
     * - Depending on the `selectionMetric`, the appropriate distance metric is used for FPS.
     * - The selected indices are stored in `indexSelection`.
     *
     * Supported selection metrics:
     * - `SelectionMetric::L1Norm`: Uses the L1 norm for distance calculation.
     * - `SelectionMetric::L2Norm`: Uses the L2 norm for distance calculation.
     * - `SelectionMetric::PolyKernelNorm`: Uses a polynomial kernel norm with a fixed power of 4.
     * - `SelectionMetric::antiPearsonNorm2`: Uses the anti-Pearson norm for distance calculation.
     *
     * @note The function outputs progress messages to the console.
     *******************************************************************************************/
    void makeFarthestPointSamplingThresholdOMP(const Size type);
    /*******************************************************************************************
     * @brief Performs farthest point sampling with a distance threshold using MPI.
     *
     * This function selects indices based on a specified distance threshold and
     * the chosen selection metric. It supports multiple distance metrics, including
     * L1 norm, L2 norm, polynomial kernel norm, and anti-Pearson norm. The function
     * is parallelized using MPI and ensures that the selection process is consistent
     * across ranks.
     *
     * @param type The type of descriptor to use for sampling, which determines the
     *             threshold and the shared memory descriptor.
     *
     * @note The function logs the start and end of the process on rank 0. It also
     *       validates the polynomial kernel norm power and throws an error if it
     *       is not supported.
     *
     * @throws std::runtime_error If the polynomial kernel norm power is not supported.
     *******************************************************************************************/
    void makeFarthestPointSamplingThresholdMPI(const Size type);
    /*******************************************************************************************
     * @brief Executes the sequential likelihood cover process for a given type.
     *
     * This function performs the following steps:
     * 1. Reorders descriptors based on the specified type.
     * 2. Creates a sequential likelihood cover for the given type.
     * 3. Adds the resulting list to Mlab.
     *
     * @param type The type identifier used to guide the sequential likelihood cover process.
     *******************************************************************************************/
    void sequentialLikelihoodCover(const Size type);
    /*******************************************************************************************
     * @brief Constructs a sequential likelihood coverage selection based on the specified selection
     *metric.
     *
     * This method selects indices based on the provided selection metric and type order descriptor.
     * It supports various metrics such as L1 norm, L2 norm, polynomial kernel norm, and
     *anti-Pearson norm. The selected indices are stored in the `indexSelection` member variable.
     *
     * Supported selection metrics:
     * - `SelectionMetric::L1Norm`: Uses the L1 norm for selection.
     * - `SelectionMetric::L2Norm`: Uses the L2 norm for selection.
     * - `SelectionMetric::PolyKernelNorm`: Uses a polynomial kernel norm with a fixed power of 4.
     * - `SelectionMetric::antiPearsonNorm2`: Uses the anti-Pearson norm for selection.
     *
     * @note The method outputs a completion message to the standard output upon finishing.
     *******************************************************************************************/
    void makeSequentialLikelihoodCover(const Size type);
    /*******************************************************************************************
     * @brief Populates atom lists by mapping atom indices to their respective structures and types.
     *
     * This function iterates over all structures and their atom types, creating mappings between
     * atom types, structure indices, and atom indices within each structure. The mappings are
     *stored in `strucIndex` and `atomIndexInStruc` for later use.
     *
     * The function performs the following steps:
     * 1. Resizes the index arrays to ensure they are ready for population.
     * 2. Iterates over all structures and retrieves the number of atoms per type and the type
     *mapping.
     * 3. For each atom type in a structure, maps the atom indices to their corresponding structure
     *    and type, while accounting for offsets within the structure.
     *
     * @note The commented-out section at the end of the function can be used for debugging purposes
     *       to print the mappings of types, structures, and atom indices.
     *******************************************************************************************/
    void makeAtomLists();
    /*******************************************************************************************
     * @brief Adds a list of selected indices to the `mlab` data structure.
     *
     * This function updates the `mlab` object with information about the selected indices
     * for a given type. It modifies the `numLrc`, `lrcStructure`, and `lrcAtom` fields
     * in the `mlab` data structure for the specified `typeIndex`.
     *
     * @param typeIndex The index representing the type for which the data is being added.
     *
     * The function performs the following steps:
     * - Updates the `numLrc` field with the size of the `indexSelection` list.
     * - Resizes the `lrcStructure` and `lrcAtom` vectors to match the size of `indexSelection`.
     * - Populates the `lrcStructure` and `lrcAtom` vectors with the corresponding structure
     *   and atom indices, incremented by 1.
     *******************************************************************************************/
    void addListToMlab(const Size typeIndex);
    /*******************************************************************************************
     * Choose data distribution for descriptor computation with MPI.
     *
     * Available distribution schemes are round robin and block distribution.
     *******************************************************************************************/
    Vec1Int choseDataDistribution(const Int& nElements);

    /*******************************************************************************************
     * Record storing incar parameters.
     *******************************************************************************************/
    ShRec incar;
    /*******************************************************************************************
     * Record ML_AB data. List of local reference configurations will be overwritten during
     * execution of the Selector class.
     *******************************************************************************************/
    ShRec mlab;
    /*******************************************************************************************
     * MPI class to controll MPI communication.
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mlmpi;
    /*******************************************************************************************
     * Record storing data needed for Selector class, as the descriptors for structures in ML_AB.
     *******************************************************************************************/
    ShRec data;
    /*******************************************************************************************
     * Basis functions used for descriptor computation.
     *******************************************************************************************/
    BasisFunctionMap basisFunctions;
    /*******************************************************************************************
     * Choose between parallel modes OpenMP and MPI for farthest point sampling
     *******************************************************************************************/
    ParallelMode parallelMode;
    /*******************************************************************************************
     * Random seed supplied from the INCAR file as ML_RANDOM_SEED
     *******************************************************************************************/
    Vec1Int randomSeed;
    /*******************************************************************************************
     * When Selector class is used with MPI then descriptors are stored in a shmem array.
     * The vector is over all types in the supplied the ML_AB file.
     * The first dimnension of the shmem array are all central atoms of a certain type and
     * the second index of the shmem array is over all entries of the descriptor
     * super vector.
     *******************************************************************************************/
    std::vector<std::shared_ptr<ShmemArray2D<Real>>> descriptorsShmem;
    /*******************************************************************************************
     * Random number generator used in random selection of local reference configurations.
     *******************************************************************************************/
    UniqueRandomizer randomUnique;
    /*******************************************************************************************
     * Used to initialize farthest point sampling approaches.
     *******************************************************************************************/
    UniformRandomIntGenerator randomUniform;
    /*******************************************************************************************
     * Structurs extracted from ML_AB file. Length of vector macthes number of
     * structures in ml_AB file.
     *******************************************************************************************/
    std::vector<Structure> structures;
    /*******************************************************************************************
     * Number of local reference configurations after reslect was done.
     *******************************************************************************************/
    Vec1Int nLocConf;
    /*******************************************************************************************
     * Threshold for frathest point sampling with treshold distance.
     *******************************************************************************************/
    Vec1Real typeTresholds;
    /*******************************************************************************************
     * In flattened descriptor array per type for selection this index tells
     * from which structure a central atom was taken.
     *******************************************************************************************/
    Vec2Int strucIndex;
    /*******************************************************************************************
     * In flattened descriptor array per type for selection this index tells
     * the atom number of this central atom in ML_AB file.
     *******************************************************************************************/
    Vec2Int atomIndexInStruc;
    /*******************************************************************************************
     * Selected indexes for a certain type. Array is used in combination with strucIndex
     * and atomIndexInStruc to tell which central atoms have been selected.
     *******************************************************************************************/
    Vec1Int indexSelection;
    /*******************************************************************************************
     * Enum to choose selection method. Available methods are
     * fpsn farthest points sampling with n samples.
     * rann random selection with n samples.
     * fpst farthest point sampling with treshold criterion.
     * slhc SequentialLikelihoodCoverage
     *******************************************************************************************/
    SelectionMethod selectionMethod;
    /*******************************************************************************************
     * Enum to choose the metric which is used during the selection process. Available
     * metrics are
     * l1norm
     * l2norm
     * lpnorm
     * pkernel4 polynomial kernel norm. Only worls for polynomial power 4.
     * pcorr Pearson correlation coefficient as metric.
     *******************************************************************************************/
    SelectionMetric selectionMetric;
    /*******************************************************************************************
     * Descriptor which has as first dimension the central atom index flattened over all
     * structures in ML_AB file. Will be set up temporarily when selection process for a certain
     * atom type starts.
     *******************************************************************************************/
    Vec2Real typeOrderDesc;
    /*******************************************************************************************
     * Descriptor sizes for central atom flattened form. String key denotes the descriptor
     * index.
     *******************************************************************************************/
    std::map<String, Int> descSizes;
    /*******************************************************************************************
     * Descriptor offsets for central atom flattened form. String key denotes the descriptor
     * index.
     *******************************************************************************************/
    std::map<String, Int> descOffset;
    /*******************************************************************************************
     * List of active descriptors.
     *******************************************************************************************/
    Vec1String activeDescriptors;
    /*******************************************************************************************
     * Enum to chose the form of the data distribution. Which can be round robin
     * or block distributed.
     *******************************************************************************************/
    DataDist dataDist;
    /*******************************************************************************************
     * Stores central atom indexes for every structure in the ML_AB file which will
     * be manipulated on current MPI rank locally.
     *******************************************************************************************/
    Vec2Int centralIndexes;
    /*******************************************************************************************
     * Determines if a certain node has to write output.
     *******************************************************************************************/
    bool doWrite;
    /*******************************************************************************************
     * Checks if all descriptors are normalized. Faster Metric functions can be used
     * for normalized descriptors.
     *******************************************************************************************/
    bool descNormalized;
};

/*******************************************************************************************
 * @brief Checks if the number of MPI ranks exceeds the number of atoms and issues an error if so.
 *
 * This function ensures that the number of MPI ranks used does not exceed the number of atoms
 * in the machine learning atomic basis (ML_AB). If the number of MPI ranks is greater than the
 * minimum number of atoms, an error message is displayed, advising the user to reduce the number
 * of MPI ranks to avoid inefficiency and unnecessary environmental impact.
 *
 * @param numberAtoms A vector containing the number of atoms for each MPI rank.
 * @param mlmpi       An instance of the MlMPI class, which provides the number of MPI ranks.
 *
 * @note If the number of MPI ranks exceeds the number of atoms, the function will terminate
 *       execution with an error message.
 *******************************************************************************************/
void mpiRankCheck(const Vec1Int& numberAtoms, MlMPI& mlmpi);

ParallelMode get_parallelMode(const Record& incar);

} //namespace vaspml

#endif
