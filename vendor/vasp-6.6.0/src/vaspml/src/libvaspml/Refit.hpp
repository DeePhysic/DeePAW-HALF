#ifndef REFIT_HPP
#define REFIT_HPP

#include "BasisFunctions.hpp"
#include "BatchMap.hpp"
#include "FitMatrix.hpp"
#include "Frame.hpp"
#include "Inference.hpp"
#include "KernelPolynomial.hpp"
#include "Structure.hpp"
#include "types.hpp"

namespace vaspml
{

class MlMPI;
struct Record;

enum class FitAlgorithm
{
    DesignMatrix,
    NormalEquation
};

void fillRegCoefficients(Vec1Real& regCoeff, const Int n);

class Refit
{
  public:
    Refit(const ShRec& incar, const ShRec& mlab, std::shared_ptr<MlMPI> mpi = nullptr);
    /*******************************************************************************************
     * @brief Prepares the fitting process by initializing necessary data structures and
     *computations.
     *
     * @param structureData Reference to the structure data object containing information about the
     *structures.
     * @param mlffData Reference to the machine learning force field data handler.
     *******************************************************************************************/
    void prepareFit();
    /*******************************************************************************************
     * @brief Performs Singular Value Decomposition (SVD) fitting for regression analysis.
     *
     * This function uses SVD to solve a least squares problem for fitting regression coefficients
     * based on a design matrix and a set of energy values extracted from the dataset. The computed
     * regression coefficients are then written to a file for further use.
     *
     * The process involves:
     * - Extracting energy values from the dataset.
     * - Resizing the regression coefficients vector to match the number of reference
     *configurations.
     * - Solving the least squares problem using SVD.
     * - Writing the computed regression coefficients to a file.
     *
     * @note This function assumes that the design matrix and dataset are properly initialized
     *       and accessible through their respective member variables.
     *
     * @throws std::runtime_error If the SVD fitting fails or file writing encounters an error.
     *******************************************************************************************/
    void SVDDesignMatrixFitting();
    /*******************************************************************************************
     * @brief Performs SVD-based fitting of the normal equations to compute regression coefficients.
     *
     * This method solves the normal equations using Singular Value Decomposition (SVD)
     * to compute the regression coefficients. The routine retrieves the design matrix
     * and target vector, applies regularization, and writes the resulting coefficients
     * to a file.
     *
     * Steps:
     * 1. Retrieves the design matrix and target vector from the fitMatrix object.
     * 2. Initializes the regression coefficients vector.
     * 3. Computes the regularization parameter using input configuration values.
     * 4. Solves the least squares problem using SVD with regularization.
     * 5. Writes the computed regression coefficients to a file.
     *
     * @note The regression coefficients are written to "RegressionCoefficients.dat".
     *******************************************************************************************/
    void SVDNormalequationFitting();
    /*******************************************************************************************
     * @brief Computes the fit based on the specified fitting algorithm.
     *
     * This function determines the fitting algorithm to use and performs the fitting process.
     * Supported algorithms are:
     * - DesignMatrix: Uses the design matrix approach for fitting.
     * - NormalEquation: Uses the normal equation approach for fitting.
     *
     * @note If an unknown fitting algorithm is specified, an error message is logged.
     *******************************************************************************************/
    void computeFit();
    void SVDFittingTest();
    void dgemmTests();
    /*******************************************************************************************
     * @brief Writes the regression coefficients to a file.
     *
     * This function writes the regression coefficients stored in the `regressionCoefficients`
     * container to a file specified by the `filename` parameter. Each coefficient is formatted
     * as a scientific notation string with 24 characters width and 16 digits precision, and
     * written to a new line in the file.
     *
     * @param filename The name of the file where the regression coefficients will be written.
     *                 This should be a valid file path as a `String`.
     *
     * @note The function uses `file_io::openFileO` to open the file for writing. Ensure that
     *       the `file_io` namespace and its `openFileO` function are properly defined and
     *       accessible.
     *
     * @throws An exception may be thrown if the file cannot be opened or written to.
     *
     * @details Example usage:
     * @code
     * Refit refitObject;
     * refitObject.writeRegressionCoefficients("output.txt");
     * @endcode
     *
     * @see file_io::openFileO
     *******************************************************************************************/
    void writeRegressionCoefficients(const String& filename) const;
    /*******************************************************************************************
     * @brief Prepares the MLFF (Machine Learning Force Field) record by populating it with
     *necessary data.
     *
     * This function initializes and configures the MLFF record with various parameters,
     *descriptors, and metadata required for machine learning-based force field calculations. It
     *includes versioning, regression coefficients, RMSE placeholders, descriptor lists, reference
     *configurations, and additional MLAB-related data.
     *******************************************************************************************/
    void prepare_mlffRecord();
    /*******************************************************************************************
     * @brief Writes the MLFF (Machine Learning Force Field) data to a file.
     *
     * This function writes the MLFF data to the specified file if the MLFF has been prepared.
     * If the MLFF is not prepared, an error message is logged, and no file is written.
     *
     * @param filename The name of the file to which the MLFF data will be written.
     *
     * @note Ensure that `prepare_mlffRecord()` has been called before invoking this function.
     *       If the MLFF is not prepared, the function will log an error and take no action.
     *******************************************************************************************/
    void write_mlffFile(const String& filename);
    void trainDataInference();
    void testFittingDesignMatrix();
    void designMatrixTest();

  private:
    /*******************************************************************************************
     * @brief Computes batched frames for a set of structures.
     *
     * This function processes structures in batches and computes frames for each batch.
     * The number of frames (`nFrames`) and the batch size (`batchSize`) determine how the
     *structures are divided into batches. For each batch, the function computes a single frame
     *batch using the specified type order.
     *
     * @details
     * - The function retrieves the type order from `mlabData`.
     * - Iterates over the number of frames (`nFrames`) to process each batch.
     * - For each batch, calculates the start and end indices of the structures to be processed.
     * - Ensures the end index does not exceed the total number of structures (`nStrucs`).
     * - Calls `computeSingleFrame` on the `frames` object for the current batch,
     *   passing the range of structures and the type order.
     * - Outputs debug information about the current batch being processed.
     *
     * @pre `mlabData` must contain a valid `types` entry of type `Vec1String`.
     * @pre `frames` must be properly initialized and have enough capacity to store results for
     *all batches.
     * @pre `structures` must contain at least `nStrucs` elements.
     *
     *******************************************************************************************/
    void computeFrames();
    /*******************************************************************************************
     * @brief Sets the reference configurations for the Refit class.
     *
     * This function initializes and configures the reference configurations
     * used in the Refit class. It retrieves necessary data from `mlabData`,
     * sets the number of structures in `strucToBatch`, and creates a shared
     * pointer to a `DescriptorCollector` containing the reference configurations.
     * Finally, it updates the `kernelPolynomial` with the created reference configurations.
     *
     * @details
     * - The function retrieves the `lrcStructure` and `lrcAtom` data from `mlabData`.
     * - It sets the number of structures in `strucToBatch` using the `numStructures` value from
     *`mlabData`.
     * - A `DescriptorCollector` for reference configurations is created using the
     *   `batch_tools::createDescriptorCollectorReferenceConfigs` function.
     * - The created reference configurations are passed to `kernelPolynomial` using the
     *`set_descriptorsRefConfs` method.
     *
     * @note This function assumes that `mlabData`, `batchFrames`, `atomBatchMap`, `strucToBatch`,
     * and `kernelPolynomial` are properly initialized before calling this function.
     *******************************************************************************************/
    void setReferenceConfigurations();
    /*******************************************************************************************
     * @brief Computes the design matrix for a batch of frames.
     *
     * This function iterates over a batch of frames and updates the kernel polynomial
     * for each frame in the batch. It prints the progress to the console for each frame
     * being processed. Optionally, the kernel matrix can be written to a file (commented out).
     *
     * @details
     * - The function assumes that `batchFrames` is a container holding the frames to be processed.
     * - The `kernelPolynomial` object is used to update the kernel for each frame.
     * - The commented-out line can be used to write the kernel matrix to a file for debugging or
     *analysis purposes.
     *
     * @note Ensure that `batchFrames` and `kernelPolynomial` are properly initialized before
     *calling this function.
     *
     * @warning This function uses `std::cout` for logging, which may impact performance for large
     *batches.
     *
     * Example usage:
     * @code
     * Refit refitObject;
     * refitObject.computeDesignMatrixBatch();
     * @endcode
     *******************************************************************************************/
    void computeFitMatrix();
    /*******************************************************************************************
     * @brief Prepares the number of descriptors for each key in the batch frames.
     *
     * This function initializes the `numberDescriptors` map, which stores the number of descriptors
     * for each key. The number of descriptors is determined based on the number of basis functions
     * and the total number of types. It assumes that all types have the same number of descriptors,
     * which is true when setting up a new refit.
     *
     * @details
     * - The descriptors are retrieved from the first frame in the `batchFrames`.
     * - The maximum number of types (`nTypes`) is fetched from the `mlabData` object.
     * - For each descriptor in the `descriptors` map, the size of the descriptor is obtained
     *   using the `get_sizeDescriptor` method.
     * - A shared pointer to a vector (`Vec1Int`) is created for each key, where the vector
     *   contains `nTypes` elements, each initialized to the size of the descriptor (`nDesc`).
     * - This function assumes that all types have the same number of descriptors, which is
     *   always true during the setup of a new refit.
     *
     * @note This function operates on the assumption that the descriptors are consistent across
     *       all batch frames and that the number of descriptors is uniform across all types.
     *
     * @pre `batchFrames` must contain at least one frame.
     * @pre `mlabData` must be properly initialized and contain the "maxTypes" parameter.
     *
     * @post The `numberDescriptors` map is populated with shared pointers to vectors containing
     *       the number of descriptors for each key.
     *******************************************************************************************/
    void prepare_numberDescriptors();
    /*******************************************************************************************
     * @brief Prepares the kernel polynomial by setting up reference configurations.
     *
     * This function initializes the kernel polynomial by calling the `setReferenceConfigurations`
     * method, which sets up the necessary reference configurations for the kernel polynomial.
     *******************************************************************************************/
    void prepare_kernelPolynomial();

    /*******************************************************************************************
     * Data which is read from the ML_AB file created during on the fly training.
     *******************************************************************************************/
    ShRec mlab;
    /*******************************************************************************************
     * Force field parameters which are read from the INCAR file.
     *******************************************************************************************/
    ShRec                  incar;
    std::shared_ptr<MlMPI> mpi;
    Inference              inference;
    /*******************************************************************************************
     * Stores generated mlff generated during refit.
     *******************************************************************************************/
    ShRec mlff;
    /*******************************************************************************************
     * Stores basis functions in which local refrence configurations are computed
     *******************************************************************************************/
    BasisFunctionMap basisFunctions;
    /*******************************************************************************************
     * Batch Size. Set to value larger than 1 to use batching. Currently there are no fitting
     * routines that support batching.
     *******************************************************************************************/
    Size batchSize = 1;
    /*******************************************************************************************
     * Number of Frame classes. This is the same as the number of batches.
     *******************************************************************************************/
    Int nFrames;
    /*******************************************************************************************
     * Number of descriptors for each key in descriptor map.
     *******************************************************************************************/
    std::map<String, ShVec1Int> numberDescriptors;
    /*******************************************************************************************
     * Structures extracted from the ML_ABdata. Are needed for Frame computation.
     *******************************************************************************************/
    std::vector<Structure> structures;
    /*******************************************************************************************
     * Stores frame classes containing local descriptors. If batchSize is larger than 1, the
     * frames will contain local reference configurations for each batch.
     *******************************************************************************************/
    std::vector<Frame> frames;
    /*******************************************************************************************
     * Number of structures in dataset which are used for fitting.
     *******************************************************************************************/
    Size nStrucs;
    /*******************************************************************************************
     * Total number of reference configurations used for fitting.
     *******************************************************************************************/
    Int nRefConfs;
    /*******************************************************************************************
     * NUmber of local reference configurations per atom type in the dataset.
     *******************************************************************************************/
    Vec1Int nLocRefConfs;
    /*******************************************************************************************
     * Maps structures to batches. This is used to compute the local reference configurations.
     *******************************************************************************************/
    std::vector<AtomBatchMap> atomBatchMap;
    /*******************************************************************************************
     * Mapping between global structure number to local batch number.
     *******************************************************************************************/
    BatchManager strucToBatch;
    /*******************************************************************************************
     * polynomial kernel used for setting up the design matrix.
     *******************************************************************************************/
    KernelPolynomial kernelPolynomial;
    /*******************************************************************************************
     * Design matrix used for regression fitting.
     *******************************************************************************************/
    FitMatrix fitMatrix;
    /*******************************************************************************************
     * @brief Stores the regression coefficients for the kernel polynomial regression model.
     *
     * This variable holds the coefficients used in a regression model, where each coefficient
     * represents the weight associated with a corresponding feature or predictor variable.
     *
     * @details
     * - The coefficients are typically computed during the training phase of the model.
     * - These values are used to make predictions by applying the regression equation.
     *******************************************************************************************/
    Vec1Real regressionCoefficients;
    /*******************************************************************************************
     * Flag to check if mlff is prepared
     *******************************************************************************************/
    bool         mlffPrepared = false;
    Real         energyShift;
    FitAlgorithm fitAlgorithm;
};

std::vector<Frame> makeBatchFrameStorage(const Int&              nFrames,
                                         const Record&           incar,
                                         const Record&           mlab,
                                         const BasisFunctionMap& basisFunctions);

std::vector<AtomBatchMap> initAtomBatchMaps(const std::vector<Structure>& structures,
                                            const Size                    nFrames,
                                            const Size                    batchSize);

} //namespace vaspml

#endif
