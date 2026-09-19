#ifndef KERNELPOLYNOMIAL_HPP
#define KERNELPOLYNOMIAL_HPP

#include "ArrayResizing.hpp"
#include "Kernel.hpp"
#include "types.hpp"

#include <tuple>  // for tuple

namespace vaspml
{

class Frame;
class MlMPI;
class TypeMap;
struct Record;

struct KernelPolynomialLookUp
{
    /*******************************************************************************************
     * Index running trough all combinations of centralAtom, type in structure and kernel entries
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * stores central atom index for all combinations of centralAtom, type in structure and kernel
     * entries
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * store types indices for all combinations of centralAtom, type in structure and
     * kernel entries
     *******************************************************************************************/
    Vec1Size typeStruc;
    /*******************************************************************************************
     * stores type indices in force field combinations of centralAtom, type in structure and
     * kernel entries
     *******************************************************************************************/
    Vec1Size typeFF;
    /*******************************************************************************************
     * Store kernel Index entries for all combinations of combinations of centralAtom, type in
     * structure and kernel entries
     *******************************************************************************************/
    Vec1Size kernelEntry;
    /*******************************************************************************************
     * computes the needed sizes of the flattened arrays
     *
     * @param nAtoms number of atoms
     * @param typeCentral type index of the central atom in the structure
     * @param typeMap mapping between force field types and structure types
     * @param nLocRefConfTypes number of local reference configurations per atom types
     *******************************************************************************************/
    Size countElements(const Size&    nAtoms,
                       const Vec1Int& typeCentral,
                       const TypeMap& typeMap,
                       const Vec1Int& nLocRefConfTypes);
    /*******************************************************************************************
     * Routine checks if a resize and a refill of the LookUp table is needed
     *
     * @param nAtoms number of atoms
     * @param typeCentral type index of central atom in structure
     * @param typeMap conversion map between force field type indices and structure type indices
     * @param nLocRefConfTypes number of local reference configurations per atom type in force field
     * @param centralAtomIndexPerType central atom index within it's type
     * @param kernelPolynomialLookUpSize structure to check if a resize and refill is needed
     *******************************************************************************************/
    void update(const Size&      nAtoms,
                const Vec1Int&   typeCentral,
                const TypeMap&   typeMap,
                const Vec1Int&   nLocRefConfTypes,
                const Vec1Int&   centralAtomIndexPerType,
                ArrayResizing1D& kernelPolynomialLookUpSize);
    /*******************************************************************************************
     * @param nAtoms number of atoms in current structure
     * @param typeCentral type index of central atom
     * @param typeMap mapping between type index in structure and force field and vice versa
     * @param nLocRefConfTypes number of local reference configurations per atom type in force field
     * @param centralAtomIndexPerType central atom index within its atom type
     *******************************************************************************************/
    void refill(const Size&    nAtoms,
                const Vec1Int& typeCentral,
                const TypeMap& typeMap,
                const Vec1Int& nLocRefConfTypes,
                const Vec1Int& centralAtomIndexPerType);
};

/*******************************************************************************************
 * @class KernelPolynomial
 * @brief Inherits from Kernel class and implements PolynomialKernel and needed derivatives
 *
 * This class inherits from Kernel and implements the functionality for a polynomial
 * Kernel. The polynomial kernel needs to implement the derivative terms related purely to
 * the polynomial part of the kernel. These derivative terms are important in force
 * computations. The polynomial kernel will be computed from a Frame class instance
 * which supplies the needed descriptors as DescriptorSHS2, DescriptorSHS3 or DescriptorSHS3Lin.
 * These are used and combined to a DescriptorCollector class from which the kernel
 * with respect to some reference configurations supplied from the ML_FF Record.
 * The KernelPolynomial can be used in the Predictor class to predict energy, forces and
 * stress tensor
 *******************************************************************************************/
class KernelPolynomial : public Kernel
{
  public:
    /*******************************************************************************************
     * Empty constructor of KernelPolynomial class.
     *
     * @note class is not usable in this state, but can be transformed into a usable state by
     * calling KernelPolynomial::init().
     *******************************************************************************************/
    KernelPolynomial(ShRec record = nullptr);
    /*******************************************************************************************
     * @brief Constructs a KernelPolynomial object for polynomial kernel computations.
     *
     * This constructor initializes the KernelPolynomial object with the provided input parameters,
     * MPI environment, execution policy, and shared record. It also sets up the polynomial kernel
     * properties and retrieves the kernel power from the input parameters.
     *
     * @param inputParameters A Record object containing the input parameters for the kernel.
     * @param mpiIn A shared pointer to an MlMPI object for managing MPI communication.
     * @param record A shared record (ShRec) for managing kernel-related data.
     *
     * @note The constructor initializes the following:
     * - `polynomialDerivative`: A vector of real numbers representing the polynomial derivative.
     * - `polynomialKernelNorm`: A vector of real numbers representing the polynomial kernel
     *normalization.
     * - `polynomialKernel`: The kernel matrix used for polynomial kernel computations.
     * - `kernelPower`: The power of the polynomial kernel, retrieved from the input parameters
     *(`ML_NHYP`).
     *******************************************************************************************/
    KernelPolynomial(const Record&                 inputParameters,
                     const std::shared_ptr<MlMPI>& mpiIn = nullptr,
                     ShRec                         record = nullptr);
    /*******************************************************************************************
     * @brief Constructs a KernelPolynomial object for machine learning kernel operations.
     *
     * This constructor initializes the KernelPolynomial object with the provided input parameters
     * and sets up the polynomial kernel and its associated properties.
     *
     * @param incar          Input record containing configuration parameters.
     * @param mlab           Record containing machine learning-related data.
     * @param mpiIn          Shared pointer to the MPI wrapper for parallel execution.
     * @param record         Shared record object, either provided or created internally.
     *
     * The constructor performs the following:
     * - Initializes the base Kernel class with the provided parameters and a record.
     * - Retrieves and stores the polynomial derivative and kernel normalization data.
     * - Initializes the polynomial kernel using the kernel matrix.
     * - Reads the kernel power (ML_NHYP) from the input configuration.
     *******************************************************************************************/
    KernelPolynomial(const Record&                 incar,
                     const Record&                 mlab,
                     const std::shared_ptr<MlMPI>& mpiIn = nullptr,
                     ShRec                         record = nullptr);
    /*******************************************************************************************
     * @brief Initializes the KernelPolynomial object with input parameters, MPI instance, and
     *execution policy.
     *
     * This method sets up the KernelPolynomial by calling the base Kernel class's initialization
     *method and retrieves the kernel power parameter from the input configuration.
     *
     * @param inputParameters A Record object containing the input parameters for initialization.
     * @param mpiIn A shared pointer to an MlMPI instance for managing MPI-related operations.
     *
     * @note The kernel power parameter is fetched using the key "ML_NHYP" from the input
     *parameters.
     *******************************************************************************************/
    void init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * update the polynomial Kernel
     *
     * compute kernelMatrix in the Parent class compute_kernelMatrix
     * allocate or resize the needed working arrays
     * compute the polynomialKernel and terms needed for derivatives
     *
     * @param frame contains the description of the structure to compute kernel of
     *******************************************************************************************/
    void updatePolynomialKernel(const Frame& frame);
    /*******************************************************************************************
     * @brief Updates the polynomial kernel with refit adjustments.
     *
     * This function updates the kernel using the provided frame, allocates
     * necessary derivative arrays for refitting, and computes the polynomial
     * kernel terms based on the frame's type map.
     *
     * @param frame The input frame containing data required for kernel updates.
     *******************************************************************************************/
    void updatePolynomialKernelRefit(const Frame& frame);
    /*******************************************************************************************
     * computing terms polynomialDerivative, polynomialKernelNorm, polynomialKernel
     *******************************************************************************************/
    void computePolynomialKernelTerms(const TypeMap& typeMap);
    /*******************************************************************************************
     * @brief Computes and refits the polynomial kernel terms based on the provided type map.
     *
     * This function calculates the polynomial kernel terms and their derivatives, refitting them
     * as necessary. The computation is executed based on the specified execution policy, which
     * can either be single-core CPU or GPU with standard library support.
     *
     * @param typeMap A mapping of types used for kernel computation.
     *
     * Execution Details:
     * - If the execution policy is `ExecutionPolicy::cpuSingleCore`, the data is prepared in a
     *   standard algorithm-compatible format and processed on the CPU.
     * - If the execution policy is `ExecutionPolicy::gpuStdLib`, the computation is performed
     *   on the GPU using parallel execution.
     *******************************************************************************************/
    void computePolynomialKernelTermsRefit(const TypeMap& typeMap);
    /*******************************************************************************************
     * Computing the derivative of the kernel with respect to descriptor vector.
     *
     * @param derivativeMatrix stores the derivative of the kernel
     * @f$K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=
     * (\hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B})^{\xi}@f$ with respect to
     * the descriptor vector @f$\hat{\mathbf{X}}_{i}@f$
     * @param forceVector is a temporary array used to compute the derivativeMatrix which
     * has to be allocated with proper size from call routine. Size is number of atoms per type
     * @f$-\sqrt{w}\Lambda'_{iB}\mathbf{w}_{B}@f$
     * @param descriptorsRefConfsWeighted reference configurations weighted with regression
     coefficients
     * @f$\hat{\mathbf{X}}_{B}^{'}=\mathbf{w}_{B}\circ \hat{\mathbf{X}}_{B}@f$
     * @param typeStruc index defining type of analysed atoms in structure
     * @param typeForceField index defining type of analysed atoms in force field
     *
     * @f[
       \mathbf{L}^{i} = \sum_{B} \Lambda_{iB} \hat{\mathbf{X}}_{B}' - \sum_{B} w_{B}
                        \Lambda'_{iB}\hat{\mathbf{X}}_{i}
       @f]
       with
       @f[\hat{\mathbf{X}}_{B}' = w_{B} \hat{\mathbf{X}}_{B}@f]
       @f[\Lambda_{iB} = \frac{\xi}{||\mathbf{X}_{i}||}
          \left(\hat{\mathbf{X}}_{i} \cdot \hat{\mathbf{X}}_{B} \right)^{\xi-1} @f]
       @f[\Lambda'_{iB}= \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i}
          \cdot \hat{\mathbf{X}}_{B} \right)^{\xi}@f]
     *******************************************************************************************/
    void computeDerivativeKernel(Vec1Real&     derivativeMatrix,
                                 Vec1Real&     forceVector,
                                 const Real*   descriptorsRefConfsWeighted,
                                 const Real*   regressionCoefficients,
                                 const Size&   typeStruc,
                                 const Size&   typeForceField,
                                 const String& key) const;
    /*******************************************************************************************
     * Seting up lookUp tables needed for the c++ std parallel algorithms
     *
     * @param typeMap mapping between force field types and structure types
     *
     * check KernelPolynomialLookUp struct for more information
     *******************************************************************************************/
    void makeLookUpTables(const TypeMap& typeMap);
    /*******************************************************************************************
     *  write polynomialDerivative to a file in a single column format (flattened)
     *
     *  @param fname file name to which the matrix will be written
     *******************************************************************************************/
    void write_polynomialDerivative(const String& fname) const;
    /*******************************************************************************************
     * write polynomialKernelNorm in single column format to file
     *
     * @param fname determines file name to which matrix is written
     *******************************************************************************************/
    void write_polynomialKernelNorm(const String& fname) const;
    /*******************************************************************************************
     * writes polynomialKernelNorm to file in single column format
     *
     * @param fname file name to which flattened matrix will be written
     *******************************************************************************************/
    void write_polynomialKernel(const String& fname) const;
    /*******************************************************************************************
     * @brief Retrieves the full polynomial derivative.
     *
     * @return A constant reference to the 2D vector (Vec2Real) representing the polynomial
     *derivative.
     *******************************************************************************************/
    const Vec2Real& get_polynomialDerivative() const;
    /*******************************************************************************************
     * @brief Retrieves the polynomial derivative for a specific type index.
     *
     * @param typeIndx The index (Int) specifying the type of polynomial derivative to retrieve.
     * @return A constant reference to the 1D vector (Vec1Real) corresponding to the specified type
     *index.
     *******************************************************************************************/
    const Vec1Real& get_polynomialDerivative(const Int& typeIndx) const;

  private:
    /*******************************************************************************************
     * Actual calculations in computing terms polynomialDerivative, polynomialKernelNorm,
     *polynomialKernel
     *
     * Single cpu version.
     *
     * @param typeMap Map from types of structure to types in force field
     * @param polynomialDerivative derivative of polynomial kernel
     *******************************************************************************************/
    void computePolynomialKernelTermsCPU(
        const TypeMap&                                            typeMap,
        std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data);
    /*******************************************************************************************
     * @brief Computes polynomial kernel terms for a given dataset on the CPU, refitting the kernel
     *values.
     *
     * This function calculates the polynomial kernel terms for a dataset by iterating over atom
     *types and their corresponding local reference configurations. The kernel values are updated
     *in-place based on the specified kernel power and normalization factors.
     *
     * @param typeMap A mapping of structure types to their corresponding indices and properties.
     * @param data A nested vector of tuples containing kernel values. Each tuple consists of:
     *             - `Real&`: The kernel value to be updated.
     *             - `Real&`: The derivative of the kernel value.
     *             - `Real&`: An auxiliary value (unused in this function).
     *
     * @note The function assumes that the `descriptorCollection` and `numberLocalRefConfs` members
     *are properly initialized before calling this function. It also relies on the
     *`math::intPowMixAlgo` utility for efficient integer power computation.
     *******************************************************************************************/
    void computePolynomialKernelTermsRefitCPU(
        const TypeMap&                                     typeMap,
        std::vector<std::vector<std::tuple<Real&, Real&>>> data);

    //    /*******************************************************************************************
    //     * Actual calculations in computing terms polynomialDerivative, polynomialKernelNorm,
    //     *polynomialKernel
    //     *
    //     *
    //     * @param typeMap Map from types of structure to types in force field
    //     * @param data is a flattened version of { polynomialKernel, polynomialDerivative and
    //     * polynomialKernelNorm.
    //     *******************************************************************************************/
    //    //void computePolynomialKernelTermsGPU(const TypeMap& typeMap);

    /*******************************************************************************************
     * computes the polynomial Kernel terms by the help of look up tables and a single flat main
     *loop
     *
     * computes polynomialKernel, polynomialDerivative, polynomialKernelNorm
     *******************************************************************************************/
    void computePolynomialKernelTerms_stdlib();
    /*******************************************************************************************
     * @brief Computes and refits the polynomial kernel terms for a given descriptor collection.
     *
     * This function recalculates the polynomial kernel terms and their derivatives based on the
     * provided kernel power and descriptor data. It uses a parallelized loop to efficiently
     * process the kernel entries and updates the polynomial kernel and its derivative values.
     *
     * @details
     * - The function iterates over the indices of the kernel polynomial lookup table.
     * - For each index, it computes the polynomial kernel value and its derivative using
     *   integer power algorithms.
     * - The results are stored in the `polynomialKernel` and `polynomialDerivative` data
     *structures.
     * - The computation is parallelized using `VASPML_PAR_UNSEQ` for improved performance.
     *
     * @note
     * - The function assumes that the `descriptorCollection`, `kernelPower`, `polynomialKernel`,
     *   and `polynomialDerivative` are properly initialized before calling this function.
     * - The `math::intPowMixAlgo` function is used for efficient integer power calculations.
     *
     * @param None
     * @return None
     *******************************************************************************************/
    void computePolynomialKernelTermsRefit_stdlib();
    /*******************************************************************************************
     * allocating the arrays where the terms needed for computing the kernel derivatives in
     * computeDerivativeKernel are allocated.
     *******************************************************************************************/
    void allocateDerivativeArrays();
    /*******************************************************************************************
     * @brief Allocates and initializes derivative arrays for kernel polynomial refitting.
     *
     * This function handles the allocation and initialization of the `polynomialDerivative`
     * and `polynomialKernel` arrays based on the current execution policy.
     *
     * - For `ExecutionPolicy::cpuSingleCore`, it resizes the `polynomialDerivative` array
     *   to match the dimensions of the `kernelMatrix` and copies the `kernelMatrix` to
     *   `polynomialKernel`.
     * - For `ExecutionPolicy::gpuStdLib`, it delegates the allocation to the GPU-specific
     *   implementation via `allocateDerivativeArraysRefitGPU()`.
     *******************************************************************************************/
    void allocateDerivativeArraysRefit();
    /*******************************************************************************************
     * allocating the arrays where the terms needed for computing the kernel derivatives in
     * computeDerivativeKernel are allocated.
     *
     * routine uses checks to avoid memory copies of the arrays from host to device and vice versa
     *******************************************************************************************/
    void allocateDerivativeArraysGPU();
    /*******************************************************************************************
     * @brief Allocates and resizes derivative arrays for GPU refitting in the Kernel Polynomial.
     *
     * This function ensures that the necessary arrays for polynomial derivatives and kernel
     *matrices are properly allocated and resized based on the number of atoms and kernel matrix
     *dimensions. It checks the current sizes of the arrays and resizes them if needed to
     *accommodate the required dimensions.
     *
     * Steps:
     * - Retrieves the number of atoms from the descriptor collection.
     * - Checks if resizing is required for the central atom index size.
     * - Checks and resizes the polynomial derivative array for 1D or 2D dimensions as needed.
     * - Updates the `polynomialKernel` to reference the `kernelMatrix` after resizing.
     *
     * @note This function is GPU-specific and assumes that the descriptor collection and kernel
     *matrix dimensions are already initialized.
     *******************************************************************************************/
    void allocateDerivativeArraysRefitGPU();
    /*******************************************************************************************
     * polynomial power @f$\xi@f$ which is taken of the kernelMatrix
     * @f$K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=\left (
     *        \hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B}\right )^{\xi}@f$
     *******************************************************************************************/
    Int kernelPower;
    /*******************************************************************************************
     * derivative of polynomial part of kernel function
     *
     * @f[ \frac{\xi}{\lVert X_{i}\rVert}K^{\xi-1}_{ij} @f]
     *******************************************************************************************/
    Vec2Real& polynomialDerivative;
    /*******************************************************************************************
     * Structure to check if a reszize of the polynomialDerivative is needed
     *******************************************************************************************/
    ArrayResizing2D polynomialDerivativeSize;
    /*******************************************************************************************
     * polynomial kernel times norm factor
     *
     * @f[ \frac{\xi}{\lVert X_{i}\rVert}K^{\xi}_{ij} @f]
     *******************************************************************************************/
    Vec2Real& polynomialKernelNorm;
    /*******************************************************************************************
     * stores the polynomial power of the kernelMatrix
     *
     * points to the memory location of the parent class
     * where kernelMatrix is stored and overwrites it with the polynomial kernel
       @f[K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=
       \left (\hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B}\right )^{\xi}@f]
     *******************************************************************************************/
    Vec2Real& polynomialKernel;
    /*******************************************************************************************
     * Object to check if the number of central atoms has changed.
     *******************************************************************************************/
    ArrayResizing1D centralAtomIndexSize;
    /*******************************************************************************************
     * stores look up tables to compute the polynomial kernel terms in flattend loop
     *******************************************************************************************/
    KernelPolynomialLookUp kernelPolynomialLookUp;
    /*******************************************************************************************
     * Structure to check if kernelPolynomialLookUp needs a recomputation and resize
     *******************************************************************************************/
    ArrayResizing1D kernelPolynomialLookUpSize;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapKernelPolynomial();

} //namespace vaspml

#endif
