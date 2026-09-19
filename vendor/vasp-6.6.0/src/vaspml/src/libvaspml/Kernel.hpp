#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "ArrayResizing.hpp"
#include "types.hpp"

namespace vaspml
{

class Descriptor;
class DescriptorCollector;
class Frame;
class MlMPI;
class TypeMap;
struct Record;
template<class T>
class ShmemArray2DVariableLen;

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> generate_ShmemArrayMap(
    const std::map<String, Real>&      weights,
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::shared_ptr<MlMPI>&      mpiIn);

/*******************************************************************************************
 * @class Kernel
 * @brief this class is a purely virtual base class for any Kernel implementation
 *
 * The class implements the minimum functions which are needed for defining a Kernel matrix
 * expressed in terms of descriptors which are derived classes of the Descriptor virtual
 * base class.
 *******************************************************************************************/
class Kernel
{
  public:
    /*******************************************************************************************
     * Empty constructor of Kernel class.
     *
     * @note Class is not usable in this state, but can be transformed into a usable state by
     * calling Kernel::init().
     *******************************************************************************************/
    Kernel(ShRec record = nullptr);
    /*******************************************************************************************
     * @brief Constructs a Kernel object with the given parameters.
     *
     * This constructor initializes the Kernel object by setting up various data structures
     * and configurations required for kernel-based machine learning computations. It handles
     * the initialization of shared memory arrays, MPI communication, and other internal
     * components necessary for the kernel's operation.
     *
     * @param inputParameters A Record object containing input parameters for the kernel setup.
     * @param mpiIn A shared pointer to an MlMPI object for MPI communication.
     * @param record A ShRec object used to initialize the kernel's data.
     *******************************************************************************************/
    Kernel(const Record&                 inputParameters,
           const std::shared_ptr<MlMPI>& mpiIn = nullptr,
           ShRec                         record = nullptr);
    /*******************************************************************************************
     * @brief Constructs a Kernel object with the specified configuration and dependencies.
     *
     * @param incar Input record containing configuration parameters.
     * @param mlab Input record containing machine learning-related parameters.
     * @param mpiIn Shared pointer to the MlMPI object for MPI communication.
     * @param record Shared record object for data initialization.
     *
     * This constructor initializes the Kernel object by:
     * - Assigning or creating a record using `assignOrMakeRecord`.
     * - Initializing MPI communication via `mpiIn`.
     * - Setting up weights for 2-body and 3-body interactions using parameters from `incar`.
     * - Loading the number of local reference configurations from `mlab`.
     * - Initializing a `DescriptorCollector` for managing descriptors.
     * - Accessing the kernel matrix from the data record.
     * - Storing the execution policy for later use.
     *******************************************************************************************/
    Kernel(const Record&                 incar,
           const Record&                 mlab,
           const std::shared_ptr<MlMPI>& mpiIn = nullptr,
           ShRec                         record = nullptr);
    /*******************************************************************************************
     * Destructor, required because this is an abstract base class.
     *******************************************************************************************/
    virtual ~Kernel() = default;
    /*******************************************************************************************
     * @brief Initializes the Kernel object with input parameters, MPI instance, and execution
     *policy.
     *
     * This method sets up the internal state of the Kernel object, including weights, MPI instance,
     * local reference configurations, descriptor maps, and the descriptor collection. It also
     *prepares the descriptor data by invoking the `fillDescriptor` method.
     *
     * @param inputParameters A reference to the input parameters containing configuration data.
     * @param mpiIn A shared pointer to the MPI instance used for parallel execution.
     *******************************************************************************************/
    void init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * updating the Kernel class by computing the kernelMatrix
     *
     * @param frame describes current atomic structure with help of descriptors DescriptorSHS2,
     * DescriptorSHS3
     *
     * routine performs the following steps:
     * |step | action                                                |
     * |:----|:------------------------------------------------------|
     * | 1   | resize kernelMatrix if necessary                      |
     * | 2   | set kernelMatrix to zero                              |
     * | 3   | compute descriptorCollector from supplied frame class |
     * | 4   | compute kernelMatrix                                  |
     *
     *******************************************************************************************/
    void updateKernel(const Frame& frame);
    /*******************************************************************************************
     * Get local reference configurations described in descriptor space of 2-body (DescriptorSHS2),
     * 3body (DescriptorSHS3),... descriptors.
     *
     * Descriptors are described in terms of DescriptorSHS2, DescriptorSHS3...
     * map can chose from SHS2-2-body SHS3-3-body; map strings are defined in constants.hpp
     *******************************************************************************************/
    const std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>>&
    get_descriptorsRefConfs() const;
    /*******************************************************************************************
     * getter for number of local reference configurations
     *******************************************************************************************/
    ShVec1Int get_numberLocalRefConfs() const;
    /*******************************************************************************************
     * getter for number of Descriptors ( basis set size )
     *******************************************************************************************/
    const std::map<String, ShVec1Int>& get_numberDescriptors() const;
    /*******************************************************************************************
     * getting the feature space size which is the product of number
     * descriptors times local reference configurations
     *******************************************************************************************/
    const std::map<String, ShVec1Int>& get_featureSpaceSize() const;
    /*******************************************************************************************
     * getter for descriptor weights ( not regression coefficients )
     *
     * map entries belong to SHS2-2-body (DescriptorSHS2), SHS3-3-body (DescriptorSHS3)
     *******************************************************************************************/
    const std::map<String, Real>& get_weights() const;
    /*******************************************************************************************
     * getter for kernelMatrix
     *
     * @note kernelMatrix will be overwritten by the kernel function defined in the child
     * classes, as for example KernelPolynomial
     *
     * dimensions of Kernel matrix are atom types as first dimension and, and a combined index
     * of atoms per type times number local reference configurations
     * @f$\mathbf{K}_{iB}=\mathbb{R}^{N_{types} \times (N_{atoms\ per\ type}\times N_{local\ ref
     *\ confs})}@f$
     *******************************************************************************************/
    const Vec2Real& get_kernelMatrix() const;
    /*******************************************************************************************
     * get the number of atoms per type belonging to structure from which kernelMatrix was computed
     *******************************************************************************************/
    const Vec1Int& get_nAtomsType() const;
    /*******************************************************************************************
     * return type map which maps force field types to structure types and vice versa
     *******************************************************************************************/
    std::shared_ptr<const TypeMap> get_typeMap() const;
    /*******************************************************************************************
     * give raw pointer to where kernelMatrix for certain atom of certain type starts
     *
     * @param type atom type for which pointer is retrieved
     * @param atomInType atom number of certain type for which pointer is retrieved
     *
     * @note the length of the slice for a certain atom of a type is determined by
     * numberLocalRefConfs[ typeMap -> toType( type ) ]
     *******************************************************************************************/
    const Real* get_kernelMatrixAtom(Size type, Size atomInType) const;
    /*******************************************************************************************
     * Computing derivative of kernel.
     *
     * formally this can be written as
     * @f[\sum_{B}w_{B}\frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}
     * {\partial \hat{\mathbf{X}}_{i} }@f]
     * for details about implementation please check the derived classes as for
     * example KernelPolynomial
     *******************************************************************************************/
    virtual void computeDerivativeKernel(Vec1Real&     derivativeMatrix,
                                         Vec1Real&     forceVector,
                                         const Real*   descriptorsRefConfsWeighted,
                                         const Real*   regressionCoefficients,
                                         const Size&   typeStruc,
                                         const Size&   typeForceField,
                                         const String& key) const = 0;
    /*******************************************************************************************
     * getter for descriptorCollection stores normalized descriptors and neighbor Lists
     *******************************************************************************************/
    const DescriptorCollector& get_descriptorCollection() const;
    /*******************************************************************************************
     * Get pointer to descriptorCollection
     *******************************************************************************************/
    const std::shared_ptr<DescriptorCollector>& get_descriptorCollectionPtr() const;
    /*******************************************************************************************
     * Write kernel matrix to file in a flattened 1D vector form.
     *******************************************************************************************/
    void write_kernelMatrix(const String& fname) const;
    /*******************************************************************************************
     *******************************************************************************************/
    void set_descriptorsRefConfs(
        const std::shared_ptr<DescriptorCollector>& referenceConfigurations);
    /*******************************************************************************************
     * @brief Sets the number descriptors for the kernel.
     *
     * This method assigns a map of number descriptors to the kernel.
     * Each descriptor is represented as a key-value pair, where the key
     * is a string identifier and the value is a shared vector of integers.
     *
     * @param numberDescriptors A map containing the number descriptors to set.
     *******************************************************************************************/
    void set_numberDescriptors(const std::map<String, ShVec1Int>& numberDescriptors);
    /*******************************************************************************************
     * @brief Generates the feature space map by computing the element-wise product of local
     *reference configurations and the number of descriptors for each descriptor key.
     *
     * This function iterates over a predefined list of descriptor keys and calculates the feature
     *space size for each key. The result is stored in the `featureSpaceSize` map, where each key is
     *associated with a shared pointer to a vector of integers representing the computed feature
     *space size.
     *
     * @details
     * - The local reference configurations (`numberLocalRefConfs`) are used as a base.
     * - The element-wise product is computed with the corresponding descriptor counts
     *(`numberDescriptors`).
     * - The descriptor keys are defined in `constants::descriptorKeyList`.
     *
     * @note Assumes that `numberLocalRefConfs` and `numberDescriptors` are properly initialized.
     *******************************************************************************************/
    void generate_featureSpaceMap();
    /*******************************************************************************************
     * @brief Creates a shared memory array map for descriptors and reference configurations.
     *
     * This function initializes the `descriptorsRefConfs` member by generating a shared memory
     * array map using the provided weights, feature space size, and MPI context.
     *
     * @note The generated shared memory array map is used to efficiently manage data
     *       across multiple processes in a parallel computing environment.
     *******************************************************************************************/
    void make_ShmemArrayMap();

  private:
    /*******************************************************************************************
     * setting up the descriptor collector class;
     *
     * This class normalizes the supplied descriptors and stores them in the proper format for the
     * Kernel class. The format is [type][ number atoms per type x number of descriptors for type ]
     *******************************************************************************************/
    void setUpDescriptorCollector(
        const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors);
    /*******************************************************************************************
     * compute the Kernel matrix for the supplied descriptors stored in descriptor collection
     * Kernel matrix is defined by a kernel function @f$ K(:,\hat{\mathbf{X}}_{B}) @f$ as
       @f[
          K_{iB} =  K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})
       @f]
     * the @f$ \hat{\mathbf{X}}_{i} @f$ denotes the description of a certain atom i
     * within it's cutoff radius by the use of the descriptors stored in the descriptorCollector
     * object.
     *******************************************************************************************/
    void compute_kernelMatrix();
    /*******************************************************************************************
     * allocating kernel matrix as [types][ atoms per type x number reference configurations per
     *type ]
     *
     * also sets matrix elements to zero by using std::fill after doing the resize
     *******************************************************************************************/
    void allocate_kernelMatrix(const Vec1Int& nAtomPerType);
    /*******************************************************************************************
     * Resize Kernel matrix for GPU usage. Copy check is first done to avoid unneccessary host
     * device copies.
     *******************************************************************************************/
    void resizeKernelMatrixGPU(const Vec1Int& nAtomPerType);
    /*******************************************************************************************
     * Compute the needed dimensions for the kernel matrix.
     *******************************************************************************************/
    void compute_kernelMatrixDim(const Vec1Int& nAtomPerType);
    /*******************************************************************************************
     * @brief Populates descriptor configurations in shared memory based on input parameters.
     *
     * This function iterates over a predefined list of descriptor keys and updates the
     * corresponding shared memory arrays with values extracted from the input parameters.
     * Only descriptors with non-zero weights are processed.
     *
     * @param inputParameters A reference to the input record containing parameter values.
     *******************************************************************************************/
    void fillDescriptor(const Record& inputParameters);

  protected:
    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /*******************************************************************************************
     * @brief Shared pointer to an instance of the MlMPI class.
     *
     * This member holds a shared ownership of an MlMPI object.
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mpi;
    /*******************************************************************************************
     * storing weights for 2-body (DescriptorSHS2),3-body(DescriptorSHS3)...
     * descriptors ( are not regression coefficients )
     *******************************************************************************************/
    std::map<String, Real> weights;
    /// number of local reference configurations per force field type
    ShVec1Int numberLocalRefConfs;
    /// storing the number of basis functions for 2-body, 3-body... terms
    std::map<String, ShVec1Int> numberDescriptors;
    /*******************************************************************************************
     * feature space size for 2-body, 3-body...
     *
     * featureSpace size is the number of local reference configurations times
     * number of basis functions for a certain descriptor. featureSpaceSize can be retrieved
     * per descriptor. The type of descriptor is defined by the string of the map
     *******************************************************************************************/
    std::map<String, ShVec1Int> featureSpaceSize;
    /*******************************************************************************************
     * local reference configurations described in descriptor space of 2-body, 3body...
     * descriptors map can chose from 2-body and 3-body... descriptors by keys
     * defined in constants.hpp
     *******************************************************************************************/
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> descriptorsRefConfs;
    /*******************************************************************************************
     * collector class to store all possible dscriptors.
     *
     * @note for the here used kernel the desciptors in the collector class
     * have to be type ordered
     *******************************************************************************************/
    std::shared_ptr<DescriptorCollector> descriptorCollection;
    /*******************************************************************************************
     * Stores kernel matrix.
     *
     * Kernel matrix is defined by a kernel function @f$ K(:,\hat{\mathbf{X}}_{B}) @f$ as
     *
     * @f[
          K_{iB} =  K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})
     * @f]
     * the @f$ \hat{\mathbf{X}}_{i} @f$ denotes the description of a certain atom i
     * within it's cutoff radius by the use of the descriptors stored in the descriptorCollector
     * object.
     *
     * The kernelMatrix has structure atom types as first index and the
     * second index is a combination of atoms per type in the structure
     * times the local reference configurations
     *******************************************************************************************/
    Vec2Real& kernelMatrix;
    /*******************************************************************************************
     * Array to store the kernel matrix sizes for the different atom types.
     *******************************************************************************************/
    Vec1Size kernelMatrixDim;
    /*******************************************************************************************
     * Object to controll if the kernel matrix is need of a resize.
     *******************************************************************************************/
    ArrayResizing2D kernelMatrixSize;
    /*******************************************************************************************
     * typeMap which assigns the force field atom types to the types in the structure and vice versa
     *******************************************************************************************/
    std::shared_ptr<const TypeMap> typeMap;
    /*******************************************************************************************
     * Local copy of neighbor list's vector containing number of atoms per type in structure.
     *******************************************************************************************/
    Vec1Int nAtomsPerType;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapKernel();

} //namespace vaspml

#endif
