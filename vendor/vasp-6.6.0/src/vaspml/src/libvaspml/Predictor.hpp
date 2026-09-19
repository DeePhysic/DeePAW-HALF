#ifndef PREDICTOR_HPP
#define PREDICTOR_HPP

#include "ArrayResizing.hpp"
#include "constants.hpp"
#include "types.hpp"

#include <tuple>

namespace vaspml
{

class DescriptorCollector;
class Kernel;
class MlMPI;
class TypeMap;
struct Record;
template <class T> class ShmemArray2DVariableLen;

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> make_descriptorsRefConfsWeighted(
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::map<String, Real>&      weights,
    const std::shared_ptr<MlMPI>&      mpiIn);

std::map<String, ShVec2Real> makeMaps2D();

std::map<String, ShVec1Real>                   makeMaps1D();
std::map<String, Vec1Size>                     makeMaps1DSize();
std::map<String, ArrayResizing1D>              makeResizeMap1D();
std::map<String, ArrayResizing2D>              makeResizeMap2D();

enum class TotalEnergyType
{
    IsolatedAtom,
    AverageTrainEnergy,
};

/*******************************************************************************************
 * functor for computing the total atomic force acting on an ion
 *
 * Create by SumFunctorForce functor( atomicForces, pairForces[ centralAtom ] ); and use
 * with a std algorithm as
 * std::for_each( neighborIndex[ centralAtom ].cbegin(), neighborIndex[ centralAtom ].cend(),
 *functor ); This procedure has to be called from an loop over central atoms. One central atom per
 *time can be processed with this routines.
 *******************************************************************************************/
class SumFunctorForce
{
  public:
    SumFunctorForce(Vec1Real& atomicForce, const Vec1Real& pairForce) :
        atomicForce(atomicForce),
        pairForce(pairForce),
        counter(0)
    {}
    void operator()(const Int& neighborIndex)
    {
        atomicForce[3 * neighborIndex] += pairForce[3 * counter];
        atomicForce[3 * neighborIndex + 1] += pairForce[3 * counter + 1];
        atomicForce[3 * neighborIndex + 2] += pairForce[3 * counter + 2];
        counter++;
    }

  private:
    /*******************************************************************************************
     * total atomic force.
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... atomicForce array will always point to the atomicForce
     * memory of the components of Predictor::atomicForce
     *******************************************************************************************/
    Vec1Real& atomicForce;
    /*******************************************************************************************
     * pair atomic force.
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of pairForces
     * memory of the components of Predictor::pairForces
     *******************************************************************************************/
    const Vec1Real& pairForce;
    /*******************************************************************************************
     * counting the number of neighbor atoms which have been processed already
     *******************************************************************************************/
    Int counter = 0;
};
/*******************************************************************************************
 * functor computing the stress tensor contribution due to a single atom
 *
 * Create by
 * SumFunctorStress functor( stressTensor, pairForces[ centralAtom ],connectionVector[ centralAtom
 *], inverseVolume ); Then use by a std algorithm as std::for_each( neighborIndex[ centralAtom
 *].cbegin(),neighborIndex[ centralAtom ].cend(), functor ); This will compute the contribution of a
 *single central atom to the stress tensor. To get total stress repeat procedure for every central
 *atom.
 *******************************************************************************************/
class SumFunctorStress
{
  public:
    SumFunctorStress(Vec1Real&       stressTensor,
                     const Vec1Real& pairForces,
                     const Vec1Real& connectionVector,
                     const Real&     inverseVolume) :
        stressTensor(stressTensor),
        pairForces(pairForces),
        connectionVector(connectionVector),
        inverseVolume(inverseVolume),
        counter(0)
    {}

    void operator()(const Size)
    {
        // compute xx term
        stressTensor[0] -= inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter];
        // compute xy term
        stressTensor[1] -=
            inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter + 1];
        // compute xz term
        stressTensor[2] -=
            inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter + 2];
        // compute yy term
        stressTensor[4] -=
            inverseVolume * pairForces[3 * counter + 1] * connectionVector[3 * counter + 1];
        // compute yz term
        stressTensor[5] -=
            inverseVolume * pairForces[3 * counter + 1] * connectionVector[3 * counter + 2];
        // compute zz term
        stressTensor[8] -=
            inverseVolume * pairForces[3 * counter + 2] * connectionVector[3 * counter + 2];
        counter++;
    }

  private:
    Vec1Real& stressTensor;
    /*******************************************************************************************
     * Pair atomic force. Stores force between pair of atoms
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of pairForces
     * memory of the components of Predictor::pairForces
     *******************************************************************************************/
    const Vec1Real& pairForces;
    /*******************************************************************************************
     * connection vector array between atom pairs
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of
     * connectionVector memory of the components of NearestNeighborNSquare::connectionVector
     *******************************************************************************************/
    const Vec1Real& connectionVector;
    /*******************************************************************************************
     * inverse volume of the current box used for the simulation
     *******************************************************************************************/
    const Real& inverseVolume;
    /*******************************************************************************************
     * counting the number of neighbor atoms which have been processed already
     *******************************************************************************************/
    Int counter = 0;
};

struct StressTensorLookUpTable
{
    /*******************************************************************************************
     * Main index to loop over for stress computation.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * Stroes central atom index for all combinations of descriptor index and central atom index
     *******************************************************************************************/
    Vec1Size centralAtomIndex;
    /*******************************************************************************************
     * Stores the descriptor index for all combinations of descriptor index and central atom
     *******************************************************************************************/
    Vec1String descriptorIndex;
    /*******************************************************************************************
     * Number of atom in the previous call. Needed to determine if size of index arrays
     * has to change.
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Update all index arrays.
     *
     * @param nAtoms number of atoms in current structure.
     *******************************************************************************************/
    void update(const Size& nAtoms);
    /*******************************************************************************************
     * Refilling the index arrays if number of central atoms changed.
     *
     * @param nAtoms number of central atoms.
     *******************************************************************************************/
    void refill(const Size& nAtoms);
};

struct AtomicForcesLookUpTable
{
    Vec1Size mainIndex;
    Vec1Size centralAtomIndex;
    Vec1Int  globalIndex;
    Size     nAtomsOld = 0;
    void     update(const Int& nAtoms);
    void     refill(const Int& nAtoms);
};

/*******************************************************************************************
 * @class Predictor
 * @brief Compute energy, forces and stress tensor from some supplied Kernel
 *
 * Predictor is a class which stores the fitting weights (regression coeffcients) obtained
 * during Kernel training. When a Kernel is supplied to the Predictor class the programmer
 * will obtain energy, forces and stresses for the structure of which the supplied
 * Kernel was computed. The Predictor class is agnostic to the type of the supplied Kernel.
 * The class is also agnostic to the types of the Descriptors from which the Kernel
 * classes are made. To keep the agnosticism any Descriptors have to inherit from Descriptor
 * and any used Kernel has to inherit from the virtual Kernel base class. When done
 * in this way the Predictor class does not be touched to implement new Kernel or
 * Descriptor functionality.
 *******************************************************************************************/
class Predictor
{
  public:
    /*******************************************************************************************
     * default construction of Predictor class.
     *
     * @note class is not usable in this state, but can be transformed into a usable state by
     * calling Predictor::init().
     *******************************************************************************************/
    Predictor(ShRec record = nullptr);
    /*******************************************************************************************
     * @brief Constructor for the Predictor class.
     *
     * Initializes the Predictor object with input parameters, MPI communication,
     * execution policy, and shared record data. Sets up various internal data
     * structures and configurations required for machine learning-based predictions.
     *
     * @param inputParameters Input parameters containing configuration and data for the predictor.
     * @param mpiIn Shared pointer to the MPI communication object for parallel execution.
     * @param record Shared record object for storing and accessing data.
     *******************************************************************************************/
    Predictor(const Record&                 inputParameters,
              const std::shared_ptr<MlMPI>& mpiIn = nullptr,
              ShRec                         record = nullptr);
    /*******************************************************************************************
     * @brief Initializes the Predictor object with input parameters, MPI instance, and execution
     *policy.
     *
     * This function sets up various internal data structures and parameters required for the
     *Predictor to perform its tasks. It processes input parameters, initializes maps, vectors, and
     *other resources, and configures the energy scaling type based on the provided input.
     *
     * @param inputParameters A Record object containing the input parameters for initialization.
     * @param mpiIn A shared pointer to an MlMPI instance for parallel processing support.
     *
     * @details
     * - Reads and processes input parameters such as maximum types, local reference configurations,
     *   descriptors, feature space size, weights, regression coefficients, and reference energies.
     * - Initializes various internal data structures, including derivative matrices, force maps,
     *   and stress tensors.
     * - Configures the energy scaling type based on the "scale-type-total-energy" parameter.
     * - Precomputes weighted descriptors for reference configurations.
     *
     * @note This function must be called before using the Predictor object for any computations.
     *******************************************************************************************/
    void init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn = nullptr);
    /*******************************************************************************************
     * Compute total energy and forces for current structure supplied by kernel.
     *
     * For the forces only the arrays centralForces and pairForces will be computed.
     * From those the total ionic force can be computed with member function
     * compute_atomicForces( const DescriptorCollector& descriptorCollection );
     *
     * @param kernel is storing the similarity measure of the current structure
     * to the reference configurations
     * kernel also has to supply the derivatives ( derivativeMatrix ) of the kernel matrix
     *******************************************************************************************/
    void update(const Kernel& kernel);
    /*******************************************************************************************
     *compute total energy of system and fill energyArray
     *
     @param kernel supplies the kernel matrix form which atomic energy is computed
    *******************************************************************************************/
    void updateEnergy(const Kernel& kernel);
    /*******************************************************************************************
     * compute force arrays centralForces and pairForces
     *
     * @param kernel supplies the kernel matrix and kernel derivatives needed
     * for computation.
     * @note the kernel has to be updated outside of the current class
     *
     * @f[ F^{i}_{\alpha}=-\sum_{B}w_{B}\frac{\partial
     *K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}{\partial r^{i}_{\alpha}} @f] centralForces
     * @f[ F^{j}_{\alpha}=\sum_{B}\sum_{i=1}^{N_{atoms}}w_{B}\frac{\partial
     *K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}{\partial r^{i}_{\alpha}} @f] pairForces
     *******************************************************************************************/
    void updateForces(const Kernel& kernel);
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Depending on the supplied execution
     * policy the function will choose between CPU and GPU execution
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForces(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * Compute total forces acting on ions and store to a global array which contains
     * all atoms not only those of current MPI island.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Depending on the supplied execution
     * policy the function will choose between CPU and GPU execution
     *
     * @param descriptorCollection contains descriptor used for kernel computation. Used to retrieve
     * neighbor lists
     * @param globalIonIndex mapping between the global ion index and the ion index on the local
     *rank
     * @param nAtomsGlob total number of atoms
     *******************************************************************************************/
    void compute_atomicForces(const DescriptorCollector& descriptorCollection,
                              const Vec1Int&             globalIonIndex,
                              const Int&                 nAtomsGlob);
    /*******************************************************************************************
     * return x component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesX(const String& key,
                                const Size    centralAtom,
                                const Size    neighborAtom) const;
    /*******************************************************************************************
     * return y component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesY(const String& key,
                                const Size    centralAtom,
                                const Size    neighborAtom) const;
    /*******************************************************************************************
     * return z component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesZ(const String& key,
                                const Size    centralAtom,
                                const Size    neighborAtom) const;
    /*******************************************************************************************
     * Return const reference to pairForces.
     *******************************************************************************************/
    const std::map<String, ShVec2Real>& get_pairForces() const;
    /*******************************************************************************************
     * return x,y,z component of pairForces for certain ion-ion pair
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     * @note for performance call with auto& [ x, y, z ] = get_pairForce( centralAtom, neighborAtom
     *)
     *******************************************************************************************/
    const std::tuple<const Real&, const Real&, const Real&>
    get_pairForces(const String& key, const Size centralAtom, const Size neighborAtom) const;
    /*******************************************************************************************
     * return x component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesX(const String& key, const Size centralAtom) const;
    /*******************************************************************************************
     * return y component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesY(const String& key, const Size centralAtom) const;
    /*******************************************************************************************
     * return z component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesZ(const String& key, const Size centralAtom) const;
    /*******************************************************************************************
     * return x,y,z component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     * @note for performance call with
     * auto& [ x, y, z ] = get_centralForces( const Size centralAtom, const Size
     *neighborAtom )
     *******************************************************************************************/
    const std::tuple<const Real&, const Real&, const Real&> get_centralForces(
        const String& key,
        const Size    centralAtom) const;
    /*******************************************************************************************
     * Get constant reference to centralForces array.
     *******************************************************************************************/
    const std::map<String, ShVec1Real>& get_centralForces() const;
    /*******************************************************************************************
     * Return x component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesX(const Size centralAtom) const;
    /*******************************************************************************************
     * Return y component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesY(const Size centralAtom) const;
    /*******************************************************************************************
     * Return z component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesZ(const Size centralAtom) const;
    /*******************************************************************************************
     * Get atomic force array in flattened 3*n form. Order is x1,y1,z1,x2,y2,z2...xn,yn,zn
     *******************************************************************************************/
    const Vec1Real& get_atomicForces() const;
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * computations are done either in
     * compute_stressTensorCPU( const DescriptorCollector& descriptorCollection, const Real& volume
     *); or compute_stressTensorGPU( const DescriptorCollector& descriptorCollection, const Real&
     *volume ); depending on the chosen execution policy
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensor(const DescriptorCollector& descriptorCollection, const Real& volume);
    /*******************************************************************************************
     * Getter for stressTensor component per descriptor.
     *
     * @param key keyword defining for which descriptor the partial stress tensor should be obtained
     * @param indx0 first component of stress tensor for which value is retrieved
     * @param indx1 first component of stress tensor for which value is retrieved
     *
     * @note the storage order is xx, xy, xz, yx, yy, yz, zx, zy, zz
     *                            00  01  02  10  11  12  20  21  22
     *******************************************************************************************/
    const Real& get_stressTensor(const String& key, const Size indx0, const Size indx1) const;
    /*******************************************************************************************
     * Getter for total stress tensor component (totalStressTensor).
     *
     * @param indx0 first component of stress tensor for which value is retrieved
     * @param indx1 first component of stress tensor for which value is retrieved
     *
     * @note the storage order is xx, xy, xz, yx, yy, yz, zx, zy, zz
     *                            00  01  02  10  11  12  20  21  22
     * the total stress tensor was summed over the descriptors
     *******************************************************************************************/
    const Real& get_totalStressTensor(const Size indx0, const Size indx1) const;
    /*******************************************************************************************
     * @brief Retrieves the total stress tensor.
     *
     * This function provides access to the total stress tensor, which is
     * represented as a `Vec1Real` object. The tensor encapsulates the
     * stress state of the system.
     *
     * @return A constant reference to the total stress tensor.
     *******************************************************************************************/
    const Vec1Real& get_totalStressTensor() const;
    /*******************************************************************************************
     * get total energy of the of considered system
     *******************************************************************************************/
    const Real& get_totalEnergy() const;
    /*******************************************************************************************
     * writing atomic forces to screen
     *******************************************************************************************/
    void write_atomicForceToScreen() const;
    /*******************************************************************************************
     * Write tempForce vector in flattened form to file.
     *
     * @param fname name of file to write to.
     *******************************************************************************************/
    void write_tempForceVector(const String& fname) const;
    /*******************************************************************************************
     * Write derivative matrix in flattened 1 column format to file.
     *
     * @param baseName base name which will be used and extended by SHS2 or SHS3 for body terms
     *******************************************************************************************/
    void write_derivativeMatrix(const String& baseName) const;
    /*******************************************************************************************
     * Write forcePreContract to file in flattened 1D column format.
     *
     * @param baseName base name which will be used and extended by SHS2 or SHS3 for body terms
     *******************************************************************************************/
    void write_forcePreContract(const String& baseName) const;
    /*******************************************************************************************
     * Write pairForces to file in flattened 1D column format.
     *
     * @param baseName base name which will be used and extended by SHS2 or SHS3 for body terms
     *******************************************************************************************/
    void write_pairForces(const String& baseName) const;
    /*******************************************************************************************
     * Write centralforces array to file in flattened 1D column format.
     *
     * @param baseName base name which will be used and extended by SHS2 or SHS3 for body terms
     *******************************************************************************************/
    void write_centralForces(const String& baseName) const;

  private:
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on CPU
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * Compute total forces acting on ions and write terms to a globally allocated array.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on CPU
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     * @param globalIonIndex mapping between mpi distributed ion indexes to index in global array
     *******************************************************************************************/
    void compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection,
                                 const Vec1Int&             globalIonIndex);
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on GPU
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * Compute total forces acting on ions and write terms to a globally allocated array.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on CPU
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     * @param globalIonIndex mapping between mpi distributed ion indexes to index in global array
     * @param nAtomsGlobal total number of atoms. Not mpi distributed number of atoms.
     *******************************************************************************************/
    void compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection,
                                 const Vec1Int&             globalIonIndex,
                                 const Size&                nAtomsGlobal);
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * This algorithm will be executed on a CPU
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensorCPU(const DescriptorCollector& descriptorCollection,
                                 const Real&                volume);
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * This algorithm will be executed on a GPU
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensorGPU(const DescriptorCollector& descriptorCollection,
                                 const Real&                volume);
    /*******************************************************************************************
     * Sum over the atom in the structure per type and normalize by 1/number ions of the kernel
     matrix.
     *
     *@param kernel used to retrieve kernel matrix @f$K_{iB}@f$
     *@param typeStruc analysed atom type in structure. Atom type of atom number i
     *@param typeForceField analysed atom type in force field. Atom type of atom number i in FF
     *@param atomsPerType number of atom per type.
     *
     *
     *@f[ K( \hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B} )=K_{iB}\in
          \mathbb{R}^{N_{atoms\ per\ type}\times N_{reference\ config}}@f]
     *@f[\mathbf{E}_{type}=\sum_{i=1}^{N_{atoms\ per\ type}}\frac{K_{iB}}{N_{atoms}}@f]
     *energyArray has storage form number of types x number reference configurations
     *@f[\mathbf{E} \in \mathbb{R}^{N_{types}\times N_{reference\ config}}@f]
     *******************************************************************************************/
    void compute_energyArray(const Kernel& kernel,
                             const Size    typeStruc,
                             const Size    typeForceField,
                             const Size    atomsPerType,
                             const Real    scaleFactor);
    /*******************************************************************************************
     * @brief Computes weighted descriptors for reference configurations.
     *
     * This function calculates the weighted descriptors for reference configurations
     * based on the provided input parameters, regression coefficients, and descriptor
     * values. The computation is performed for each descriptor key in the predefined
     * descriptor key list. The results are stored in the `descriptorsRefConfsWeighted`
     * data structure.
     *
     * @param inputParameters A `Record` object containing the input parameters, including
     *                        regression coefficients and descriptor values.
     *
     * The function performs the following steps:
     * - Iterates over all descriptor keys in `constants::descriptorKeyList`.
     * - Checks if the weight for the current descriptor key is greater than zero.
     * - Retrieves the regression coefficients and descriptor values for the current key.
     * - Computes the weighted descriptors for each type, local reference configuration,
     *   and descriptor count.
     * - Updates the `descriptorsRefConfsWeighted` data structure with the computed values.
     * - Synchronizes data across shared memory and inter-node communication if `use_shmem`
     *   is enabled.
     *
     * @note The function assumes the presence of shared memory (`use_shmem`) and MPI
     *       communication for distributed environments. Synchronization and data
     *       distribution are handled accordingly.
     *******************************************************************************************/
    void compute_descriptorsRefConfsWeighted(const Record& inputParameters);
    /*******************************************************************************************
     * Compute the total reference energy of considered system.
     *
     * Result is stored to referenceEnergyTotal
     *
     * @f[E_{ref}=\sum_{n=0}^{number\ types}N_{n}*E_{n}^{ref}@f]
     *******************************************************************************************/
    void compute_referenceEnergyTotal(const Vec1Int& atomsPerType);
    /*******************************************************************************************
     * allocating the energy array
     *
     * first dimension of energy array will be the number of types
     * the second dimension is the number of local reference configurations
     *
     * @param numberTypesStruc number of types in structure
     *******************************************************************************************/
    void allocate_energyArray(const Size numberTypesStruc);
    /*******************************************************************************************
     * allocating the energy array for CPU using the resize functionality of the std::vector
     *
     * @param numberTypesStruc number of types in structure
     *******************************************************************************************/
    void allocate_energyArrayCPU(const Size numberTypesStruc);
    /*******************************************************************************************
     * allocating the energy array for GPU using hand crafted resize objects to avoid host
     * device data copies if not needed.
     *
     * @param numberTypesStruc number of types in structure
     *******************************************************************************************/
    void allocate_energyArrayGPU(const Size numberTypesStruc);
    /*******************************************************************************************
     * Compute the dimensions of the energy arrays. Energy arrays store temporary results for
     * total energy computation.
     *******************************************************************************************/
    void compute_energyArrayDim(const Size numberTypesStruc);
    /*******************************************************************************************
     * allocating the derivativeMatrix which stores the derivative of the specific kernel
     *
     * the dimension of the derivativeMatrix is the number of types in the first dimension
     * the second dimension is the number of atoms per type times the number of descriptors
     * the number of descriptors changes faster in the second dimension
     *
     * @param numberAtomType number of atoms per type
     *******************************************************************************************/
    void allocate_derivativeMatrix(const Vec1Int& numberAtomType);
    /*******************************************************************************************
     * Allocate derivative matrix on CPU. Use resize function of std::vector if number of atoms
     * per type changes.
     *
     * @param numberAtomType number of atoms per type.
     *******************************************************************************************/
    void allocate_derivativeMatrixCPU(const Vec1Int& numberAtomType);
    /*******************************************************************************************
     * Allocate derivative matrix on CPU. Use resize function of std::vector if number of atoms
     * per type changes.
     *
     * Hand craftet resize function is used to avoid unneeded host device copies
     *
     * @param numberAtomType number of atoms per type.
     *******************************************************************************************/
    void allocate_derivativeMatrixGPU(const Vec1Int& numberAtomType);
    /*******************************************************************************************
     * Determine dimensions of derivativeMatrix
     *******************************************************************************************/
    void compute_derivativeMatrixDim(const Vec1Int& numberAtomType);
    /*******************************************************************************************
     * allocating the derivativeMatrixDescriptorDerivativeProduct
     *
     * the first dimension will be the number of atom in the supplied structure
     * the second dimension is the number of descriptors for the specific atom
     *
     * @param descriptorCollection collection of descriptors used for kernel computation
     * @param numberAtoms number of central atoms.
     *******************************************************************************************/
    void allocate_forcePreContract(const DescriptorCollector& descriptorCollection,
                                   const Size                 numberAtoms);
    /*******************************************************************************************
     * allocating the derivativeMatrixDescriptorDerivativeProduct
     * the first dimension will be the number of atom in the supplied structure
     * the second dimension is the number of descriptors for the specific atom
     *
     * @param numberAtoms total number of central atoms.
     *******************************************************************************************/
    void allocate_forcePreContractCPU(const Size numberAtoms);
    /*******************************************************************************************
     * allocating the derivativeMatrixDescriptorDerivativeProduct
     * the first dimension will be the number of atom in the supplied structure
     * the second dimension is the number of descriptors for the specific atom
     *
     * Resize check will be done to avoid unneccessary host device copies.
     *
     * @param numberAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_forcePreContractGPU(const Size numberAtoms);
    /*******************************************************************************************
     * Compute the size of the dimensions for the array forcePreContract
     *
     * Size is equal to the DescriptorSHS2 size per atom
     *
     * @param descriptorCollection descriptor collection to get descriptor classes
     * @param numberAtoms total number of central atoms
     *******************************************************************************************/
    void compute_forcePreContractDim(const DescriptorCollector& descriptorCollection,
                                     const Size                 numberAtoms);
    /*******************************************************************************************
     * allocating the arrays for the pair forces
     *
     * @param descriptorCollection retrieve neighbor list to get number of neighbors
     * per descriptor and central atom
     * size is chosen in the same way as the connectionVector in neighborList
     * first index is neighbor second index is combination of neighbor xyz
     * neighbor 1 x
     * neighbor 1 y
     * neighbor 1 z
     * neighbor 2 x
     * neighbor 2 y
     * neighbor 2 z
     * .
     * .
     * .
     *******************************************************************************************/
    void allocate_pairForces(const DescriptorCollector& descriptorCollection,
                             const Size&                numberAtoms);
    /*******************************************************************************************
     * Allocating the arrays for the pair forces for CPU variant.
     *
     * @param numberAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_pairForcesCPU(const Size& numberAtoms);
    /*******************************************************************************************
     * Allocating the arrays for the pair forces for CPU variant.
     *
     * Avoiding GPU CPU copies of data array
     *
     * @param numberAtoms total number of central atoms.
     *******************************************************************************************/
    void allocate_pairForcesGPU(const Size& numberAtoms);
    /*******************************************************************************************
     * Computing size for pair forces.
     *
     * Size of pair forces for given central atom is 3xnumber neighbor atoms
     *
     * @param descriptorCollection collection of descriptors to retrieve neigbhor lists.
     *******************************************************************************************/
    void compute_pairForcesDim(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * allocate central force arrays
     *
     * allocating in size number central atoms times 3 (for xyz)
     *******************************************************************************************/
    void allocate_centralForces(const Int numberAtoms);
    /*******************************************************************************************
     * Allocate central force arrays for CPU variant.
     *
     * Just uses std::vector resize functionality
     *
     * @param numberAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_centralForcesCPU(const Int numberAtoms);
    /*******************************************************************************************
     * Allocate central force arrays for CPU variant.
     *
     * Just uses std::vector resize functionality
     *
     * @param numberAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_centralForcesGPU(const Int numberAtoms);
    /*******************************************************************************************
     * Allocate centralAtomIndex containing a iota from 0 to number of atoms.
     *
     * @param numberAtoms total number of atoms.
     *******************************************************************************************/
    void allocate_centralAtomIndex(const Size numberAtoms);
    /*******************************************************************************************
     *allocate stressTensor and totalStressTensor
     *
     * allocation will be done as flattened 3x3 array of size 9
     *******************************************************************************************/
    void allocate_stressTensor();
    /*******************************************************************************************
     *allocate stressTensor and totalStressTensor for CPU usage
     *******************************************************************************************/
    void allocate_stressTensorCPU();
    /*******************************************************************************************
     *allocate stressTensor and totalStressTensor for GPU usage
     *******************************************************************************************/
    void allocate_stressTensorGPU();
    /*******************************************************************************************
     * allocate atomicForces
     *
     * allocation will be done as flattened Natoms x 3 array
     *******************************************************************************************/
    void allocate_atomicForces();
    /*******************************************************************************************
     * allocate atomicForces for CPU version
     *
     * allocation will be done as flattened Natoms x 3 array
     *******************************************************************************************/
    void allocate_atomicForcesCPU();
    /*******************************************************************************************
     * allocate atomicForces for GPU version
     *
     * allocation will be done as flattened Natoms x 3 array
     *******************************************************************************************/
    void allocate_atomicForcesGPU();
    /*******************************************************************************************
     * Allocate atomicForces in global form.
     *
     * Allocation will be done for all ions in system not only those present on current mpi rank.
     * Allocation will be done as flattened Natoms x 3 array.
     *
     *@param nAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_atomicForces(const Size& nAtoms);
    /*******************************************************************************************
     * Allocate atomicForces in global form for CPU version.
     *
     * Allocation will be done for all ions in system not only those present on current mpi rank.
     * Allocation will be done as flattened Natoms x 3 array.
     *
     *@param nAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_atomicForcesCPU(const Size& nAtoms);
    /*******************************************************************************************
     * Allocate atomicForces in global form for GPU version.
     *
     * Allocation will be done for all ions in system not only those present on current mpi rank.
     * Allocation will be done as flattened Natoms x 3 array.
     *
     *@param nAtoms total number of central atoms
     *******************************************************************************************/
    void allocate_atomicForcesGPU(const Size& nAtoms);
    /*******************************************************************************************
     * Allocate tempForceVector needed for computation of derivativeMatrix.
     *
     * @param nAtomsType number of atoms per type
     *******************************************************************************************/
    void allocate_tempForceVector(const Vec1Int& nAtomsType);
    /*******************************************************************************************
     * Allocate tempForceVector needed for computation of derivativeMatrix in CPU version.
     *
     * @param nAtomsType number of atoms per type
     *******************************************************************************************/
    void allocate_tempForceVectorCPU();
    /*******************************************************************************************
     * Allocate tempForceVector needed for computation of derivativeMatrix in GPU version.
     *
     * @param nAtomsType number of atoms per type
     *******************************************************************************************/
    void allocate_tempForceVectorGPU();
    /*******************************************************************************************
     * Determine size of tempForceVector. Is number of atoms per type.
     *
     * @param nAtomsType number of atoms per type
     *******************************************************************************************/
    void compute_tempForceVectorDim(const Vec1Int& nAtomsType);
    /*******************************************************************************************
     * Computing the energy contribution to the total energy of a specific atom type.
     *
     * computes the dot product between the regression coefficients and energyArray
     @f$\mathbf{E}_{type}@f$
     * @f[\sum\limits_{B}w_{B}\left [ \sum\limits_{i} \frac{1}{N_{\mathrm{atom}}} K(\mathbf{X}_{i},
         \mathbf{X}_{B}) \right ]=\sum\limits_{B}w_{B}E_{type}^{B}@f]
     *******************************************************************************************/
    Real computeEnergyPerType(const Size typeStruc, const Size typeForceField);
    /*******************************************************************************************
     * Adding reference energies to the total energies.
     *
     * depending on the tag TotalEnergyType the energy will be scaled to the isolated atoms
     * or the average energy of the system
     *******************************************************************************************/
    void finalizeEnergyComputation(const Vec1Int& atomsPerType);
    /*******************************************************************************************
     * Computing the derivative of the kernel matrix
     *
     * derivative matrix @f$\mathbf{L}^{i}=\{L_{ln_{1}n_{2}}^{iJ_{1}J_{2}}\}@f$
     * @f[
        \mathbf{L}^{i} = -\sum_{B} w_{B}
               \frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{i})}
                    {\partial \hat{\mathbf{X}}_{i}}
       @f]

     * @note for details of the computation of the derivative matrix please check
     * the documentation of the used Kernel class, as for example KernelPolynomial.
     * Function which implements functionality will be called computeDerivativeKernel
     *******************************************************************************************/
    void compute_derivativeMatrix(const Kernel& kernel,
                                  Vec1Real&     tempForceVector,
                                  const Size    typeStruc,
                                  const Size    typeForceField);
    /*******************************************************************************************
     * Computes the derivative matrix times the derivatives of the descriptors w.r.t SHS2 and SHS3.
     *
     * @f[\mathbf{\tilde{L}}=\sum\limits_{i}\mathbf{L}_{i} \frac{d\mathbf{X}_{i}}{dp_{nn'l}^{iJJ'}}
     * \frac{dp_{nn'l}^{iJJ'}}{dc_{n''lm}^{iJ''}}@f]
     *
     * dimensions of derivativeMatrix @f$\mathbf{\tilde{L}}\in \mathbb{R}^{N_{atoms}\times lnm }@f$
     *
     * @note implementation is descriptor specifix and the implementing code can be found
     * in the descriptor classes DescriptorSHS2 and DescriptorSHS2
     *******************************************************************************************/
    void compute_forcePreContract(const DescriptorCollector& descriptorCollection,
                                  const TypeMap&             typeMap);
    /*******************************************************************************************
     * Computing force arrays pairForces and centralForces
     *
     * centralForces
     * @f[-\frac{\partial K(\mathbf{X}_{i},\mathbf{X}_{B})}{\partial \mathbf{X}_{i}}*
     *    \frac{d\ \mathbf{X}_{i}}{d\ \mathbf{r}_{i}}@f]
     *
     * pairForces
     * @f[-\sum_{k\in \mathcal{N}_{i}}\frac{\partial K(\mathbf{X}_{i},\mathbf{X}_{B})}{\partial
     *\mathbf{X}_{i}}
     * \frac{d\ \mathbf{X}_{i}}{d\ \mathbf{r}_{k}}@f]
     *
     * @note implementation details can be found in DescriptorSHS2. There the final step is done
     * which is the contraction of forcePreContract and the pair descriptors.
     *******************************************************************************************/
    void computeForceArrays(const DescriptorCollector& descriptorCollection);

    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /// number of different atom types stored in force field
    Int numberTypesForceField;
    /*****************make2DRandReal**************************************************************************
     * number of local reference configurations per force field type
     *******************************************************************************************/
    ShVec1Int numberLocalRefConfs;
    /*******************************************************************************************
     * storing the number of descriptors for 2-body, 3-body... terms
     *******************************************************************************************/
    std::map<String, ShVec1Int> numberDescriptors;
    /*******************************************************************************************
     *  feature space size for 2-body, 3-body... descriptors;
     *  feature space is number of local configurations times number of basis functions
     *******************************************************************************************/
    std::map<String, ShVec1Int> featureSpaceSize;
    /*******************************************************************************************
     * storing weights ( not regression coefficients ) for 2-body,3-body... descriptors
     *******************************************************************************************/
    std::map<String, Real> weights;
    /*******************************************************************************************
     * regression coefficients stored in format (types,local-reference-configurations)
     *******************************************************************************************/
    std::shared_ptr<ShmemArray2DVariableLen<Real>> regressionCoefficients;
    /*******************************************************************************************
     * reference energies of different atom types.
     *
     * @note this can be considered as the energies of a certain atom type in vaccuum
     *******************************************************************************************/
    ShCVec1Real referenceEnergyPerType;
    /*******************************************************************************************
     * average total energy which was ecountered during training of the force field
     *******************************************************************************************/
    Real averageTrainEnergy;
    /*******************************************************************************************
     * SmartEnum to determine if total energy is expressed as isolated atom energy or
     * with respect to average train energy
     *******************************************************************************************/
    TotalEnergyType energyType;
    /*******************************************************************************************
     * total reference energy for current system defined by the supplied Kernel class
     *******************************************************************************************/
    Real referenceEnergyTotal;
    /*******************************************************************************************
     * stores total energy of system after the update routine finished successfully
     *******************************************************************************************/
    Real totalEnergy;
    /*******************************************************************************************
     * typeMap which assigns the force field atom types to the types in the structure and vice versa
     *******************************************************************************************/
    std::shared_ptr<const TypeMap> typeMap;
    /*******************************************************************************************
     * Storing the local reference configurations expressed as descriptors times the regression
     *coefficients.
     *
     * @f[\mathbf{X}_{i_{B}} \in \mathbb{R}^{number-types\times feature-space-size}@f]
     * @f[\mathbf{w}_{i_{B}} \in \mathbb{R}^{number-types\times local-reference-configurations}@f]
     * @f[\mathbf{X}^{'}_{B}=\mathbf{X}_{B} \circ \mathbf{w}_{B}@f]
     *
     * the storage order is types, local reference configurations x numberDescriptors.
     * descriptorsRefConfsWeighted is @f$\mathbf{X}^{'}_{B}\in\mathbb{R}^{N_{types}\times
     * N_{local\ ref\ conf\times N_{descriptors}}}=\mathbb{R}^{N_{types}\times
     *N_{feature\ space}}@f$ numberDescriptors changes fastest in the squashed index
     *******************************************************************************************/
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> descriptorsRefConfsWeighted;
    /*******************************************************************************************
     * Array to store the contraction of the kernel matrix over the atoms times 1/number-ions.
     *
     * @f[
        \sum\limits_{B}w_{B}\sum\limits_{i} \frac{1}{N_{\mathrm{atom}}}
     K(\mathbf{X}_{i},\mathbf{X}_{B})
       @f]

     * storage order energyArray: first index types second index is reference configurations
     * @f$\mathbf{E}\in \mathbb{R}^{N_{types}\times N_{local\ ref\ confs}}@f$
     *******************************************************************************************/
    Vec2Real& energyArray;
    /*******************************************************************************************
     * Stores dimension of energy array which is number of local reference configurations.
     *******************************************************************************************/
    Vec1Size energyArrayDim;
    /*******************************************************************************************
     * Object helping to determine if an resize of the energy array is needed.
     *******************************************************************************************/
    ArrayResizing2D energyArraySize;
    /*******************************************************************************************
     * Array to store first part of the Kernel derivative w.r.t to descriptor SHS2/SHS3.
     *
     * for the polynomial Kernel this looks like
     *
     * @f[
        \sum_{B} \Lambda_{iB} \hat{\mathbf{X}}_{B}' - \sum_{B}
     w_{B}\Lambda'_{iB}\hat{\mathbf{X}}_{i}
       @f]
     * with
     * @f[\hat{\mathbf{X}}_{B}' = w_{B} \hat{\mathbf{X}}_{B}@f]
     * @f[\Lambda_{iB} = \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i} \cdot
     \hat{\mathbf{X}}_{B} \right)^{\xi-1}
       @f]
     * @f[\Lambda'_{iB}= \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i} \cdot
     \hat{\mathbf{X}}_{B} \right)^{\xi}@f]
     *
     * @note order of array is type index central atom and then squashed index of central atom,
     type1, type2, l, n0, n1
     * is computed in compute_derivativeMatrix
     * @note this matrix is used to compute the forcePreContract array
     *******************************************************************************************/
    std::map<String, ShVec2Real> derivativeMatrix;
    /*******************************************************************************************
     * Store dimenions of derivativeMatrix. Size is number atoms per type times number descriptors
     *******************************************************************************************/
    std::map<String, Vec1Size> derivativeMatrixDim;
    /*******************************************************************************************
     * Object to chek if the derivativeMatrix needs a resize.
     *******************************************************************************************/
    std::map<String, ArrayResizing2D> derivativeMatrixSize;
    /*******************************************************************************************
     * matrix storing the derivative of the kernel times the derivatives of the descriptors
     *
     * storage order is central atom as first index; second index neighbor type times
     * angular number times radial number x second times spherical harmonic parameter m.
     * As symbol @f$\mathbf{\tilde{L}}\in \mathbb{R}^{N_{atoms}\times lnm }@f$
     * @note the derivative of the descriptor with respect to the position is not included
     *******************************************************************************************/
    std::map<String, ShVec2Real> forcePreContract;
    /*******************************************************************************************
     * Store the dimensions of forcePreContract. Dimension is number SHS2 descriptors.
     *******************************************************************************************/
    std::map<String, Vec1Size> forcePreContractDim;
    /*******************************************************************************************
     * Object to determine if forcePreContract has to be resized.
     *******************************************************************************************/
    std::map<String, ArrayResizing2D> forcePreContractSize;
    /*******************************************************************************************
     * stores the pairwise forces obtained by
     *
     * array is computed in routine computeForceArrays( const DescriptorCollector&
     descriptorCollection )
     *
     * @f[
        \sum_{j\in
     \mathcal{N}_{i}}\sum_{B}w_{B}\frac{K(\mathbf{X}_{i},\mathbf{K}_{B})}{d\mathbf{r}_{j}}
       @f]
     *
     * forces are stored in order first index is central atom index;
     * second index is
     * 1st neighbor x
     * 1st neighbor y
     * 1st neighbor z
     * 2nd neighbor x
     * 2nd neighbor y
     * 2nd neighbor z
     * ...
     *******************************************************************************************/
    std::map<String, ShVec2Real> pairForces;
    /*******************************************************************************************
     * Store dimension of pairForces. Length of vectors is number of central atoms. Value stored
     * in vector is number of neighbors x 3
     *******************************************************************************************/
    std::map<String, Vec1Size> pairForcesDim;
    /*******************************************************************************************
     * Object to check if the pairForce array needs resize.
     *******************************************************************************************/
    std::map<String, ArrayResizing2D> pairForcesSize;
    /*******************************************************************************************
     * stores forces obatined by central derivative
     *
     * @f[
       -\sum_{B}w_{B}\frac{K(\mathbf{X}_{i},\mathbf{K}_{B})}{d\mathbf{r}_{i}}
       @f]
     *
     * ordering is atom1 x atom1 y atom1 z, atom2 x atom2 y atom2 z...
     *******************************************************************************************/
    std::map<String, ShVec1Real> centralForces;
    /*******************************************************************************************
     * Object to check if the central atom array has to be resized.
     *******************************************************************************************/
    std::map<String, ArrayResizing1D> centralForcesSize;
    /*******************************************************************************************
     * stress tensor obtained by  virial theorem
     *
     * @f[
          \tau_{\alpha,\beta}=\frac{1}{\Omega}\sum_{i=1}\sum_{j\in \mathcal{N}_{i}}\sum_{B}w_{B}
                              \frac{d\ K(\mathbf{X}_{i},\mathbf{K}_{B})}{dr_{j,\alpha}} *
     r^{ij}_{\beta} =
                              \frac{1}{\Omega}\sum_{i=1}\sum_{j\in \mathcal{N}_{i}}
                              F^{ij}_{\alpha}r^{ij}_{\beta}
       @f]
     *
     * stored in order xx, xy, xz, yx, yy, yz, zx, zy, zz
     *******************************************************************************************/
    std::map<String, ShVec1Real> stressTensor;
    /*******************************************************************************************
     * stress tensor obtained by  virial theorem and summed over various descriptors
     *
     * is sum of stressTensor variable over descriptors
     *
     * stored in order xx, xy, xz, yx, yy, yz, zx, zy, zz
     *******************************************************************************************/
    Vec1Real totalStressTensor;
    /*******************************************************************************************
     * Determines if totalStressTensor has to be resized. Will never happen except during
     * initialization.
     *******************************************************************************************/
    ArrayResizing1D totalStressTensorSize;
    /*******************************************************************************************
     * stores the total force acting on a single ion
     *
     * the order of the array is atom1 x, atom1 y, atom1 z, atom2 x, atom2 y, atom2 z,...
     *******************************************************************************************/
    Vec1Real atomicForces;
    /*******************************************************************************************
     * Object to determine if number of central atoms changed and atomicForce array has to be
     * resized.
     *******************************************************************************************/
    ArrayResizing1D atomicForcesSize;
    /*******************************************************************************************
     * check if force arrays pairForces and centralforces were computed.
     *
     * this is a necessary condition if the atomicForces array,
     * the stressTensor or the totalStressTensor are computed
     *******************************************************************************************/
    bool forceArraysComputed;
    /*******************************************************************************************
     * Temporary array needed for the computation of the derivativeMatrix.
     *******************************************************************************************/
    Vec2Real tempForceVector;
    /*******************************************************************************************
     * Object to determine if tempForceVector has to be resized.
     *******************************************************************************************/
    ArrayResizing2D tempForceVectorSize;
    /*******************************************************************************************
     * Stores dimensions of force vectors which is the number of atoms per type.
     *******************************************************************************************/
    Vec1Size tempForceVectorDim;
    /*******************************************************************************************
     * Look up table class to flatten loops for total stress tensor computation.
     *******************************************************************************************/
    StressTensorLookUpTable stressTensorLookUpTable;
    /*******************************************************************************************
     * Look up table for the atomic force computation.
     *******************************************************************************************/
    AtomicForcesLookUpTable atomicForcesLookUpTable;
    /*******************************************************************************************
     * Store inverse volume for stress computation.
     *******************************************************************************************/
    Real inverseVolume;
    /*******************************************************************************************
     * Entries for stress tensor computation. Are the components of the stress tensor which
     * allow to determine remaining components by symmetry.
     *******************************************************************************************/
    Vec1Size entries = {0, 1, 2, 4, 5, 8};
    /*******************************************************************************************
     * Make local copy of constants namespace descriptorKeyList to make them accesible in
     * parallel std algo regions
     *******************************************************************************************/
    const Vec1String descriptorKeyList = constants::descriptorKeyList;
    /*******************************************************************************************
     * Make local copy of constants namespace FUNIT to make them accesible in
     * parallel std algo regions. Needed for unit conversion of forces.
     *******************************************************************************************/
    const Real      FUNIT = constants::FUNIT;
    Vec1Size        centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapPredictor();

} //namespace vaspml

#endif
