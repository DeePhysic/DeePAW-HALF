#ifndef DESCRIPTORSHS3_HPP
#define DESCRIPTORSHS3_HPP

#include "ArrayResizing.hpp"
#include "Descriptor.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

namespace vaspml
{

class DescriptorSHS2;
class NearestNeighborNSquare;

class LookUpSparseMatrixspectrumVasp
{
  public:
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index for all index combinations
     *******************************************************************************************/
    Vec1Size centralIndex;
    /*******************************************************************************************
     * stores the descriptor indices of the two body descriptor from which the spectrumVasp is
     * computed
     *******************************************************************************************/
    Vec1Size nDesc;
    /*******************************************************************************************
     * Storing type indices of the central atom in the force field
     *******************************************************************************************/
    Vec1Size typeCentralFF;
    /*******************************************************************************************
     * number of atoms in the previous MD step. Needed to check if resize is needed
     *******************************************************************************************/
    Size nAtomsOld;
    /*******************************************************************************************
     * Count the number of elements in the flattened index arrays.
     *
     * @param nAtoms number of central atoms
     * @param descriptorList list of active descriptors.
     * @param typeMap map which assigns structure types to force field types and vice versa
     * @param typeIndexCentral type indices of central atoms in structure
     *******************************************************************************************/
    Size countElements(const Size&     nAtoms,
                       const Vec2Size& descriptorList,
                       const TypeMap&  typeMap,
                       const Vec1Int&  typeIndexCentral);
    /*******************************************************************************************
     * Count number of elements and update index arrays
     *
     * @param nAtoms number of central atoms.
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     * @param lookUpSize defines if look-up tables need a resize or a refill
     *******************************************************************************************/
    void update(const Size&     nAtoms,
                const Vec2Size& descriptorList,
                const TypeMap&  typeMap,
                const Vec1Int&  typeIndexCentral);
    /*******************************************************************************************
     * Filling all the idex arrays defined in this structure
     *
     * @param nAtoms number of central atoms.
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    void refill(const Size&     nAtoms,
                const Vec2Size& descriptorList,
                const TypeMap&  typeMap,
                const Vec1Int&  typeIndexCentral);
    LookUpSparseMatrixspectrumVasp();
};

class LookUpForcePreContract
{
  public:
    /*******************************************************************************************
     * Main index to loop over in std par version. Is a combination index over alll cetral atoms
     * and the corresponding descriptors.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * Central atom index
     *******************************************************************************************/
    Vec1Size centralIndex;
    /*******************************************************************************************
     * The central atom index in a type resolved form. Atoms are counted within its type.
     *******************************************************************************************/
    Vec1Size centralIndexPerType;
    /*******************************************************************************************
     * Type index of the central atom in the current structure.
     *******************************************************************************************/
    Vec1Size typeIndexCentral;
    /*******************************************************************************************
     * Type index of the central atom in the used machine learning force field.
     *******************************************************************************************/
    Vec1Size typeIndexCentralFF;
    /*******************************************************************************************
     * Index in the DescriptorSHS2 index of the current central atom and descriptor index.
     *******************************************************************************************/
    Vec2Size shs2Index;
    /*******************************************************************************************
     * Start index of the descriptor from which the contraction is started.
     *******************************************************************************************/
    Vec2Size shs2Start;
    /*******************************************************************************************
     * End index of the descriptor at which the contraction is terminated.
     *******************************************************************************************/
    Vec2Size shs2End;
    /*******************************************************************************************
     * Index in the derivative matrix which has to be multiplied by the with SHS2 during contraction
     *******************************************************************************************/
    Vec2Size indexDerivativeMatrix;
    /*******************************************************************************************
     * Stores the entries of the active descriptors for every force field type.
     *******************************************************************************************/
    Vec2Size unique_lnmCombiIndex;
    /*******************************************************************************************
     * The index of the descriptor in a combination of lnm. Does not contain central atom.
     *******************************************************************************************/
    Vec1Size lnmCombiCentral;
    /*******************************************************************************************
     * Array to store prefactors which come from analytical derivatives of kernel matrix
     *******************************************************************************************/
    Vec2Real preFactor;
    /*******************************************************************************************
     * Number of atoms in previous MD step. Needed to check if the index arrays have to be
     *recomputed
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Variable to check if the class was already initialized.
     *******************************************************************************************/
    bool initialized = false;
    /*******************************************************************************************
     * Initialize the class. Will only be called once when initialized = false.
     *
     * @param preFactorPreContract prefactors which are stored in preFactor for all combinations
     * @param sparseMap_derivativeMatrix map of activate descriptors in derivative matrix
     * @param sparseMap_SHS2 elements are active descriptors of 2-body DescriptorSHS2
     * @param sparseMap_central sparse map within the central atom.
     * @param typeMap gives mapping between force field types and structure types
     *******************************************************************************************/
    void init(const Vec2Real& preFactorPreContract,
              const Vec2Size& sparseMap_derivativeMatrix,
              const Vec2Size& sparseMap_SHS2,
              const Vec2Size& sparseMap_central,
              const TypeMap&  typeMap);
    /*******************************************************************************************
     * Count number of elements in the index arrays. Determine needed sizes for index arrays.
     *
     * @param nAtoms number of atoms in current structure.
     * @param typeIndexCentral type indexes of central atoms in structure
     * @param typeMap gives mapping between force field types and structure types
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Vec1Int& typeIndexCentral, const TypeMap& typeMap);
    /*******************************************************************************************
     * Update function for the whole object. If object was not initialized yet the function will
     * also do initialization.
     *
     * @param nAtoms number of atoms in current structure.
     * @param typeIndexCentral type index of central atom in current structure.
     * @param centralAtomIndexPerType index of central atom counted within the types.
     * @param sparseMap_central sparse map of central atoms
     * @param typeMap mapping between force field types and central atom types.
     *******************************************************************************************/
    void update(const Size&    nAtoms,
                const Vec1Int& typeIndexCentral,
                const Vec1Int& centralAtomIndexPerType,
                const TypeMap& typeMap);
    /*******************************************************************************************
     * Refilling the index arrays in case the number of central atoms changed.
     *
     * @param nAtoms number of atoms in current structure
     * @param typeIndexCentral type index of central atom in current structure.
     * @param centralAtomIndexPerType index of central atom counted within the types.
     * @param sparseMap_central sparse map of central atoms
     * @param typeMap mapping between force field types and central atom types.
     *******************************************************************************************/
    void refill(const Size&    nAtoms,
                const Vec1Int& typeIndexCentral,
                const Vec1Int& centralAtomIndexPerType,
                const TypeMap& typeMap);
};

class DescriptorSHS3 : public Descriptor
{
  public:
    /***************************************************************************************
     * @class DescriptorSHS3
     * @brief Three-body descriptor and it's derivative relevant stuff.
     *
     * @descriptorList list containing the used descriptors. This list is obtained from the
     * vasl ML_FF file
     * @angularFilterOn switching on and off angular filtering
     * @angularFilterType determine type of angular filtering function
     * @angularFilterScaling scaling parameter of angular filtering function
     * @param nRoots number of roots of modified spherical Bessel functions for various angular
     *parameter
     * @maxOrder maximal order of spherical Bessel functions/ spherical harmonics
     * @param spectrum_vasp optional shared_ptr where the SHS3 in vasp format is stored
     * @param spectrum_pair pair wise SHS3 spectrum where atom resolution is included
     *
     * This class is used to create a three-body descriptor. This descriptor
     * is derived from the two-body descriptor.
     * This class can for example be created with as part of a descriptor collector
     * DescriptorMap which is a std::map.
     * @code
     * DescriptorMap descriptors;
     * @endcode
     * It can then be accessed with the key "SHS3-3-body", but compared to the
     * regular three-body descriptor the shared pointer needs to be called
     * with the keyword "LinearScalingDescriptor".
     * The filling would look like the following (in the Frame class):
     * @code
     * descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3>(
     *    inputParameters["SHS3-3-body-descriptor-list"].dcget<ShVec2Int>(),
     *    inputParameters["SHS3-3-body-angular-filter-on"].cget<bool>(),
     *    inputParameters["SHS3-3-body-angular-filter-type"].cget<Int>(),
     *    inputParameters["SHS3-3-body-angular-filter-scale"].cget<Real>(),
     *    basis3Body.get_nValidRoots(),
     *    inputParameters["SHS3-3-body-max-angular-number"].cget<Int>(),
     *    inputParameters["SHS3-3-body-weight"].cget<Real>(),
     *    inputParameters["SHS3-3-body-is-normalized"].cget<bool>(),
     *    std::get<ShVec2Real>(arrayData.at("3-body-SHS3-spectrum")),
     *    std::get<ShVec2Real>(arrayData.at("3-body-SHS3-pair-spectrum")));
     * @endcode
     *
     ***************************************************************************************/
    DescriptorSHS3(const Vec2Int&       descriptorList,
                   const bool           angularFilterOn,
                   const Int            angularFilterType,
                   const Real           angularFilterScaling,
                   const Vec1Size&      nRoots,
                   const Int            maxOrder,
                   const Real           weight = 1.0,
                   const bool           isNormalized = false,
                   const DescriptorType descriptorType = DescriptorType::none,
                   ShRec                record = nullptr);

    DescriptorSHS3(const Size&          nTypes,
                   const bool           angularFilterOn,
                   const Int            angularFilterType,
                   const Real           angularFilterScaling,
                   const Vec1Size&      nRoots,
                   const Int            maxOrder,
                   const Real           weight = 1.0,
                   const bool           isNormalized = false,
                   const DescriptorType descriptorType = DescriptorType::none,
                   ShRec                record = nullptr);
    DescriptorSHS3(ShRec record = nullptr);
    /****************************************************************************************
     * Compute SHS3 (3-body spectrum ) for a given neighbor list configuration.
     *
     * @param shec spherical harmonics expansion coefficients for a given neighbor list
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *
     ****************************************************************************************/
    void computeSHS3(const std::shared_ptr<DescriptorSHS2>& shec, const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute SHS3 (3-body spectrum ) for a given neighbor list configuration by the use of
     * lookUp lists
     *******************************************************************************************/
    void computeSHS3GPULookUp();
    /*******************************************************************************************
     * Getter for 3-body descriptor (spectrum) in vasp format for a given central atom.
     *
     * @centralAtom central atom index which mactches neighbor list and expansion
     * coefficients shec
     *******************************************************************************************/
    const Vec1Real& get_SHS3Atom_vasp(const Int centralAtom) const;
    /*******************************************************************************************
     * Function returns same thing as get_SHS3Atom_vasp, but implements the virtual function of
     * Descriptor used in the DescriptorCollector class.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    const Vec1Real& get_descriptor(const Size centralAtom) const;
    /*******************************************************************************************
     * Returns the size of 3-body descriptor.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    Int get_sizeDescriptor(const Size centralAtom) const;
    /*******************************************************************************************
     * @brief Gets the number of elements in the MLFF (Machine Learning Force Field).
     *
     * @return The number of elements in the MLFF as a size type.
     *******************************************************************************************/
    Size get_numberElementsMLFF() const;
    /*******************************************************************************************
     * @brief Gets the number of elements in the structure.
     *
     * @return The number of elements in the structure as a size type.
     *******************************************************************************************/
    Size get_numberElementsStructure() const;
    /*******************************************************************************************
     * @brief Checks whether the angular filter is enabled.
     *
     * @return True if the angular filter is enabled, false otherwise.
     *******************************************************************************************/
    bool get_angularFilterOn() const;
    /*******************************************************************************************
     * @brief Retrieves the type of the angular filter.
     *
     * @return The angular filter type as an integer.
     *******************************************************************************************/
    Int get_angularFilterType() const;
    /*******************************************************************************************
     * @brief Retrieves the scaling factor for the angular filter.
     *
     * @return The scaling factor for the angular filter as a real number.
     *******************************************************************************************/
    Real get_angularFilterScaling() const;
    /*******************************************************************************************
     * @brief Retrieves the VASP spectrum for the entire system.
     *
     * @return A reference to a 2D vector of real numbers representing the VASP spectrum.
     *******************************************************************************************/
    const Vec2Real& get_spectrumVasp() const;
    /*******************************************************************************************
     * @brief Retrieves the VASP spectrum for a specific atom.
     *
     * @param atomIdx The index of the atom for which the VASP spectrum is requested.
     * @return A reference to a 1D vector of real numbers representing the VASP spectrum for the
     *specified atom.
     *******************************************************************************************/
    const Vec1Real& get_spectrumVasp(const Size atomIdx) const;
    Int             get_totalNumberBasisFunctions() const;
    /*******************************************************************************************
     * Rescale 3-body descriptor (spectrum).
     *
     * @param centralAtom central atom index for which descriptor is rescaled
     * @param centralAtom central atom index for which descriptor is rescaled
     *******************************************************************************************/
    void rescale_descriptor(const Size centralAtom, const Real scaleFactor);
    /*******************************************************************************************
     * Compute the derivative of a kernel matrix times the derivatives of the descriptors of this
     *class.
     *
     * @param derivativeMatrix kernel matrix derivative
     * @param forcePreContract on output contains derivativeMatrix times descriptor derivative
     * @param typeMap translates atom types from current structure to types stored in force field
     *
     * @f[ -\sum\limits_{i}\mathbf{L}_{i} \frac{d\mathbf{X}_{i}}{dp_{nn'l}^{iJJ'}}
     *      \frac{dp_{nn'l}^{iJJ'}}{dc_{n''lm}^{iJ''}}@f]
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Return the size of the second dimension (number of Cnlm) of the array forcePreContract.
     *
     * @param atomIndex index of central atom for which memory allocation is done
     *
     * This routine is used in the Predictor class for allocation (resize) of the
     * forcePreContract array. The size of forcePreContract is
     * [central atom ][ neighborType x l x nRadial x m ].
     *******************************************************************************************/
    Size get_forcePreContractSize(const Size atomIndex) const;
    /*******************************************************************************************
     * Compute centralForces and pairForces which are needed to compute total force on ion.
     *
     * @param forcePreContract array which stores the derivative of the
     *        kernelMatrix times the descriptor derivative
     * @param pairForces on output storing pair forces in same format as neighbor list
     * @param centralForces force acting on central atom by making derivative wrt to atom itself
     *
     * Routine is called from Predictor class to be able to predict total force
     * in predictor class. Routine implements a wrapper to the DescriptorSHS2 class
     * which implements the functionality.
     * @note For more implementation details check
     * compute_forcePreContractSparse. Currently
     * only compute_forcePreContractSparse is used.
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Setting up look up tables for parallel standard algorithms.
     *
     * @param typeMap mapping between force field types and structure types
     *******************************************************************************************/
    void make_lookUpTables(const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute the Look up tables needed for computeSHS3GPULookUp
     *
     * @param typeMap conversion map from structure types to force field types and vice versa
     *******************************************************************************************/
    void make_lookUpSpectrumVasp(const TypeMap& typeMap);
    /*******************************************************************************************
     *******************************************************************************************/
    void make_lookUpForcePreContract(const TypeMap& typeMap);
    /*******************************************************************************************
     * Setting parameters needed for the computation of spectrumVasp
     *
     * @param descriptorSHS2 two body descriptor from which the power spectrum is computed
     * @param typeMap conversion map from structure types to force field types
     *******************************************************************************************/
    void set_parameters(const std::shared_ptr<DescriptorSHS2>& descriptorSHS2,
                        const TypeMap&                         typeMap);
    /*******************************************************************************************
     * write spectrumVasp to some file in a single column format.
     *
     * @param fname file name to which spectrum will be written
     *******************************************************************************************/
    void write_spectrumVasp(const String& fname) const;
    /*******************************************************************************************
     * @brief Sets the parameters for the DescriptorSHS3 object.
     *
     * This method initializes the internal parameters of the DescriptorSHS3 object
     * with the provided values. These parameters control the behavior and configuration
     * of the descriptor, including its weight, normalization, type, execution policy,
     * and angular filtering settings.
     *
     * @param weight The weight parameter used in the descriptor calculations.
     * @param isNormalized A boolean flag indicating whether the descriptor should be normalized.
     * @param descriptorType The type of descriptor to be used (e.g., SHS3, MLFF, etc.).
     * @param numberElementsMLFF The number of elements used in the MLFF (Machine Learning Force
     *Field) calculations.
     * @param numberElementsStructure The number of elements used in the structural calculations.
     * @param angularFilterOn A boolean flag indicating whether angular filtering is enabled.
     * @param angularFilterType The type of angular filter to be applied (e.g., Gaussian, cosine,
     *etc.).
     * @param angularFilterScaling The scaling factor for the angular filter.
     *******************************************************************************************/
    void setParameters(const Real           weight,
                       const bool           isNormalized,
                       const DescriptorType descriptorType,
                       const Size           numberElementsMLFF,
                       const Size           numberElementsStructure,
                       const bool           angularFilterOn,
                       const Int            angularFilterType,
                       const Real           angularFilterScaling);
    /*******************************************************************************************
     * @brief Sets the spectrum data for the DescriptorSHS3 object.
     *
     * This function assigns the provided spectrum data to the `spectrumVasp` member variable
     * of the DescriptorSHS3 object. The `spectrumVasp` represents the spectrum data obtained
     * from VASP (Vienna Ab initio Simulation Package) or similar sources.
     *
     * @param spectrumVasp A reference to a `Vec2Real` object containing the spectrum data.
     *
     * @note The `spectrumPair` assignment is commented out and not currently used in this function.
     *******************************************************************************************/
    void setData(const Vec2Real& spectrumVasp);
    /*******************************************************************************************
     * @brief Extends the internal `spectrumVasp` data by appending the contents of the given
     *`spectrumVasp` vector.
     *
     * This function takes a vector of real values (`spectrumVasp`) and appends its elements to the
     *internal `spectrumVasp` member of the `DescriptorSHS3` class. The insertion is performed using
     *iterators to ensure efficient copying of the data.
     *
     * @param spectrumVasp A constant reference to a vector of real values (`Vec2Real`) that
     *contains the data to be appended to the internal `spectrumVasp` member.
     *
     * @note The commented-out code for `spectrumPair` suggests that similar functionality might be
     *intended for another member variable, but it is currently disabled.
     *******************************************************************************************/
    void extendData(const Vec2Real& spectrumVasp);
    /*******************************************************************************************
     * @brief Adds a new element to the DescriptorSHS3 object.
     *
     * This function appends a given spectrum vector to the `spectrumVasp` member of the
     *DescriptorSHS3 object. The `spectrumVasp` is a container that stores vectors of real numbers,
     *representing spectral data.
     *
     * @param spectrumVasp A vector of real numbers representing the spectral data to be added.
     *
     * @note The `spectrumPair` functionality is currently commented out and not used in this
     *implementation.
     *******************************************************************************************/
    void addElement(const Vec1Real& spectrumVasp);
    /*******************************************************************************************
     * @brief Retrieves the descriptor list.
     *
     * This function returns a constant reference to the `descriptorList` stored in the underlying
     *data.
     *
     * @return A constant reference to a `Vec2Int` object representing the descriptor list.
     *******************************************************************************************/
    const Vec2Int& get_descriptorList() const;
    void           compute_RefitDescDerivatives(const Int                     kAtom,
                                                const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Computes the derivative of the power spectrum (dPS) for a given atom.
     *
     * This function calculates the derivative of the power spectrum (dPS) for a specified atom
     * based on its nearest neighbors. The computation involves iterating over the neighbors
     * and their associated types, angular momentum orders, and other parameters to compute
     * the contributions to the dPS. The results are stored in the `dPS` member variable.
     *
     * @param kAtom The index of the central atom for which the dPS is computed.
     * @param maxNeighborList The list of nearest neighbors for the central atom, including
     *                        their types and global indices.
     *
     * The computation involves:
     * - Allocating memory for the dPS array.
     * - Accessing the second-order descriptor (DescriptorSHS2) for required data.
     * - Iterating over neighbors and their types to compute contributions to the dPS.
     * - Handling both the central atom's contributions and its neighbors' contributions.
     *
     * @note This function assumes that the `descriptorSHS2` member is already initialized
     *       and contains the necessary data for computation.
     * @note The `dPS` array is updated in-place with the computed values.
     *******************************************************************************************/
    void compute_dPS(const Int kAtom, const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Computes the derivative of the power spectrum head (dPSHead) for a given atom.
     *
     * This function calculates the derivative of the power spectrum head for a specific atom
     * based on its type and the associated neighbor list. The computation involves iterating
     * over the angular and radial components of the descriptor and applying the weight factors.
     * The results are stored in the `dPS` data structure.
     *
     * @param kAtom The index of the central atom for which the derivative of the power spectrum
     *head is computed.
     *
     * The computation is divided into two parts:
     * 1. **Power Spectrum Calculation**:
     *    - Iterates over the type pairs and angular components.
     *    - Computes the power spectrum using the `clnmVasp` values.
     *    - Stores the result in `dPS[kType][nDesc]`.
     *
     * 2. **Derivative Calculation**:
     *    - Iterates over the Cartesian directions (x, y, z).
     *    - Computes the derivative terms using `clnmVasp` and `clnmVaspDerivative`.
     *    - Stores the result in `dPS[kType][nDesc + shift]`.
     *
     * Dependencies:
     * - Neighbor list for retrieving the type index of the central atom.
     * - DescriptorSHS2 for accessing `clnmVasp` and `clnmVaspDerivative`.
     * - Type mappings, angular lists, and weight factors for the computation.
     *
     * Notes:
     * - The function assumes that `allocate_dPSHead()` has been called to allocate memory for
     *`dPS`.
     * - The computation is parallelized using `std::for_each` for efficiency.
     *******************************************************************************************/
    void            compute_dPSHead(const Int kAtom);
    const Vec2Real& get_refitHeadTerm() const;

  private:
    /*******************************************************************************************
     * @brief Allocates and initializes the `dPS` data structure based on the maximum number of
     *neighbors.
     *
     * This function resizes the `dPS` vector to accommodate the maximum number of neighbors in the
     *provided `maxNeighborList`. Each slice of `dPS` is further resized to match the required
     *dimensions based on the number of types and the size of `type0List[0]`. The values in `dPS`
     *are initialized to zero.
     *
     * @param maxNeighborList A reference to the `NearestNeighborNSquare` object that provides the
     *maximum neighbor list size.
     *
     * @note This function only supports CPU single-core execution. GPU execution is not implemented
     *and will throw an error.
     *******************************************************************************************/
    void allocate_dPS(const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Allocates and initializes the `dPS` data structure for the head of the descriptor.
     *
     * This function resizes the `dPS` vector to match the number of types in the `neighborList`.
     *Each slice of `dPS` is further resized to match the required dimensions based on the number of
     *types and the size of `type0List[0]`. The values in `dPS` are initialized to zero.
     *
     * @note This function only supports CPU single-core execution. GPU execution is not implemented
     *and will throw an error.
     *******************************************************************************************/
    void allocate_dPSHead();
    /*******************************************************************************************
     * Loop over atoms for computeSHS3SingleAtom.
     *
     * Single CPU version, using a for loop.
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *******************************************************************************************/
    void computeSHS3CPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * Loop over atoms for computeSHS3SingleAtom.
     *
     * GPU version, using standard library algorithms.
     *
     * @shec shared pointer to two-body descriptor
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *******************************************************************************************/
    void computeSHS3GPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor times a scalar.
     *
     * Single CPU version, uses c++ standard algorithms for computations.
     *
     * @param derivativeMatrix stores derivative of kernel
     *matrix.@f[L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f]
     * @param forcePreContract stores the result. The product of the derivativeMatrix and 2 body
     *descriptor
     * @param typeMap gives the ability to convert between structure and force field types
     *
     * @f[
     *f^{iJ_{1}}_{ln_{1}m}=\sum_{J_{2}}\sum_{n_{2}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     *n_{2}m} @f]
     * @f[
     *\tilde{f}^{iJ_{2}}_{ln_{2}m}\sum_{J_{1}}\sum_{n_{1}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l
     *n_{1}m} @f] these terms are then summed up.
     * @f[ F^{iJ}_{lnm}=f^{iJ}_{lnm}+f^{iJ}_{lnm}@f]
     *
     * @note The algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses an index list representation to compute both terms in one loop.
     *******************************************************************************************/
    void compute_forcePreContractSparseSingleCPU(const Vec2Real& derivativeMatrix,
                                                 Vec2Real&       forcePreContract,
                                                 const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor times a scalar.
     *
     * GPU version, uses c++ standard algorithms for computations.
     *
     * @param derivativeMatrix stores derivative of kernel
     *matrix.@f[L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f]
     * @param forcePreContract stores the result. The product of the derivativeMatrix and 2 body
     *descriptor
     * @param typeMap gives the ability to convert between structure and force field types
     *
     * @f[
     *f^{iJ_{1}}_{ln_{1}m}=\sum_{J_{2}}\sum_{n_{2}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     *n_{2}m} @f]
     * @f[
     *\tilde{f}^{iJ_{2}}_{ln_{2}m}\sum_{J_{1}}\sum_{n_{1}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l
     *n_{1}m} @f] these terms are then summed up.
     * @f[ F^{iJ}_{lnm}=f^{iJ}_{lnm}+f^{iJ}_{lnm}@f]
     *
     * @note The algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses an index list representation to compute both terms in one loop.
     *******************************************************************************************/
    void compute_forcePreContractSparse_stdpar(const Vec2Real& derivativeMatrix,
                                               Vec2Real&       forcePreContract) const;
    /*******************************************************************************************
     * Compute index lists for compute_forcePreContract array. This is a necessary prerequisite
     * for calling compute_forcePreContractSparse.
     *
     * computes: \n
     * sparseMap_derivativeMatrix \n
     * sparseMap_SHS2 \n
     * sparseMap_central \n
     * preFactorPreContract \n
     * To implement the product of the matrix @f$L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$ the
     * following selection rules have to be obeyed:
     * @f[ \frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}=
     *    a_{l}\left[ \delta_{n_{1}n_{3}} \delta_{J_{1}J_{3}}c^{iJ_{2}}_{ln_    {2}m_{1}} +
     *    \delta_{n_{2}n_{3}} \delta_{J_{2}J_{3}}c^{iJ_{1}}_{ln_{1}m} \right] @f]
     *
     * The lists are implemented as sparse lists for
     * @f$\frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     *******************************************************************************************/
    void prepareSparseforcePreContract(const TypeMap& typeMap);
    /*******************************************************************************************
     * Make list of the active descriptors used in the current force field.
     *
     * @param descriptorList list containing the used descriptors. This list is obtained from the
     * vasp ML_FF file
     * @param nRoots the zeros of the modified spherical Bessel functions as computed
     * in desriptorRadialSpline
     * @param maxOrder maximal order of the angular parameter of the spherical harmonics
     *******************************************************************************************/
    void make_sparseList(const Vec2Int& descriptorList, const Vec1Size& nRoots, const Int maxOrder);
    /*******************************************************************************************
     * Make list of the active descriptors used in the current force field.
     *
     * @param descriptorList list containing the used descriptors. This list is obtained from the
     * vasp ML_FF file
     * @param nRoots the zeros of the modified spherical Bessel functions as computed
     * in desriptorRadialSpline
     * @param maxOrder maximal order of the angular parameter of the spherical harmonics
     *******************************************************************************************/
    void make_sparseListInference(const Vec2Int&  descriptorList,
                                  const Vec1Size& nRoots,
                                  const Int       maxOrder);
    /*******************************************************************************************
     * @brief Constructs sparse lists for descriptors based on the given parameters.
     *
     * This function initializes and populates several lists (`type0List`, `type1List`, `type2List`,
     * `angularList`, `n0List`, `n1List`, and `weightFactor`) that are used to represent sparse
     * descriptors. The lists are populated based on the number of supported types (`nTypesSup`),
     * the number of roots (`nRoots`), and the maximum order (`maxOrder`). The function iterates
     * through combinations of types and orders to generate the descriptor data.
     *
     * @param nTypesSup The number of supported types. Determines the size of the lists.
     * @param nRoots A vector specifying the number of roots for each angular order.
     *               The size of this vector should correspond to the maximum order + 1.
     * @param maxOrder The maximum angular order to consider. Determines the range of angular
     *orders.
     *
     * @details
     * - The function clears and resizes the lists to match the number of supported types
     *(`nTypesSup`).
     * - For each combination of `type0`, `type1`, `type2`, `orderL`, `order0`, and `order1`,
     *   the lists are populated with corresponding values.
     * - The `weightFactor` list is calculated based on whether `order0` equals `order1`.
     *   If they are equal, the weight factor is scaled by `angularFilter[orderL]`.
     *   Otherwise, it is scaled by `constants::SQRT2 * angularFilter[orderL]`.
     * - The function uses nested loops to iterate through all possible combinations of types and
     *orders.
     *
     * @note
     * - The `angularFilter` array and `constants::SQRT2` are assumed to be defined elsewhere in the
     *code.
     * - The function assumes that `nRoots` provides valid values for each angular order up to
     *`maxOrder`.
     * - The lists (`type0List`, `type1List`, `type2List`, `angularList`, `n0List`, `n1List`,
     *`weightFactor`) are member variables of the `DescriptorSHS3` class.
     *******************************************************************************************/
    void make_sparseListRefit(const Size& nTypesSup, const Vec1Size& nRoots, const Int maxOrder);
    /*******************************************************************************************
     * Precompute anfgular filtering functions.
     *
     * @param maxOrder maximal order of angular parameter of spherical harmonics
     *
     * @note more information can be found https://www.vasp.at/wiki/index.php/ML_IAFILT2
     *******************************************************************************************/
    void make_angularFilter(const Int maxOrder);
    /*******************************************************************************************
     * Allocate storage arrays for spectrumVasp. Calls resizeArraysGPU when std parallelism is
     * used
     *
     *@param typeMap map class which associates machine-learning force field types with types in
     * current structure
     *
     *@param shec spherical harmonics expansion coefficients precomputed for current neighbor list
     * instance
     *******************************************************************************************/
    void allocateArrays(const TypeMap& typeMap);
    /*******************************************************************************************
     * Allocating storage for GPU code. The array sizes are checked by hand to avoid
     * unescessary device to host and host to device copies.
     *******************************************************************************************/
    void resizeArraysGPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute the 3-body spectrum (stored in spectrumVasp) for given central atom.
     *
     * @param shec  spherical harmonics expansion coefficients from which the 3-body spectrum is
     *computed
     * @param atom central atom index ordered as in neighhbor list or in shec
     *******************************************************************************************/
    void computeSHS3SingleAtom(const DescriptorSHS2& shec, const Size atom);
    /*******************************************************************************************
     * Compute the 3-body spectrum (stored in spectrumVasp) for given central atom by the use
     * of Look up tables and a single loop.
     *
     * @param centralType_ff force field type index in force field
     * @param nDesc descriptor index
     * @param clnmVasp pair coefficients from which the power spectrum is computed
     * @param spectrumVasp variable to which power spectrum is stored for given coefficients
     *******************************************************************************************/
    void computeSHS3SingleAtomLookUp(const Size&     centralType_ff,
                                     const Size&     nDesc,
                                     const Vec1Real& clnmVasp,
                                     Real&           spectrumVasp);
    /*******************************************************************************************
     * Number of elements which are present in force field.
     *******************************************************************************************/
    Size numberElementsMLFF;
    /*******************************************************************************************
     * Number of elements which are present in structure.
     *******************************************************************************************/
    Size numberElementsStructure;
    /*******************************************************************************************
     * Contains the SHS3 (3-body descriptor).
     *
     * First index is the central atom index.
     * Second index is a combined index of
     * type1, type2, l, n0, n1.
     *******************************************************************************************/
    Vec2Real& spectrumVasp;
    /*******************************************************************************************
     * Structure which helps to determine the sizes of the array spectrumVasp.
     *
     * This is used to avoid not needed host to device and device to host transfers when GPUs
     * are used.
     *******************************************************************************************/
    ArrayResizing2D spectrumVaspSize;
    /*******************************************************************************************
     * Pair wise storage of SHS3 (3-body descriptor).
     *
     * Not used at the moment.
     *******************************************************************************************/
    Vec2Real& spectrumPair;
    /*******************************************************************************************
     * Number of atoms for which power spectrum is stored.
     *******************************************************************************************/
    Int numberAtoms;
    /*******************************************************************************************
     * List of centram atom types.
     *
     * @Note The type lists here need a converter
     * from force field type to structure type.
     *******************************************************************************************/
    Vec2Size type0List;
    /*******************************************************************************************
     * List of neighbor types 1.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     *******************************************************************************************/
    Vec2Size type1List;
    /*******************************************************************************************
     * List of neighbor types 2.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     *******************************************************************************************/
    Vec2Size type2List;
    /*******************************************************************************************
     * List of angular number of spherical harmonic l.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     *******************************************************************************************/
    Vec2Size angularList;
    /*******************************************************************************************
     * Tist of first radial number ( number of roots ) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     *******************************************************************************************/
    Vec2Size n0List;
    /*******************************************************************************************
     * Tist of second radial number ( number of roots ) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     *******************************************************************************************/
    Vec2Size n1List;
    /*******************************************************************************************
     * Weight factors to compute n0List n1List terms.
     *
     * The weight factor will be 1 for n0List == n1List
     * and will be sqrt(2) otherwise. This alllows to compute
     * only the upper triangular matrix of n0 n1 terms.
     *
     *******************************************************************************************/
    Vec2Real weightFactor;
    /*******************************************************************************************
     * Controls whether the angular filtering is switched on or not.
     *
     * For more information check the vasp wiki
     * https://www.vasp.at/wiki/index.php/ML_IAFILT2
     *
     *******************************************************************************************/
    bool angularFilterOn;
    /*******************************************************************************************
     * Determine which type of angular filtering is used.
     *******************************************************************************************/
    Int angularFilterType;
    /*******************************************************************************************
     * Drray stores the precomputed angular filter.
     *
     * Dimension is maxOrder + 1, filter value for every l is stored.
     *******************************************************************************************/
    Vec1Real angularFilter;
    /*******************************************************************************************
     * Scaling factor to determine the width
     * of the distribution with which the angular filtering is done.
     *******************************************************************************************/
    Real angularFilterScaling;
    /*******************************************************************************************
     * Store SHS2 descriptor which was used to build the SHS3 descriptor.
     *******************************************************************************************/
    std::shared_ptr<const DescriptorSHS2> descriptorSHS2;
    /*******************************************************************************************
     * Check if the sparse maps for sparsePreContract were set up.
     *******************************************************************************************/
    bool isSparseforcePreContractReady;
    /*******************************************************************************************
    * Sparse map which shows the entries in the derivativeMatrix which
    * are needed for the computation of the array forcePreContract.
    *
    * The derivativeMatrix is given by
    * @f[ L^{iJ_{1}J_{2}}_{ln_{1}n_{2}} =
          \frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}
          {\partial \hat{\mathbf{X}}_{i}} @f]
    *@f[ i @f] central atom index
    *@f[J_{1} @f] index denoting first neighbor type
    *@f[J_{2} @f] index denoting second neighbor type
    *@f[ l @f] indexing angular parameter
    *@f[ n_{1} @f] radial index belonging to neighbor type 1
    *@f[ n_{2} @f] radial index belonging to neighbor type 2
    * The order in which the indices are described denotes storage order.
    * The sparse list stores those elements of the derivativeMatrix for which the equation is
    non-zero.
    *
    * The actual derivative matrix is not stored in this class, but for example in Predictor.hpp.
    *******************************************************************************************/
    Vec2Size sparseMap_derivativeMatrix;
    /*******************************************************************************************
     * Sparse map which shows the entries in the SHS2-descriptor which
     * are needed for the computation of the array forcePreContract.
     *
     * forcePreContract is computed by the following equation. Einstein sum convention is used
     * over @f[J_{2}@f] and @f[ n_{2}@f] in first term. And in the second term the sums are taken
     over
     * @f[J_{1}@f] and @f[n_{2}@f]. The central atom index i is summed over in both terms
     * @f[
       a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     n_{2}m}+a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l n_{1}m}
       @f]
     *@f[ L^{iJ_{1}J_{2}}_{ln_{1}n_{2}} @f] derivativeMatrix
     *@f[ c^{iJ_{2}}_{l n_{2}m} @f] descriptorSHS2
     *@f[ a_{l} @f] preFactorPreContract
    *******************************************************************************************/
    Vec2Size sparseMap_SHS2;
    /*******************************************************************************************
    * Sparse map which shows the entries in the forcePreContract array (central atom index)
    * which are needed for the computation of the array forcePreContract
    *
    * forcePreContract is computed by the following equation. Einstein sum convention is used
    * over @f[J_{2}@f] and @f[ n_{2}@f] in first term. And in the second term the sums are taken
    over
    * @f[J_{1}@f] and @f[n_{2}@f]
    * @f[
      a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
    n_{2}m}+a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l n_{1}m}
      @f]
    * The result will therefore be a matrix @f[ Y^{kJ}_{lnm}@f] this map stores to which entry of
    * Y the result is written.
    *******************************************************************************************/
    Vec2Size sparseMap_central;
    /*******************************************************************************************
    * Pre-factors for the forcePreContract array computation
    * which are needed for the computation of the array forcePreContract.
    *
    * The prefactors can be expressed with Kronecker-Deltas:
    * @f[
      f^{iJ_{1}J_{2}}_{ln_{1}n_{2}} =
      \begin{cases}
      2 * a_{l} &\text{if } J_{1}=J_{2} \land n_{1}=n_{2} \\
      a_{l} &\text{if } J_{1}\neq J_{2} \land n_{1}=n_{2}\\
      \sqrt{2}a_{l} &\text{else}
      \end{cases}
    @f]
    * the factor @f[a_{l}@f] is stored in angularFilter.
    *******************************************************************************************/
    Vec2Real preFactorPreContract;
    /*******************************************************************************************
     * Stores indices needed for computation of spectrumVasp
     *******************************************************************************************/
    LookUpSparseMatrixspectrumVasp lookUpSpectrumVasp;
    /*******************************************************************************************
     * Look up table to compute force pre contraction
     *******************************************************************************************/
    LookUpForcePreContract lookUpForcePreContract;
    /*******************************************************************************************
     * Object to check if number of atoms has changed
     *******************************************************************************************/
    ArrayResizing1D centralAtomIndexSize;
    /*******************************************************************************************
     * Store size of descriptor, which is equal to the number of active descriptors in the sparse
     * list of n0Lis. First index is type index in force field and second intex is descriptor index
     *******************************************************************************************/
    Vec1Size descriptorSize;
    /*******************************************************************************************
     * local type map used to avoid copies between cpu and gpu of object
     *******************************************************************************************/
    TypeMap typeMapLoc;

    Vec2Real dPS;
    Vec2Real dPSHead;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapDescriptorSHS3();

/*******************************************************************************************
 * @brief Copies the neighbor list from one DescriptorSHS3 object to another.
 *
 * This function copies the neighbor list from the source descriptor (`descA`)
 * to the target descriptor (`descB`) by using the `get_neighborList_ptr()`
 * method of `descA` and setting it in `descB` using `set_neighborList()`.
 *
 * @param[in] descA The source DescriptorSHS3 object from which the neighbor list is copied.
 * @param[out] descB The target DescriptorSHS3 object to which the neighbor list is copied.
 *******************************************************************************************/
void copyNeighborList(const DescriptorSHS3& descA, DescriptorSHS3& descB);
/*******************************************************************************************
 * @brief Copies various parameters from one DescriptorSHS3 object to another.
 *
 * This function copies multiple parameters from the source descriptor (`descA`)
 * to the target descriptor (`descB`) using the respective getter methods of `descA`
 * and the `setParameters()` method of `descB`.
 *
 * @param[in] descA The source DescriptorSHS3 object from which the parameters are copied.
 * @param[out] descB The target DescriptorSHS3 object to which the parameters are copied.
 *******************************************************************************************/
void copyParameters(const DescriptorSHS3& descA, DescriptorSHS3& descB);
/*******************************************************************************************
 * @brief Copies spectrum data from one DescriptorSHS3 object to another.
 *
 * This function copies the spectrum data from the source descriptor (`descA`)
 * to the target descriptor (`descB`) using the `get_spectrumVasp()` method of `descA`
 * and the `setData()` method of `descB`.
 *
 * @param[in] descA The source DescriptorSHS3 object from which the spectrum data is copied.
 * @param[out] descB The target DescriptorSHS3 object to which the spectrum data is copied.
 *******************************************************************************************/
void copyData(const DescriptorSHS3& descA, DescriptorSHS3& descB);
/*******************************************************************************************
 * @brief Extends the spectrum data of one DescriptorSHS3 object using another.
 *
 * This function extends the spectrum data of the target descriptor (`descB`)
 * by appending the spectrum data from the source descriptor (`descA`) using
 * the `get_spectrumVasp()` method of `descA` and the `extendData()` method of `descB`.
 *
 * @param[in] descA The source DescriptorSHS3 object whose spectrum data is used for extension.
 * @param[out] descB The target DescriptorSHS3 object whose spectrum data is extended.
 *******************************************************************************************/
void extendData(const DescriptorSHS3& descA, DescriptorSHS3& descB);
/*******************************************************************************************
 * @brief Adds element-specific spectrum data from one DescriptorSHS3 object to another.
 *
 * This function adds the spectrum data of a specific element (identified by `atomIdx`)
 * from the source descriptor (`descA`) to the target descriptor (`descB`) using
 * the `get_spectrumVasp(atomIdx)` method of `descA` and the `addElement()` method of `descB`.
 *
 * @param[in] descA The source DescriptorSHS3 object from which the element-specific spectrum data
 *is copied.
 * @param[out] descB The target DescriptorSHS3 object to which the element-specific spectrum data is
 *added.
 *******************************************************************************************/
void addElementData(const DescriptorSHS3& descA, DescriptorSHS3& descB, const Size atomIdx);

} // namespace vaspml

#endif
