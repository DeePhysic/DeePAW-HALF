#ifndef DESCRIPTORSHS3REDUCEDLINELEM_HPP
#define DESCRIPTORSHS3REDUCEDLINELEM_HPP

#include "ArrayResizing.hpp"
#include "Descriptor.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

namespace vaspml
{

class DescriptorSHS2;

class computeSHS3clnmNoElementGPULookUp
{
  public:
    Vec1Size mainIndex;
    Vec1Size basisSetIndex;
    Vec1Size atomIndex;
    Vec1Size typeIndex;
};

class DescriptorSHS3ReducedLinElem : public Descriptor
{
  public:
    /***************************************************************************************
     * @class DescriptorSHS3ReducedLinElem
     * @brief Reduced (linear-scaling) three-body descriptor and it's derivative relevant stuff.
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
     * This class is used to create a reduced three-body descriptor. This descriptor
     * is derived from the two-body descriptor.
     * This class can for example be created with as part of a descriptor collector
     * DescriptorMap which is a std::map.
     * @code
     * DescriptorMap descriptors;
     * @endcode
     * It can then be accessed with the key "SHS3-3-body".
     * The filling would look like the following (in the Frame class):
     * @code
     * descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3ReducedLinElem>(
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
     ****************************************************************************************/
    DescriptorSHS3ReducedLinElem(const Vec2Int&       descriptorList,
                                 const bool           angularFilterOn,
                                 const Int            angularFilterType,
                                 const Real           angularFilterScaling,
                                 const Vec1Size&      nRoots,
                                 const Int            maxOrder,
                                 const Real           weight = 1.0,
                                 const bool           isNormalized = false,
                                 const DescriptorType descriptorType = DescriptorType::none,
                                 ShRec                record = nullptr);
    DescriptorSHS3ReducedLinElem(const Size&          nTypes,
                                 const bool           angularFilterOn,
                                 const Int            angularFilterType,
                                 const Real           angularFilterScaling,
                                 const Vec1Size&      nRoots,
                                 const Int            maxOrder,
                                 const Real           weight = 1.0,
                                 const bool           isNormalized = false,
                                 const DescriptorType descriptorType = DescriptorType::none,
                                 ShRec                record = nullptr);
    DescriptorSHS3ReducedLinElem(ShRec record = nullptr);
    /***************************************************************************************
     * Compute reduced SHS3 (3-body spectrum ) for a given neighbor list configuration.
     *
     * @param shec spherical harmonics expansion coefficients for a given neighbor list
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     ****************************************************************************************/
    void computeSHS3(const std::shared_ptr<DescriptorSHS2>& shec, const TypeMap& typeMap);
    /***************************************************************************************
     * Get reduced 3-body descriptor computed from the elements of the DescriptorSHS2 class.
     *
     * The returned descriptor is already contracted over atoms of the same type.
     *
     * Function implements the virtual function of
     * Descriptor used in the DescriptorCollector class.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     ****************************************************************************************/
    const Vec1Real& get_descriptor(const Size centralAtom) const;
    /***************************************************************************************
     * Returns the size of reduced 3-body descriptor.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     ****************************************************************************************/
    Int get_sizeDescriptor(const Size centralAtom) const;
    /***************************************************************************************
     * Rescale descriptor.
     *
     * @param centralAtom central atom index for which descriptor is rescaled
     * @param centralAtom central atom index for which descriptor is rescaled
     ****************************************************************************************/
    void rescale_descriptor(const Size centralAtom, const Real scaleFactor);
    /***************************************************************************************
     * Get index within the element agnostic Cnlm (intermediate descriptor).
     *
     * @param l is the angular quantum number
     * @param lOffset offset for l quantum number
     * @param n is the radial number
     * @param m is the magnetic quantum number
     ****************************************************************************************/
    Size getIndexNoElement(const Size l, const Size n, const Size m) const;
    /***************************************************************************************
     * Get vector of element-agnostic Cnlm (intermediate descriptor) for a certain central_alom).
     *
     * @param central_atom central atom for Cnlm
     * @param l is the angular quantum number
     * @param lOffset offset for l quantum number
     * @param n is the radial number
     * @param m is the magnetic quantum number
     ****************************************************************************************/
    const Vec1Real& get_clnmVaspNoElement(const Size central_atom) const;
    Int             get_totalNumberBasisFunctions() const;
    /*******************************************************************************************
     * Computes the derivative of a kernel matrix times the derivative of the descriptor of this
     *class.
     *
     * @param derivativeMatrix kernel matrix derivative
     * @param forcePreContract on output contains derivativeMatrix times descriptor derivative
     * @param typeMap translates atom types from current structure to types stored in force field
     *
     * @f[ -\sum\limits_{i}\mathbf{L}_{i} \frac{d\mathbf{X}_{i}}{dp_{nn'l}^{iJ}}
     *      \frac{dp_{nn'l}^{iJ}}{dc_{n''lm}^{iJ''}} @f]
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Return the size of the second dimension of the array forcePreContract.
     *
     * @param atomIndex index of central atom for which memory allocation is done
     *
     * This routine is used in the Predictor class for allocation of the
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
     * which implements the functionallity.
     * @note For more implementation details check
     * compute_forcePreContractSparse. Currently
     * only compute_forcePreContractSparse is used.
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Getter for three-body descriptor (spectrum) in vasp format for a given central atom.
     *
     * @centralAtom central atom index which mactches neighbor list and expansion
     * coefficients shec
     *******************************************************************************************/
    const Vec1Real& get_SHS3Atom_vasp(const Int centralAtom) const;

  private:
    /***************************************************************************************
     * compute type linear SHS3 descriptor for every atom in simulation box
     *
     * non parallel CPU version of computeSHS3
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     ****************************************************************************************/
    void computeSHS3CPU(const TypeMap& typeMap);

    /***************************************************************************************
     * compute type linear SHS3 descriptor for every atom in simulation box
     *
     * GPU version computeSHS3 using the parallel std library
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     ****************************************************************************************/
    void computeSHS3GPU(const TypeMap& typeMap);

    /***************************************************************************************
     * Make list of the active descriptors used in the current force field.
     *
     * @descriptorList list containing the used descriptors. This list is obtained from the
     * vasp ML_FF file
     * @nRoots the zeros of the modified spherical Bessel functions as computed
     * in desriptorRadialSpline
     * @maxOrder maximal order of the angular parameter of the spherical harmonics
     ****************************************************************************************/
    void make_sparseList(const Vec2Int& descriptorList, const Vec1Size& nRoots, const Int maxOrder);
    /*******************************************************************************************
     * @brief Constructs and refits sparse descriptor lists for reduced linear elements.
     *
     * This method initializes and populates various descriptor lists based on the
     * provided number of supported types, root counts, and maximum angular order.
     * It generates combinations of types, angular orders, and root indices, while
     * applying weight factors and filters for descriptor calculations.
     *
     * @param nTypesSup Number of supported types.
     * @param nRoots Vector containing the number of roots for each angular order.
     * @param maxOrder Maximum angular order to consider.
     *
     * - Clears and resizes internal lists (`type0List`, `type1List`, `type2List`,
     *   `angularList`, `n0List`, `n1List`, `weightFactor`) to match the number of
     *   supported types.
     * - Iterates over combinations of types, angular orders, and root indices to
     *   populate the lists.
     * - Applies weight factors based on angular filters and root indices.
     * - Updates the `descriptorList` in the associated data object for external use.
     *******************************************************************************************/
    void make_sparseListRefit(const Size& nTypesSup, const Vec1Size& nRoots, const Int maxOrder);
    /***************************************************************************************
     * Precompute angular filtering functions.
     *
     * @param maxOrder maximal order of angular parameter of spherical harmonics
     *
     * @note More information can be found https://www.vasp.at/wiki/index.php/ML_IAFILT2
     ****************************************************************************************/
    void make_angularFilter(const Int maxOrder);
    /***************************************************************************************
     * Allocating storage arrays.
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *
     * @shec spherical harmonics expansion coefficients precomputed for current neighbor list
     *instance
     ****************************************************************************************/
    void allocateArrays(const TypeMap& typeMap);
    void resizeArraysGPU(const TypeMap& typeMap);
    /***************************************************************************************
     * Compute element agnostic Clnm (intermediate descriptors) by summing over elements.
     * The Clnm are stored in clnmVaspNoElement.
     ****************************************************************************************/
    void computeSHS3clnmNoElement();
    void computeSHS3clnmNoElementGPU();
    /***************************************************************************************
     * Compute the reduced three-body spectrum (SHS3) for a given central atom.
     *
     * @shec  spherical harmonics expansion coefficients from which the 3-body spectrum is computed
     * @map class which associates machine-learning force field types with types in
     * current structure
     * @atom central atom index ordered as in neighhbor list or in shec
     ****************************************************************************************/
    void computeSHS3SingleAtom(const TypeMap& typeMap, const Size atom);
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor
     *
     * Single CPU version, uses c++ standard algorithms without parallel execution policy
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
     * @note the algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses a index list representation to compute both terms in one loop
     *******************************************************************************************/
    void compute_forcePreContractSparseCPU(const Vec2Real& derivativeMatrix,
                                           Vec2Real&       forcePreContract,
                                           const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor
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
     * @note the algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses a index list representation to compute both terms in one loop
     *******************************************************************************************/
    void compute_forcePreContractSparseGPU(const Vec2Real& derivativeMatrix,
                                           Vec2Real&       forcePreContract,
                                           const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Compute index lists for compute_forcePreContract array. This is a necessary prerequisite
     * for calling compute_forcePreContractSparse.
     *
     * computes: \n
     * sparseMap_derivativeMatrix \n
     * sparseMap_SHS2 \n
     * sparseMap_central \n
     * preFactorPreContract \n
     * This routine computes the product of the matrix @f$L^{iJ_{1}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     * The following derivative is implemented via sparse lists:
     * @f\frac{\mathrm{d}p^{iJ}_{nn'l}}{\mathrm{d}c^{iJ''}_{n''lm}} = \sqrt{\frac{8\pi^2}{2l+1}}
     * \sum_{m=-l}^{l} \left[\left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) + c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right) \right].@f]
     *******************************************************************************************/
    void prepareSparseforcePreContract(const TypeMap& typeMap);

    /************************************************************
     * Check if the sparse maps for sparsePreContract were set up.
     ************************************************************/
    bool isSparseforcePreContractReady;
    /***************************************************************************************
     * Number of elements which are present in force field.
     ****************************************************************************************/
    Size numberElementsMLFF;
    /***************************************************************************************
     * Number of elements which are present in current structure.
     ****************************************************************************************/
    Size numberElementsStructure;
    /***************************************************************************************
     * Contains the reduced three-body (SHS3) descriptor.
     *
     * First index is the central atom index.
     * Second index is a combined index of
     * l, type1, type2, l, n0, n1.
     ****************************************************************************************/
    Vec2Real&       spectrumVasp;
    ArrayResizing2D spectrumVaspSize;
    /***************************************************************************************
     * Element agnostic Cnlm (intermediate descriptors) summed over types.
     *
     * First index is central atom.
     * Second index is combined index without types.
     * All lnm
     * .
     * .
     * .
     ****************************************************************************************/
    Vec2Real&       clnmVaspNoElement;
    ArrayResizing2D clnmVaspNoElementSize;
    /***************************************************************************************
     * Number of atoms for which the descriptor is stored.
     ****************************************************************************************/
    Int numberAtoms;
    /***************************************************************************************
     * List of centram atom types.
     *
     * @Note The type lists here need a converter
     * from force field type to structure type.
     ****************************************************************************************/
    Vec2Size type0List;
    /***************************************************************************************
     * List of neighbor types 1.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ****************************************************************************************/
    Vec2Size type1List;
    /***************************************************************************************
     * List of neighbor types 2.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ****************************************************************************************/
    Vec2Size type2List;
    /***************************************************************************************
     * List of angular number of spherical harmonic l.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ****************************************************************************************/
    Vec2Size angularList;
    /***************************************************************************************
     * List of first radial numbers (number of roots) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ****************************************************************************************/
    Vec2Size n0List;
    /***************************************************************************************
     * List of second radial number (number of roots) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ****************************************************************************************/
    Vec2Size n1List;
    /***************************************************************************************
     * weight factors to compute n0List n1List terms.
     *
     * The weight factor will be 1 for n0List == n1List
     * and will be sqrt(2) otherwise. This alllows to compute
     * only the upper triangular matrix of n0 and n1 terms.
     ****************************************************************************************/
    Vec2Real weightFactor;
    /***************************************************************************************
     * Controls whether the angular filtering is switched on or not.
     *
     * For more information check the vasp wiki
     * https://www.vasp.at/wiki/index.php/ML_IAFILT2
     ****************************************************************************************/
    bool angularFilterOn;
    /***************************************************************************************
     * Determine which type of angular filtering is used.
     ****************************************************************************************/
    Int angularFilterType;
    /***************************************************************************************
     * Array stores the precomputed angular filter.
     *
     * Dimension is maxOrder + 1, filter value for every l is stored.
     ****************************************************************************************/
    Vec1Real angularFilter;
    /***************************************************************************************
     * Scaling factor to determine the width
     * of the distribution with which the angular filtering is done.
     ****************************************************************************************/
    Real angularFilterScaling;
    /***************************************************************************************
     * Store SHS2 descriptor which was used to build the reduced SHS3 descriptor.
     ****************************************************************************************/
    std::shared_ptr<const DescriptorSHS2> descriptorSHS2;
    /*******************************************************************************************
     * Sparse map which shows the entries in the derivativeMatrix which
     * are needed for the computation of the array forcePreContract.
     *
     * the derivativeMatrix is given by
     * @f[ L^{iJ_{1}}_{ln_{1}n_{2}} =
     *     \frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}
     *     {\partial \hat{\mathbf{X}}_{i}} @f]
     *@f[ i @f] central atom index
     *@f[J_{1} @f] index denoting first neighbor type
     *@f[ l @f] indexing angular parameter
     *@f[ n_{1} @f] radial index belonging to neighbor type 1
     *@f[ n_{2} @f] radial index belonging to neighbor type 2
     * The order in which the indices are described denotes storage order.
     * The sparse list stores those elements of the derivativeMatrix for which the equation is
     *non-zero.
     *
     * Two lists are needed because the derivative has two terms
     * for the reduced descriptor.
     * We need to compute the product of the matrix @f$L^{iJ_{1}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     * For that the following derivative is required:
     * @f[ \frac{\mathrm{d}p^{iJ}_{nn'l}}{\mathrm{d}c^{iJ''}_{n''lm}} = \sqrt{\frac{8\pi^2}{2l+1}}
     * \sum_{m=-l}^{l} \left[\left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) + c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right) \right].@f]
     * We define
     * @f[ A = \left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) @f]
     * and
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f]
     * sparseMap_derivativeMatrix and other arrays without the postfix "B" implement spars list
     * related stuff for the @f[ A @f] term.
     *******************************************************************************************/
    Vec2Size sparseMap_derivativeMatrix;
    /*******************************************************************************************
     * Same as sparseMap_derivativeMatrix, except that spars list related stuff is stored for
     * the second term
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f].
     *******************************************************************************************/
    Vec2Size sparseMap_derivativeMatrix_B;
    /**************************************************************
     * Sparse map which shows the entries in the SHS2-descriptor which
     * are needed for the computation of the array forcePreContract.
     *
     * Two lists are needed because the derivative has two terms
     * for the reduced descriptor.
     * We need to compute the product of the matrix @f$L^{iJ_{1}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     * For that the following derivative is required:
     * @f[ \frac{\mathrm{d}p^{iJ}_{nn'l}}{\mathrm{d}c^{iJ''}_{n''lm}} = \sqrt{\frac{8\pi^2}{2l+1}}
     * \sum_{m=-l}^{l} \left[\left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) + c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right) \right].@f]
     * We define
     * @f[ A = \left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) @f]
     * and
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f]
     * sparseMap_derivativeMatrix and other arrays without the postfix "B" implement spars list
     * related stuff for the @f[ A @f] term.
     **************************************************************/
    Vec2Size sparseMap_SHS2;
    /*******************************************************************************************
     * Same as sparseMap_SHS2, except that spars list related stuff is stored for
     * the second term
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f].
     *******************************************************************************************/
    std::vector<Vec2Size> sparseMap_SHS2_B;
    /*******************************************************************************************
     * Sparse map which shows the entries in the forcePreContract array
     * which are needed for the computation of the array forcePreContract.
     *
     * Two lists are needed because the derivative has two terms
     * for the reduced descriptor.
     * We need to compute the product of the matrix @f$L^{iJ_{1}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     * For that the following derivative is required:
     * @f[ \frac{\mathrm{d}p^{iJ}_{nn'l}}{\mathrm{d}c^{iJ''}_{n''lm}} = \sqrt{\frac{8\pi^2}{2l+1}}
     * \sum_{m=-l}^{l} \left[\left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) + c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right) \right].@f]
     * We define
     * @f[ A = \left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) @f]
     * and
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f]
     * sparseMap_derivativeMatrix and other arrays without the postfix "B" implement spars list
     * related stuff for the @f[ A @f] term.
     *******************************************************************************************/
    Vec2Size sparseMap_central;
    /*******************************************************************************************
     * Same as sparseMap_central, except that spars list related stuff is stored for
     * the second term
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f].
     *******************************************************************************************/
    Vec2Size sparseMap_central_B;
    /*******************************************************************************************
     * Pre-factors for the forcePreContract array computation
     * which are needed for the computation of the array forcePreContract.
     *
     * Two lists are needed because the derivative has two terms
     * for the reduced descriptor.
     * We need to compute the product of the matrix @f$L^{iJ_{1}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     * For that the following derivative is required:
     * @f[ \frac{\mathrm{d}p^{iJ}_{nn'l}}{\mathrm{d}c^{iJ''}_{n''lm}} = \sqrt{\frac{8\pi^2}{2l+1}}
     * \sum_{m=-l}^{l} \left[\left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) + c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right) \right].@f]
     * We define
     * @f[ A = \left(\frac{\mathrm{d}c^{iJ}_{nlm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right)
     * \left(\sum\limits_{J'} c^{iJ'}_{n'lm} \right) @f]
     * and
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f]
     * sparseMap_derivativeMatrix and other arrays without the postfix "B" implement spars list
     * related stuff for the @f[ A @f] term.
     *******************************************************************************************/
    Vec2Real preFactorPreContract;
    /*******************************************************************************************
     * Same as preFactorPreContract, except that spars list related stuff is stored for
     * the second term
     * @f[ B = c^{iJ}_{nlm}\left(\sum\limits_{J'}
     * \frac{\mathrm{d}c^{iJ'}_{n'lm}}{\mathrm{d}c^{iJ''}_{n''lm}} \right). @f].
     *******************************************************************************************/
    Vec2Real        preFactorPreContractB;
    Vec1Size        centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
    Vec1Size        descriptorSize;
    TypeMap         typeMapLoc;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapDescriptorSHS3ReducedLinElem();

} // namespace vaspml

#endif
