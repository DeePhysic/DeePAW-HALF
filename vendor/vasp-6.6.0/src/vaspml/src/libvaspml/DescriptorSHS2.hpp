#ifndef DESCRIPTORSHS2_HPP
#define DESCRIPTORSHS2_HPP

#include "ArrayResizing.hpp"
#include "Descriptor.hpp"
#include "types.hpp"

namespace vaspml
{

class AtomBatchMap;
class BasisFunctions;
class NearestNeighborNSquare;
class ThermodynIntegration;
class TypeMap;

struct LookUpClnmVasp
{
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index for all index combinations
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * Angular index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size ll;
    /*******************************************************************************************
     * Stores the root Index of the radialBasisFunction and radialBasisFunction derivative for
     * all possible combinations.
     *******************************************************************************************/
    Vec1Size nRoot;
    /*******************************************************************************************
     * Magnetic index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size mm;
    /*******************************************************************************************
     * Count the number of elements in the index vectors
     *
     * @param nAtoms number of central atoms.
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
    /*******************************************************************************************
     * Count number of elements and update index arrays
     *
     * @param nAtoms number of central atoms.
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     * @param lookUpSize defines if look-up tables need a resize or a refill
     *******************************************************************************************/
    void update(const Size&      nAtoms,
                const Size&      lmax,
                const Vec1Size&  nRoots,
                ArrayResizing1D& lookUpSize);
    /*******************************************************************************************
     * Filling all the idex arrays defined in this structure
     *
     * @param nAtoms number of central atoms.
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    void refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
};

// forward declare DescriptorSHS2 to precompute indices of forcePreContract
class DescriptorSHS2;

struct ForcePreContractLookUp
{
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index which can be assigned by the mainIndex
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * Atom type index in the structure which can be accesed with mainIndex
     *******************************************************************************************/
    Vec1Size centralType;
    /*******************************************************************************************
     * Central atom index within its type
     *******************************************************************************************/
    Vec1Size atomIndexPerType;
    /*******************************************************************************************
     * Structure type index of the neighbor atoms beloging to current central atom. Access with
     * mainIndex
     *******************************************************************************************/
    Vec1Size neighborTypeStruc;
    /*******************************************************************************************
     * Entry of derivative matrix
     *******************************************************************************************/
    Vec1Size matrixIndex;
    /*******************************************************************************************
     * Entry index in forcePreContract of DescriptorSHS2
     *******************************************************************************************/
    Vec1Size forcePreIndx;
    /*******************************************************************************************
     * previous number of central atoms to check if resize is needed
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Count the number of linearised arrays
     *
     * @param nAtoms number of central atoms
     * @param nTypeStruc number of types in current structure
     * @param total number of basis functions in DescriptorSHS2
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Size& nTypeStruc, const Size& basisSetSize);
    /*******************************************************************************************
     * Update the look up tables in case the number of central atoms changed
     *
     * @param nAtoms number of atoms
     * @param nTypesStruc number of types in structure
     * @param basisSetSize total number of basis functions in DescriptorSHS2
     * @param typeMap mapping between force field type indices and structure type indices
     * @param typeIndexCentral central atom type index
     * @param atomIndexInType central atom index counted within the corresponding type
     * @param descriptorSHS2 DescriptorSHS2 class for which look up is set up
     *******************************************************************************************/
    void update(const Size&           nAtoms,
                const Size&           nTypesStruc,
                const Size&           basisSetSize,
                const TypeMap&        typeMap,
                const Vec1Int&        typeIndexCentral,
                const Vec1Int&        atomIndexInType,
                const DescriptorSHS2& descriptorSHS2);
    /*******************************************************************************************
     * Refill the index arrays in case the number of atoms changed
     *
     * @param nAtoms number of atoms
     * @param nTypesStruc number of types in structure
     * @param basisSetSize total number of basis functions in DescriptorSHS2
     * @param typeMap mapping between force field type indices and structure type indices
     * @param typeIndexCentral central atom type index
     * @param atomIndexInType central atom index counted within the corresponding type
     *******************************************************************************************/
    void refill(const Size&           nAtoms,
                const Size&           nTypesStruc,
                const Size&           basisSetSize,
                const TypeMap&        typeMap,
                const Vec1Int&        typeIndexCentral,
                const Vec1Int&        atomIndexInType,
                const DescriptorSHS2& descriptorSHS2);
};

struct LookUpRadialSplineNN
{
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index for all index combinations
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * Angular index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size ll;
    /*******************************************************************************************
     * Stores the root Index of the radialBasisFunction and radialBasisFunction derivative for
     * all possible combinations.
     *******************************************************************************************/
    Vec1Size nRoot;
    /*******************************************************************************************
     * Combination index of neighbor atom and ll and nRoot to find in which entry pair description
     * of certain atom pair is stored.
     *******************************************************************************************/
    Vec1Size combiIndex;
    /*******************************************************************************************
     * storing the index for the radial basis functions and its derivatives
     *******************************************************************************************/
    Vec1Size radEntry;
    /*******************************************************************************************
     * Multipicative offset for neighbor atoms.
     *******************************************************************************************/
    Size offset;
    /*******************************************************************************************
     * Number of atoms in old/ previous configuration.
     *
     * This variable is used to check if a resize/refill of the look up table has to be done.
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Count the number of elements in the index vectors
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
    /*******************************************************************************************
     * Count number of elements and update index arrays
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     * @param lookUpSize defines if look-up tables need a resize or a refill
     *******************************************************************************************/
    void update(const Size&      nAtoms,
                const Size&      lmax,
                const Vec1Size&  nRoots,
                ArrayResizing1D& lookUpSize);
    /*******************************************************************************************
     * Filling all the idex arrays defined in this structure
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    void refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
};

struct LookUpYlmNN
{
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index for all index combinations
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * Angular index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size ll;
    /*******************************************************************************************
     * Magnetic index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size mm;
    /*******************************************************************************************
     * Entry index for certain atom pair and set of l,m to which result has to be written in ylm
     *******************************************************************************************/
    Vec1Size ylmEntry;
    /*******************************************************************************************
     * Entry index for direction xyz certain atom pair and set of l,m to which result has to
     * be written in ylmd
     *******************************************************************************************/
    Vec1Size ylmdEntry;
    /*******************************************************************************************
     * Multipicative offset for neighbor atoms.
     *******************************************************************************************/
    Size offset;
    /*******************************************************************************************
     * Number of atoms in old/ previous configuration.
     *
     * This variable is used to check if a resize/refill of the look up table has to be done.
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Count the number of elements in the index vectors
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Size& lmax);
    /*******************************************************************************************
     * Count the number of elements in the index vectors and refill those
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param lookUpSize defines if look-up tables need a resize or a refill
     *******************************************************************************************/
    void update(const Size& nAtoms, const Size& lmax, ArrayResizing1D& lookUpSize);
    /*******************************************************************************************
     * Filling all the idex arrays defined in this structure
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     *******************************************************************************************/
    void refill(const Size& nAtoms, const Size& lmax);
};

struct LookUpclnmPairNN
{
    /*******************************************************************************************
     * index array over which parallelization is done. Index stored in this array has to be used
     * to access remaining arrays.
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * central atom index for all index combinations
     *******************************************************************************************/
    Vec1Size centralAtom;
    /*******************************************************************************************
     * Angular index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size ll;
    /*******************************************************************************************
     * Stores the root Index of the radialBasisFunction and radialBasisFunction derivative for
     * all possible combinations.
     *******************************************************************************************/
    Vec1Size nRoot;
    /*******************************************************************************************
     * Magnetic index of spherical harmonics or Bessel functions for all index combinations
     *******************************************************************************************/
    Vec1Size mm;
    /*******************************************************************************************
     * Stores index to which result for certain atom pair and given lnm has to be stored in clnmPair
     *******************************************************************************************/
    Vec1Size clnmPairEntry;
    /*******************************************************************************************
     * Stores index to which result for certain atom pair and xyz and given lnm has to be stored in
     * clnmPairDerivative
     *******************************************************************************************/
    Vec1Size clnmPairDOffset;
    /*******************************************************************************************
     * Entry index of radial basis functions which are used to compute clnmPair and
     * clnmPairDerivative
     *******************************************************************************************/
    Vec1Size radialEntry;
    /*******************************************************************************************
     * Entry index of spherical harmonics which are used to compute clnmPair and
     * clnmPairDerivative
     *******************************************************************************************/
    Vec1Size ylmEntry;
    /*******************************************************************************************
     * Entry index of spherical harmonics derivative which are used to compute clnmPair and
     * clnmPairDerivative
     *******************************************************************************************/
    Vec1Size ylmdEntry;
    /*******************************************************************************************
     * Number of atoms in old/ previous configuration.
     *
     * This variable is used to check if a resize/refill of the look up table has to be done.
     *******************************************************************************************/
    Size nAtomsOld = 0;
    /*******************************************************************************************
     * Multipicative offset for neighbor atoms in the spline arrays.
     *******************************************************************************************/
    Size offsetRadial;
    /*******************************************************************************************
     * Multipicative offset for neighbor atoms in the ylm arrays.
     *******************************************************************************************/
    Size offsetYLM;
    /*******************************************************************************************
     * Multipicative offset for neighbor atoms in the clnm arrays.
     *******************************************************************************************/
    Size offsetCLNM;
    /*******************************************************************************************
     * Count the number of elements in the index vectors
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    Size countElements(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
    /*******************************************************************************************
     * Count number of elements and update index arrays
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     * @param lookUpSize defines if look-up tables need a resize or a refill
     *******************************************************************************************/
    void update(const Size&      nAtoms,
                const Size&      lmax,
                const Vec1Size&  nRoots,
                ArrayResizing1D& lookUpSize);
    /*******************************************************************************************
     * Filling all the idex arrays defined in this structure
     *
     * @param nAtoms number of central atoms.
     * @param nNeighbors number of neighbors for every central atom
     * @param lmax maximal angular quantum number + 1
     * @param nRoots stores the number of zeros or roots in the Bessel function for a certain l
     *******************************************************************************************/
    void refill(const Size& nAtoms, const Size& lmax, const Vec1Size& nRoots);
};

struct LookUpForceTerms
{
    /*******************************************************************************************
     * Main index which is a combination index of central atom and xyz
     *******************************************************************************************/
    Vec1Size mainIndex;
    /*******************************************************************************************
     * Central atom index which can be retrieved by the mainIndex
     *******************************************************************************************/
    Vec1Size centralAtomIndex;
    /*******************************************************************************************
     * Direction x y z index wich can be retrieved by mainIndex
     *******************************************************************************************/
    Vec1Size xyzIndex;
    /*******************************************************************************************
     * Stores number of atoms to check if the look up tables need resiszing
     *******************************************************************************************/
    Size nAtomsOld = 0;

    void update(const Size& nAtoms);
    void refill(const Size& nAtoms);
};

enum class ClnmPairOrder
{
    neighborMajor,
    paramMajor,
};

class DescriptorSHS2 : public Descriptor
{
  public:
    /*******************************************************************************************
     * @class DescriptorSHS3
     * @brief Two-body descriptor and it's derivative relevant stuff.
     *
     * @param weight descriptor weight when computing normalization constant
     * @param isNormalized determines if descriptor will partake in normalization
     * @param clnmPair supply if memory for clnmPair is supplied from outside
     * @param clnmPairDerivative supply if memory for clnmPair is supplied from outside
     * @param clnmVasp supply if memory for clnmPair is supplied from outside
     * @param clnm_clnmDerivativeCentralVasp supply if memory for clnmPair is supplied from outside
     *
     * This class is used to create a two-body (SHS2) descriptor.
     * This class can for example be created as part of a descriptor collector
     * DescriptorMap which is of std::map type.
     * @code
     * DescriptorMap descriptors;
     * @endcode
     * It can then be accessed with the key "SHS2-2-body" for the two-body descriptors
     * or by "SHS2-3-body" for the two-body part of the thre-body descriptors.
     * The filling for the "SHS2-2-body" descriptors would look like the following (in the Frame
     *class):
     * @code
     * descriptors["SHS2-2-body"] = std::make_shared<DescriptorSHS2>(
     * inputParameters["SHS2-2-body-weight"].cget<Real>(),
     * inputParameters["SHS2-2-body-is-normalized"].cget<bool>(),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair")),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair-derivative")),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS-vasp")),
     * std::get<ShVec2Real>(arrayData.at("2-body-derivative-SHS2-central-vasp")));
     * @endcode
     *
     *******************************************************************************************/
    DescriptorSHS2(const Real           weight = (Real)1.0,
                   const bool           isNormalized = false,
                   const DescriptorType descriptorType = DescriptorType::none,
                   ShRec                record = nullptr);
    /****************************************************************************************
     * Updating the spherical harmonic pair coefficients
     * with the current neighbor list and a supplied basis set.
     *
     * function is only wrapper to decide what kind of algorithm execution will be choosen.
     * For implementation check functions updatePairCoefficientsSingleCPU,
     * updatePairCoefficients_stdparNN
     *
     * @param nn_list nearest neighbor list
     ****************************************************************************************/
    void updatePairCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list);
    /*******************************************************************************************
     * Computes clnmPair for a single atom pair ij and for a single set of parameters.
     *
     * This routine uses lookup tables and is coded in a way to obtain most
     * parallel performance for the c++ standard algorithms.
     *
     * @param atom index of central atom
     * @param radIndx radial basis function index given in the radialBasisfunction array
     * @param ylmIndx spherical harmonics index in ylm which is a combination of neighbor, l and m
     * @param ylmdIndx spherical harmonics derivative index in ylmd which is a combination
     * of xyz,neighbor, l and m
     * @param normed_x normalized x component of vector connecting atom pair i and j
     * @param normed_y normalized y component of vector connecting atom pair i and j
     * @param normed_z normalized z component of vector connecting atom pair i and j
     * @param clnmPair component of @f[ c^{ij}_{nlm} @f] to which result is written
     * @param clnmPairDerivativeX x component @f[ \frac{dc^{ij}_{nlm}}{dx_{j}} @f] to store result
     * @param clnmPairDerivativeY y component @f[ \frac{dc^{ij}_{nlm}}{dy_{j}} @f] to store result
     * @param clnmPairDerivativeZ z component @f[ \frac{dc^{ij}_{nlm}}{dz_{j}} @f] to store result
     *******************************************************************************************/
    void updateSinglePairCoefficients_stdpar(const Size& atom,
                                             const Size& radIndx,
                                             const Size& ylmIndx,
                                             const Size& ylmdIndx,
                                             const Real& normed_x,
                                             const Real& normed_y,
                                             const Real& normed_z,
                                             Real&       clnmPair,
                                             Real&       clnmPairDerivativeX,
                                             Real&       clnmPairDerivativeY,
                                             Real&       clnmPairDerivativeZ);
    /*******************************************************************************************
     * Performs a update step for the DescriptorSHS2 pair coefficients.
     *
     * This routine computes clnmPair and clnmPairderivative. The function exploits parallel
     * standard algorithms with the use of look up-tables.
     *******************************************************************************************/
    void updatePairCoefficients_stdparNN();
    /*******************************************************************************************
     * Performs a update step for the DescriptorSHS2 pair coefficients.
     *
     * This routine computes clnmPair and clnmPairderivative. The function exploits parallel
     * standard algorithms with the use of look up-tables. Routine stores arrays in neighbor
     * major order.
     *******************************************************************************************/
    void updatePairCoefficients_stdparNO();
    /*******************************************************************************************
     * Prepare arrays and look-up tables which will be used when executing the parallel
     * standard algorithms.
     *******************************************************************************************/
    void arrayResizing_stdpar();
    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * This routine selects which algorithm to use (CPU or openACC GPU).
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficients();
    /*******************************************************************************************
     * @brief Computes VASP refit coefficients from precomputed VASP coefficients.
     *
     * This function calculates the VASP refit coefficients based on the VASP coefficients
     * that must have been computed beforehand. It ensures that the required data is available
     * and allocates necessary arrays for the computation. The computation is executed based
     * on the specified execution policy (CPU single-core or GPU). Currently, GPU execution
     * is not implemented.
     *
     * @note The function assumes that `computeVaspCoefficientsFromPairCoefficients()` has
     *       been called prior to this function to compute the required VASP coefficients.
     *       If not, an error is raised.
     *
     * @throws std::runtime_error If the required VASP coefficients are not computed or if
     *         GPU execution is attempted (not implemented).
     *******************************************************************************************/
    void computeVaspRefitCoefficientsFromVaspCoefficients();
    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * CPU implementation using standard library algorithms.
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficientsCPU();
    /*******************************************************************************************
     * @brief Computes the refit coefficients for VASP (Vienna Ab initio Simulation Package)
     *        based on the existing VASP coefficients using CPU processing.
     *
     * This function iterates over all atoms and their respective types to calculate the
     * refit coefficients for VASP. The computation involves subtracting the pairwise
     * derivatives from the central derivatives for each basis set dimension. The results
     * are stored in the `clnmDerivativeRefitVasp` data structure.
     *
     * @note The function assumes that the `neighborList`, `clnmPair`, `clnmPairDerivative`,
     *       `clnmDerivativeCentralVasp`, and `clnmDerivativeRefitVasp` data structures
     *       are properly initialized and populated before calling this function.
     *
     * @warning The code contains a potential infinite loop due to the double semicolon
     *          (`;;`) in the second `for` loop. This should be corrected to avoid runtime issues.
     *
     * @pre `neighborList` must provide valid atom and type counts.
     * @pre All input data structures must have consistent dimensions and sizes.
     *
     * @post The `clnmDerivativeRefitVasp` data structure is updated with the computed refit
     *coefficients.
     *******************************************************************************************/
    void computeVaspRefitCoefficientsFromVaspCoefficientsCPU();
    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * GPU implementation using standard library algorithms.
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficientsNN_stdpar();
    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * GPU implementation using standard library algorithms. The routine uses the
     * neighbor major layout for the pair arrays. The clnmVasp and clnmVaspDerivative
     * will have the usual format.
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficientsNO_stdpar();
    /****************************************************************************************
     * Compute thermodynamic integration coefficients.
     * @param thermodynVars contains storage arrays for thermodynamic integration
     *******************************************************************************************/
    void computeThermodynamicIntegration(ThermodynIntegration& thermodynVars);
    /****************************************************************************************
     * Get the clnmPair vector for a certain centralAtom and all its neighbors.
     * @param centralAtom central atom index for which cnlm pairs are obtained
     *
     * @note The order of the clnm pair is neighbor,l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmPair(const Size centralAtom) const;
    const Vec2Real& get_clnmPair() const;
    /****************************************************************************************
     * Obtain the clnmPairDerivative vector for a certain centralAtom and all its neighbors.
     * @param centralAtom central atom index for which cnlm pairs are obtained
     *
     * @note The order of the clnm pair is neighbor,xyz,l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmPairDerivative(const Size centralAtom) const;
    const Vec2Real& get_clnmPairDerivative() const;
    /****************************************************************************************
     * Get the intermediate expansion coefficients clnmVasp for a certain centralAtom.
     *
     * This is a version of cnlm_pair which was summed over neighbors with same types.
     *
     * @param centralAtom atom index for which the coefficients are returned
     *
     * @note format is atom type, l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmVasp(const Size centralAtom) const;
    const Vec2Real& get_clnmVasp() const;
    /****************************************************************************************
     * Same as get_clnmVasp, but implements the virtual function of
     * Descriptor used in the DescriptorCollector class.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    const Vec1Real& get_descriptor(const Size centralAtom) const;
    /*******************************************************************************************
     * Return the whole descriptor array as const reference. Returns clnmVasp.
     *******************************************************************************************/
    const Vec2Real& get_descriptor() const;
    /****************************************************************************************
     * Get size of spherical harmonics expansion coefficient for a certain central atom.
     *
     * @param centralAtom central atom index for which size is returned
     *******************************************************************************************/
    Int get_sizeDescriptor(const Size centralAtom) const;
    /****************************************************************************************
     * Rescale spherical harmonic expansion coeffcients for a certain central atom.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     * @param scaleFactor scaling factor by which the coefficients are rescaled
     *******************************************************************************************/
    void rescale_descriptor(const Size centralAtom, const Real scaleFactor);
    /****************************************************************************************
     * Get the clnmDerivativeCentralVasp vector for a certain centralAtom.
     *******************************************************************************************/
    const Vec1Real& get_clnmDerivativeCentralVasp(const Size centralAtom) const;
    const Vec2Real& get_clnmDerivativeCentralVasp() const;
    /*******************************************************************************************
     * @brief Retrieves the clnmDerivativeRefitVasp vector for a specific central atom.
     *
     * @param centralAtom The index of the central atom.
     * @return A constant reference to the Vec1Real corresponding to the specified central atom.
     *******************************************************************************************/
    const Vec1Real& get_clnmDerivativeRefitVasp(const Size centralAtom) const;
    const Vec2Real& get_clnmDerivativeRefitVasp() const;
    /****************************************************************************************
     * Return clnmVasp for a certain central atom and a nearest neighbor.
     *
     * @param central_atom central index as nearest neighbor list index
     * @param neighbor     nearest neighbor index as in the NearestNeighborNSquare
     * @param l            angular quantum number of spherical haromnics
     * @param n            determines the nth root of the modified spherical Bessel function
     * @param m            magnetic quantum number of spherical harmonics
     *******************************************************************************************/
    const Real& get_clnmVasp(const Size central_atom,
                             const Size type_index,
                             const Size l,
                             const Size n,
                             const Size m) const;
    /****************************************************************************************
     *
     * Get maximal order of spherical harmonics expansion coefficients.
     *
     * @return maxOrder-1
     *******************************************************************************************/
    Size get_maxOrder() const;
    /****************************************************************************************
     *
     * Get number of roots/radial basis functions for given l.
     *
     * @param l order of spherical harmonic
     *
     * @return number of roots for order l
     *
     *******************************************************************************************/
    const Size& get_nRootsOrder(const Size l) const;
    /****************************************************************************************
     *
     * Get number of roots of radial part of basis functions.
     *
     * @return const referecce to number roots of basis function
     *
     *******************************************************************************************/
    const Vec1Size& get_nRoots() const;
    /****************************************************************************************
     *
     *Return the index for a specifig @[fc^{i}_{lnm}@f].
     *
     *@atom atom index
     *@l angular index (index of angular basis function (spherical harmonic))
     *@n radial index (index of radial basis function)
     *@m second index to define the spherical harmonic
     *******************************************************************************************/
    Size get_Index(const Size atom, const Size l, const Size n, const Size m) const;
    /*******************************************************************************************
     * @brief Computes the X-component index for a given atom and basis set parameters.
     *
     * @param atom The index of the atom.
     * @param l The angular momentum quantum number.
     * @param n The radial quantum number.
     * @param m The magnetic quantum number.
     * @return The computed index for the X-component.
     *******************************************************************************************/
    Size get_IndexX(const Size atom, const Size l, const Size n, const Size m) const;
    /*******************************************************************************************
     * @brief Computes the Y-component index for a given atom and basis set parameters.
     *
     * @param atom The index of the atom.
     * @param l The angular momentum quantum number.
     * @param n The radial quantum number.
     * @param m The magnetic quantum number.
     * @return The computed index for the Y-component.
     *******************************************************************************************/
    Size get_IndexY(const Size atom, const Size l, const Size n, const Size m) const;
    /*******************************************************************************************
     * @brief Computes the Z-component index for a given atom and basis set parameters.
     *
     * @param atom The index of the atom.
     * @param l The angular momentum quantum number.
     * @param n The radial quantum number.
     * @param m The magnetic quantum number.
     * @return The computed index for the Z-component.
     *******************************************************************************************/
    Size get_IndexZ(const Size atom, const Size l, const Size n, const Size m) const;
    /*******************************************************************************************
     * @brief Computes the index for a specific atom, direction, and basis set parameters.
     *
     * @param atom The index of the atom.
     * @param xyz The spatial direction index (0 for x, 1 for y, 2 for z).
     * @param l The angular momentum quantum number.
     * @param n The radial quantum number.
     * @param m The magnetic quantum number.
     * @return The computed index as a Size.
     *******************************************************************************************/
    Size get_IndexDir(const Size atom,
                      const Size xyz,
                      const Size l,
                      const Size n,
                      const Size m) const;
    //    /*******************************************************************************************
    //     * @brief Computes the index for refitting based on the given parameters.
    //     *
    //     * @param xyz        The spatial coordinate index (x, y, z).
    //     * @param atomOrType The atom or type index.
    //     * @param l          The angular momentum quantum number.
    //     * @param n          The radial quantum number.
    //     * @param m          The magnetic quantum number.
    //     * @return The computed index for refitting.
    //     *
    //     * This function calculates the index for refitting using the provided
    //     * spatial coordinate, atom/type, and quantum numbers. The calculation
    //     * incorporates the basis set size and precomputed offsets for angular
    //     * momentum (lOffset).
    //    *******************************************************************************************/
    //    Size get_IndexRefit(const Size xyz, const Size atomOrType, const Size l, const Size n,
    //    const Size m) const;
    /****************************************************************************************
     *
     *Return the basis set size of the two-body descriptor.
     *
     *******************************************************************************************/
    Size get_basisSetSize() const;
    /*******************************************************************************************
     * @brief Retrieves the basis functions associated with the descriptor.
     *
     * This function returns a shared pointer to the `BasisFunctions` object
     * used by the `DescriptorSHS2` instance.
     *
     * @return std::shared_ptr<BasisFunctions> Shared pointer to the basis functions.
     *******************************************************************************************/
    std::shared_ptr<BasisFunctions> get_basisFunctions() const;
    /****************************************************************************************
     *
     * Get offset within clnm for a given l quantum number.
     *
     * Function gives back variable lOffset[l].
     *
     *@param l angular quantum number
     *******************************************************************************************/
    Size get_lOffset(const Size l) const;
    /*******************************************************************************************
     * @brief Retrieves the number of radial basis functions.
     *
     * @return Size The total number of radial basis functions.
     *******************************************************************************************/
    Size get_nRadialBasis() const;
    /*******************************************************************************************
     * @brief Checks if the pairwise descriptors have been completed.
     *
     * @return True if pairwise computations are completed, false otherwise.
     *******************************************************************************************/
    bool get_pairsComputed() const;
    Int  get_totalNumberBasisFunctions() const;
    /*******************************************************************************************
     * @brief Checks if the VASP format descriptors have been completed.
     *
     * @return True if VASP format computations are completed, false otherwise.
     *******************************************************************************************/
    bool get_vaspFormatComputed() const;
    /*******************************************************************************************
     * @brief Checks if the basis functions have been set.
     *
     * @return True if the basis functions are set, false otherwise.
     *******************************************************************************************/
    bool get_basisFunctionsSet() const;
    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     *this class. Overlaying function, will choose between compute_forcePreContractCPU and
     * compute_forcePreContractGPUOpenACC.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     *derivative will be computed in this function and has to be zero on input
     * @param typeMap used to convert  force field types to structure types
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     * this class. CPU version which uses standard algorithms.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     * derivative will be computed in this function and has to be zero on input
     * @param typeMap used to convert  force field types to structure types
     *******************************************************************************************/
    void compute_forcePreContractCPU(const Vec2Real& derivativeMatrix,
                                     Vec2Real&       forcePreContract,
                                     const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     * this class. CPU version which uses standard algorithms.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     * derivative will be computed in this function and has to be zero on input
     *******************************************************************************************/
    void compute_forcePreContract_stdpar(const Vec2Real& derivativeMatrix,
                                         Vec2Real&       forcePreContract) const;
    /*******************************************************************************************
     * Compute the size of the array where the pre contraction of the derivativeMatrix and the
     * descriptor is stored. Routine is needed in the class Predictor.
     *******************************************************************************************/
    Size get_forcePreContractSize(const Size atomIndex) const;
    /*******************************************************************************************
     * Returns the lOffsets used in index computation routines.
     * Takes into account, that only upper triangular part for nRadial1 and nRadial2 is stored.
     *
     * call with Vec1Size& variableName = get_lOffset();
     *******************************************************************************************/
    const Vec1Size& get_lOffset() const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * Routine decides whether CPU or GPUOpenACC version is called.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     VASPML_PALGO_GPU* the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * CPU implementation using standard library algorithms.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTermsCPU(const Vec2Real& forcePreContract,
                              Vec2Real&       pairForces,
                              Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * GPU implementation using standard library algorithms.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTermsGPU(const Vec2Real& forcePreContract,
                              Vec2Real&       pairForces,
                              Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * GPU implementation using standard library algorithms by the use of look up tables.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTerms_stdpar(const Vec2Real& forcePreContract,
                                  Vec2Real&       pairForces,
                                  Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * setting the basis functions which are used for computation of local atom descriptors.
     *
     * @param basisFunctions sets the class basis functions to thois input
     *******************************************************************************************/
    void set_basisFunctions(const std::shared_ptr<BasisFunctions>& basisFunctions);
    /*******************************************************************************************
     * Write clnmPair to file which is specified by fname
     *
     * @param fname name of outputfile to which data is written
     *******************************************************************************************/
    void write_clnmPair(const String& fname) const;
    /*******************************************************************************************
     * Write clnmVasp to file which is specified by fname
     *
     * @param fname name of outputfile to which data is written
     *******************************************************************************************/
    void write_clnmVasp(const String& fname) const;
    /*******************************************************************************************
     * Write clnmDerivativeCentralVasp to file which is specified by fname
     *
     * @param fname name of outputfile to which data is written
     *******************************************************************************************/
    void write_clnmDerivativeCentralVasp(const String& fname) const;
    /*******************************************************************************************
     * Write clnmPairDerivative to file which is specified by fname
     *
     * @param fname name of outputfile to which data is written
     *******************************************************************************************/
    void write_clnmPairDerivative(const String& fname) const;
    /*******************************************************************************************
     * Prepare an array containing the atom types of the central atoms.
     *******************************************************************************************/
    void allocate_typeIndexArray(const Size& numberTypes);
    /*******************************************************************************************
     * Resize look up tables for compute_forcePreContract and computeForceTermsGPU
     *
     * @param typeMap mapping between structure type and atom type
     *******************************************************************************************/
    void resize_LookUpTablesForceCompute(const TypeMap& typeMap);
    /*******************************************************************************************
     * @brief Fills the descriptor data for a specific atom in the SHS2 model.
     *
     * This function updates the descriptor data for a given atom by assigning the provided
     * values for pairwise contributions, their derivatives, VASP contributions, and their
     * central derivatives.
     *
     * @param atom The index of the atom for which the descriptor data is being updated.
     * @param clnmPair The pairwise contributions to the descriptor.
     * @param clnmPairDerivative The derivatives of the pairwise contributions.
     * @param clnmVasp The VASP contributions to the descriptor.
     * @param clnmDerivativeCentralVasp The central derivatives of the VASP contributions.
     *******************************************************************************************/
    void fillDescriptorRefConfSHS2(const Size&     atom,
                                   const Vec1Real& clnmPair,
                                   const Vec1Real& clnmPairDerivative,
                                   const Vec1Real& clnmVasp,
                                   const Vec1Real& clnmDerivativeCentralVasp);
    /*******************************************************************************************
     * @brief Copies parameters from one DescriptorSHS2 object to another.
     *
     * This function takes two DescriptorSHS2 objects as input: `descA` and `descB`.
     * It copies all relevant parameters from `descA` to `descB` by calling the
     * `setParameters` method of `descB` with the values retrieved from `descA`.
     *
     * @param[in] descA The source DescriptorSHS2 object from which parameters are copied.
     * @param[out] descB The target DescriptorSHS2 object to which parameters are copied.
     *
     * The parameters copied include:
     * - Weight (`get_weight`)
     * - Normalization status (`get_isNormalized`)
     * - Descriptor type (`get_descriptorType`)
     * - Basis set size (`get_basisSetSize`)
     * - Number of roots (`get_nRoots`)
     * - Number of radial basis functions (`get_nRadialBasis`)
     * - Maximum order (`get_maxOrder`)
     * - Offset for l (`get_lOffset`)
     * - Computed pairs status (`get_pairsComputed`)
     * - VASP format computation status (`get_vaspFormatComputed`)
     * - Basis functions set (`get_basisFunctionsSet`)
     * - Basis functions (`get_basisFunctions`)
     *
     * Example usage:
     * @code
     * DescriptorSHS2 descA, descB;
     * // Initialize descA with parameters
     * copyParameters(descA, descB);
     * // Now descB has the same parameters as descA
     *******************************************************************************************/
    void setParameters(const Real                             weight,
                       const bool                             isNormalizede,
                       const DescriptorType                   descriptorType,
                       const Size                             basisSetSize,
                       const Vec1Size&                        nRoots,
                       const Size                             nRadialBasis,
                       const Size                             maxOrder,
                       const Vec1Size&                        lOffset,
                       const bool                             pairsComputed,
                       const bool                             vaspFormatComputed,
                       const bool                             basisFunctionsSet,
                       const std::shared_ptr<BasisFunctions>& basisFunctions);
    /*******************************************************************************************
     * @brief Sets the data for the DescriptorSHS2 object.
     *
     * This method initializes the internal data members of the DescriptorSHS2 object
     * with the provided input values. It also calculates the number of atoms based
     * on the size of the `clnmPair` vector.
     *
     * @param clnmPair A vector of real numbers representing the pairwise data.
     * @param clnmPairDerivative A vector of real numbers representing the derivative of the
     *pairwise data.
     * @param clnmVasp A vector of real numbers representing the VASP data.
     * @param clnmDerivativeCentralVasp A vector of real numbers representing the derivative of the
     *central VASP data.
     *
     * @note The size of the `clnmPair` vector is used to determine the number of atoms (`nAtoms`).
     *******************************************************************************************/
    void setData(const Vec2Real& clnmPair,
                 const Vec2Real& clnmPairDerivative,
                 const Vec2Real& clnmVasp,
                 const Vec2Real& clnmDerivativeCentralVasp);
    /*******************************************************************************************
     * @brief Extends the internal data structures with new data.
     *
     * This function appends the provided data vectors to the corresponding internal data
     *structures. It also updates the number of atoms (`nAtoms`) based on the size of the `clnmPair`
     *vector.
     *
     * @param clnmPair A 2D vector containing pairwise data to be appended.
     * @param clnmPairDerivative A 2D vector containing derivative data for pairs to be appended.
     * @param clnmVasp A 2D vector containing VASP-related data to be appended.
     * @param clnmDerivativeCentralVasp A 2D vector containing derivative data for central VASP to
     *be appended.
     *******************************************************************************************/
    void extendData(const Vec2Real& clnmPair,
                    const Vec2Real& clnmPairDerivative,
                    const Vec2Real& clnmVasp,
                    const Vec2Real& clnmDerivativeCentralVasp);
    /*******************************************************************************************
     * @brief Adds a single element to the internal data structures.
     *
     * This function appends individual 1D vectors to the corresponding internal data structures.
     *
     * @param clnmPair A 1D vector containing pairwise data to be added.
     * @param clnmPairDerivative A 1D vector containing derivative data for pairs to be added.
     * @param clnmVasp A 1D vector containing VASP-related data to be added.
     * @param clnmDerivativeCentralVasp A 1D vector containing derivative data for central VASP to
     *be added.
     *******************************************************************************************/
    void addElement(const Vec1Real& clnmPair,
                    const Vec1Real& clnmPairDerivative,
                    const Vec1Real& clnmVasp,
                    const Vec1Real& clnmDerivativeCentralVasp);
    /*******************************************************************************************
     * @brief Reorders the internal data structures based on a batch map.
     *
     * This function reorders the internal data structures (`clnmPair`, `clnmPairDerivative`,
     *`clnmVasp`, and `clnmDerivativeCentralVasp`) using the mapping provided by the `batchMap`. The
     *batch map specifies the new order of atoms relative to their original order.
     *
     * @param batchMap An object of type `AtomBatchMap` that provides the mapping from the new order
     *                 to the original order of atoms.
     *******************************************************************************************/
    void reorderSHS2DataWithBatchMap(const AtomBatchMap& batchMap);
    /*******************************************************************************************
     * @brief Precomputes the derivatives of the refitted descriptor for a given atom.
     *
     * This function calculates the necessary derivatives for the refitted descriptor
     * by invoking the appropriate internal computations.
     *
     * @param kAtom The index of the atom for which the derivatives are to be precomputed.
     * @param maxNeighborList The nearest neighbor list containing information about the atom's
     *neighbors.
     *******************************************************************************************/
    void compute_RefitDescDerivatives(const Int                     kAtom,
                                      const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Computes the dC00 matrix for a given atom based on its neighbors.
     *
     * This function calculates the dC00 matrix, which is used in the computation of
     * symmetry functions for a given central atom. The computation involves iterating
     * over the neighbors of the central atom and their respective types, and populating
     * the dC00 matrix with values derived from the clnmVasp and clnmPairDerivative arrays.
     *
     * @param kAtom The index of the central atom for which the dC00 matrix is computed.
     * @param maxNeighborList The neighbor list containing information about the neighbors
     *                        of the central atom, including their types and global indices.
     *
     * @note The function assumes that the `allocate_dC00()` method has been implemented
     *       to allocate memory for the dC00 matrix. It also assumes that the `neighborList`
     *       object and other required data structures (e.g., `clnmVasp`, `clnmPairDerivative`)
     *       are properly initialized.
     *
     * @details
     * - The function first allocates memory for the dC00 matrix.
     * - It retrieves the type index of the central atom and its neighbors.
     * - For each neighbor, it computes the contributions to the dC00 matrix based on the
     *   symmetry function basis set and the derivatives of the pairwise interactions.
     * - The computation involves iterating over the number of types, roots, and indices
     *   for the symmetry function basis.
     * - Higher orders of angular momentum (l and m) can be added by extending the loops
     *   as indicated in the comments.
     *
     * @warning Ensure that the `neighborList` and `maxNeighborList` objects are consistent
     *          and contain valid data for the given `kAtom`.
     *******************************************************************************************/
    void compute_dC00(const Int kAtom, const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Computes the dC00Head matrix for a given atom based on its type and nearest neighbors.
     *
     * This function calculates and populates the `dC00Head` matrix for the specified atom (`kAtom`)
     * using the provided nearest neighbor list (`maxNeighborList`). The computation involves
     * contributions from the central atom's type and its neighbors' types, as well as their
     * derivatives in the x, y, and z directions.
     *
     * @param kAtom The index of the central atom for which the `dC00Head` matrix is computed.
     * @param maxNeighborList The nearest neighbor list containing information about the central
     *                        atom and its neighbors.
     *
     * @note
     * - The function assumes that higher orders of `l` are not required. If needed, loops for `l`
     *   and `m` should be added where indicated in the code.
     * - The `allocate_dC00Head()` function is called to ensure memory allocation for `dC00Head`.
     * - The `get_Index`, `get_IndexX`, `get_IndexY`, and `get_IndexZ` functions are used to
     *   retrieve the appropriate indices for the computation.
     * - The `clnmVasp` and `clnmPairDerivative` arrays are used as input data for the computation.
     *******************************************************************************************/
    void            compute_dC00Head(const Int kAtom);
    const Vec2Real& get_refitHeadTerm() const;
    void            reorderElements(const Vec1Int& newOrder);

  private:
    /*******************************************************************************************
     * @brief Allocates and initializes the `dC00` data structure based on the maximum number of
     *neighbors.
     *
     * This function resizes the `dC00` vector to accommodate the maximum number of neighbors in the
     *system. Each slice of `dC00` is further resized to match the required size based on the number
     *of types in the neighbor list and the basis set size. All elements are initialized to zero.
     *
     * @param maxNeighborList A reference to the `NearestNeighborNSquare` object that provides the
     *maximum number of neighbors for each atom.
     *
     * @note This function only supports CPU single-core execution. GPU execution is not implemented
     *and will trigger an error.
     *******************************************************************************************/
    void allocate_dC00(const NearestNeighborNSquare& maxNeighborList);
    /*******************************************************************************************
     * @brief Allocates and initializes the `dC00Head` data structure.
     *
     * This function resizes the `dC00Head` vector to match the number of types in the neighbor
     *list. Each slice of `dC00Head` is further resized to the required size based on the number of
     *types in the neighbor list and the basis set size. All elements are initialized to zero.
     *
     * @note This function only supports CPU single-core execution. GPU execution is not implemented
     *and will trigger an error.
     *******************************************************************************************/
    void allocate_dC00Head();
    /*******************************************************************************************
     * Reorder the clnmPairDerivativeNO array and reorder it into clnmPairDerivative
     *
     * Reorder clnmPairDerivativeNO from [centralAtom][l, n, m, neighbor xyz] order to
     * clnmPairDerivative[ centralAtom ][ neighbor x l n m, neighbor y l n m, nighbor z l n m ]
     * order.
     *******************************************************************************************/
    void reorder_clnmPairDerivative_stdpar();
    /*******************************************************************************************
     * Check if a resize of the look up table for forcePreContract is needed.
     *
     * @param typeMap mapping between structure types and force field types
     *******************************************************************************************/
    void resize_forcePreContractLookUp(const TypeMap& typeMap);
    /*******************************************************************************************
     * Check if resize is needed for look up tables of pairForce and centralForce computation
     *******************************************************************************************/
    void resize_lookUpForceTerms();
    /****************************************************************************************
     * Updating the spherical harmonic pair coefficients
     * with the current neighbor list and a supplied basis set.
     *
     * Algorithm will run on CPU. Function is called from updatePairCoefficients
     *
     * @param basisFunctions basis functions for descriptor
     ****************************************************************************************/
    void updatePairCoefficientsCPU(const std::shared_ptr<BasisFunctions>& basisFunctions);
    /*******************************************************************************************
     * Compute two body descriptor expressed as spline interpolated bessel function
     *
     * Routine computes the two body descriptor for a single atom pair and a single basis
     * radial basis function denoted by some index n which is a combination index
     * of the order l of the Bessel function and the root index of the corresponding function.
     *
     * @param distance distance between atom pair ij
     * @param atom central atom index which is one atom of the pair
     * @param neighbor index of neighbor atom in the neighbor list array
     * @param ll order of the used Bessel function
     * @param nRoot root index of the used Bessel function
     * @param combiIndex a combined index of neighbor ll and root to define position to which
     * result is written.
     *******************************************************************************************/
    void computeSplineSinglePair(const Real& distance,
                                 const Size& atom,
                                 const Size& neighbor,
                                 const Size& ll,
                                 const Size& nRoot,
                                 const Size& combiIndex);
    /*******************************************************************************************
     * Compute Offsets for indexing of the storage arrays.
     *
     * generates lOffset shifts between increments of angular quantum number l and stores
     * them to lOffset
     *******************************************************************************************/
    void computeOffsets();
    /*******************************************************************************************
     * Compute pair coefficients clnmPair and the derivatives clnmPairDerivative
     * for a certain central atom.
     *
     * @param index central atom index under consideration
     * @param basisFunctions class containing basis functions to get expansion coefficients
     * clnmPair and the derivatives clnmPairDerivative
     *******************************************************************************************/
    void computePairCoefficientSingleAtom(const Size                             index,
                                          const std::shared_ptr<BasisFunctions>& basisFunctions,
                                          const NearestNeighborNSquare&          neighborList);
    /*******************************************************************************************
     * Compute VASP coeffcients from pair coefficients for a single atom
     *
     * @param index central atom index
     * @param radialBasisFunction radial basis functions for given atom
     * @param radialBasisFunctionDerivative radial basis function derivative for given atom
     *******************************************************************************************/
    void computeVaspCoefficientSingleAtom(const Size index,
                                          Vec1Real&  radialBasisFunction,
                                          Vec1Real&  radialBasisFunctionDerivative);
    /*******************************************************************************************
     * Compute basisSetSize, the total number of basis functions denoted by l,n,m.
     *
     * @param lmax maximum angular quantum number of spherical harmonic basis function
     *******************************************************************************************/
    Size compute_BasisSetSize(const Size lmax, const Vec1Size& nRoots);
    /****************************************************************************************
     * Allocating the arrays cnlm_pair and cnlm_pair_derivative according
     * to the nearest neighbor list topology.
     *
     * @param[ in ] ldim number of spherical harmonics
     *******************************************************************************************/
    void allocatePairCoefficientArrays(const Size& ldim);
    /*******************************************************************************************
     * Prepare an array containing the atom types of the central atoms for CPU code.
     *******************************************************************************************/
    void resize_typeIndexArrayCPU(const Size& numberTypes);
    /*******************************************************************************************
     * Prepare an array containing the atom types of the central atoms for GPU code.
     *******************************************************************************************/
    void resize_typeIndexArrayGPU(const Size& numberTypes);
    /****************************************************************************************
     * Allocate arrays clnmVasp and clnmDerivativeCentralVasp.
     *******************************************************************************************/
    void allocateArraysVaspFormat();
    /*******************************************************************************************
     * @brief Allocates and initializes arrays for VASP refit format based on the execution policy.
     *
     * This function allocates memory for the `clnmDerivativeRefitVasp` array and initializes its
     * elements to zero when the execution policy is set to `cpuSingleCore`. The size of the array
     * is determined by the number of atoms in the neighbor list and the basis set size.
     * If the execution policy is set to `gpuStdLib`, an error message is triggered indicating
     * that the functionality is not implemented.
     *
     * @note This function currently supports only the `cpuSingleCore` execution policy.
     *
     * @throws std::runtime_error If the execution policy is `gpuStdLib`, an error is triggered.
     *******************************************************************************************/
    void allocateArraysVaspRefitFormat();
    /****************************************************************************************
     * Compute the index pairs to assign a specific descriptor derivative term to an entry in the
     * derivative of the kernel.
     *******************************************************************************************/
    Size computeIndexDerivativeMatrix(Size centralAtom,
                                      Size typeIndexForceField,
                                      Size radialIndex,
                                      Size totalShift);
    /*******************************************************************************************
     * Resizes the clnmPair and clnmPairDerivative array. This routine is used when std::par
     * is used.
     *
     * This function helps to avoid unnecessary host to device and device to host memory
     * copies in case a resize is needed
     *******************************************************************************************/
    void resizePairArraysGPU(const Size& ldim);
    /*******************************************************************************************
     * Resizes the clnmVasp and clnmDerivativeCentralVasp array. This routine is used when std::par
     * is used.
     *
     * This function helps to avoid unnecessary host to device and device to host memory
     * copies in case a resize is needed
     *******************************************************************************************/
    void resizeArraysVaspFormatGPU();
    /*******************************************************************************************
     * @brief Resizes the `centralAtomIndex` array to match the number of atoms in the neighbor
     *list.
     *
     * This function checks if the size of the `centralAtomIndex` array needs to be updated
     * based on the number of atoms in the neighbor list. If resizing is required, the array
     * is resized to the appropriate size and populated with sequential indices starting from 0.
     *
     * @details
     * - The number of atoms is retrieved from the `neighborList` object using its
     *`get_numberAtoms()` method.
     * - The `checkResize()` method of `centralAtomIndexSize` determines whether resizing is
     *necessary.
     * - If resizing is required:
     *   - The `centralAtomIndex` array is resized to match the number of atoms.
     *   - The array is populated with sequential indices using `std::iota`.
     *
     * @pre The `neighborList` object must be properly initialized and contain valid atom data.
     * @post If resizing occurs, the `centralAtomIndex` array will contain indices ranging from 0 to
     *`nAtoms - 1`.
     *
     * @note This function is efficient and avoids unnecessary resizing if the current size of
     * `centralAtomIndex` matches the number of atoms.
     *******************************************************************************************/
    void resize_centralAtomIndex();
    /*******************************************************************************************
     *  Prepare the radial basis functions such that the basis function can be stored for every
     *  atom pair and radial basis function index.
     *******************************************************************************************/
    void resize_radialBasisFunction();
    /*******************************************************************************************
     * Prepare ylm and ylmd in way to store all descriptors for every atom pair and ylm index
     *******************************************************************************************/
    void resize_ylm(const Size& ldim);
    /****************************************************************************************
     * clnmPair summed over nearest neighbors.
     *
     * type resolution is kept
     * first index is central atom
     * second index is combined index
     * type 1 all lnm
     * type 2 all lnm
     * .
     * .
     * .
     *******************************************************************************************/
    Vec2Real& clnmVasp;
    /*******************************************************************************************
     * Structure decides if the clnmVasp array needs a resize. Resize is only done if array grows.
     *
     * If array shrinkgs no resize is done only the new array dimension is stored in the object.
     *******************************************************************************************/
    ArrayResizing2D clnmVaspSize;
    /*******************************************************************************************
     * Look up table to flatten loop to compute the clnmVASP.
     *
     * The look up table will flatten the outer loop over central atom, l,nRoot, and m
     *******************************************************************************************/
    LookUpClnmVasp lookUpClnmVasp;
    /*******************************************************************************************
     * Object to check if size of the look up table for clnmVASP computation has to be resized
     *******************************************************************************************/
    ArrayResizing1D lookUpClnmVaspSize;
    /****************************************************************************************
     * clnmPairDerivative summed over nearest neighbors.
     *
     * type resolution is kept
     * first index is central atom
     * second index is combined index
     * x component type 1 all lnm
     * y component type 1 all lnm
     * z component type 1 all lnm
     * x component type 2 all lnm
     * y component type 2 all lnm
     * z component type 2 all lnm
     * .
     * .
     * .
     *******************************************************************************************/
    Vec2Real& clnmDerivativeCentralVasp;
    Vec2Real& clnmDerivativeRefitVasp;
    /*******************************************************************************************
     * Structure decides if the clnmDerivativeCentralVaspSize array needs a resize.
     *
     * Resize is only done if array grows.
     * If array shrinkgs no resize is done only the new array dimension is stored in the object.
     *******************************************************************************************/
    ArrayResizing2D clnmDerivativeCentralVaspSize;
    /****************************************************************************************
     * Store expansion coefficients of the nearest neighbor environment.
     *
     *  index 1 is the central atom, index 2 is
     *  neighbor 1 all lnm
     *  neighbor 2 all lnm
     *  neighbor 3 all lnm...
     *******************************************************************************************/
    Vec2Real& clnmPair;
    /*******************************************************************************************
     * Structure decides if the clnmPair array needs a resize. Resize is only done if array grows.
     *
     * If array shrinkgs no resize is done only the new array dimension is stored in the object.
     *******************************************************************************************/
    ArrayResizing2D clnmPairSize;
    /****************************************************************************************
     * Derivatives of the expansion coefficients of the nearest neighbor environment.
     *
     * index 1 central atom
     * index 2 neighbor atom x all nlm, y all nlm, z all nlm
     *******************************************************************************************/
    Vec2Real& clnmPairDerivative;
    /****************************************************************************************
     * Derivatives of the expansion coefficients of the nearest neighbor environment.
     *
     * index 1 central atom
     * index 2 l0 n0 m0 neighbor atom 0 x y z, l0 n0 m0 neighbor atom 1 x, y, z ...
     * l1 n1 m1 neighbor atom 0 x, y, z
     *******************************************************************************************/
    Vec2Real& clnmPairDerivativeNO;
    /*******************************************************************************************
     * Structure decides if the clnmPairDerivative array needs a resize. Resize is only done if
     *array grows.
     *
     * If array shrinkgs no resize is done only the new array dimension is stored in the object.
     *******************************************************************************************/
    ArrayResizing2D clnmPairDerivativeSize;
    /****************************************************************************************
     * Total number of basis functions (n,l,m).
     *******************************************************************************************/
    Size basisSetSize;
    /****************************************************************************************
     * Number of roots of radial function.
     *
     * @todo Maybe this should become a shared pointer at some point, also in the descriptorRadial
     * spline routine, like this a copy constructor will be called.
     *******************************************************************************************/
    Vec1Size nRoots;
    /****************************************************************************************
     * Total number of radial basis functions.
     *******************************************************************************************/
    Size nRadialBasis;
    /****************************************************************************************
     * Number of atoms in the nearest neighbor list that was supplied to construct the object.
     *******************************************************************************************/
    Size nAtoms;
    /****************************************************************************************
     * Max Order of spherical haromincs, stores max angular momentum + 1 for loop writing
     *convenience.
     *
     *@note When returned with getter the the plus 1 will be removed.
     *******************************************************************************************/
    Size maxOrder;
    /****************************************************************************************
     * Array offset for indexing.
     *******************************************************************************************/
    Vec1Size lOffset;
    /****************************************************************************************
     * To check if pair coefficients are computed.
     *******************************************************************************************/
    bool pairsComputed;
    /****************************************************************************************
     * To check if calculations in vasp format were carried out.
     *******************************************************************************************/
    bool vaspFormatComputed;
    /*******************************************************************************************
     * Array used to store the interpolated radial basis functions for every atom neighbor pair.
     *
     * First index is central atom index. Second index is a combindation of neighbor, besserl order
     * and root index of the radial basis functions. This array is only used in combination
     * with std parallel execution policy.
     *******************************************************************************************/
    Vec2Real radialBasisFunction;
    /*******************************************************************************************
     * Array used to store derivative of radial basis functions for every atom neighbor pair.
     *
     * First index is central atom index. Second index is a combindation of neighbor, besserl order
     * and root index of the radial basis functions. This array is only used in combination
     * with std parallel execution policy.
     *******************************************************************************************/
    Vec2Real radialBasisFunctionDerivative;
    /*******************************************************************************************
     * Array to check if a resize is necessary for radialBasisFunction or
     *radialBasisFunctionDerivative
     *
     * only used in parallel std version and only important for the GPU version to prevent
     * unnecessary memory copies from host to device and back
     *******************************************************************************************/
    ArrayResizing2D radialBasisFunctionSize;
    /*******************************************************************************************
     * Variable checks if basis functions were already supplied to the object.
     *
     * If no basisFunctions are supplied the object can not compute descriptor since those
     * have to be expressed in some basis set.
     *******************************************************************************************/
    bool basisFunctionsSet;
    /*******************************************************************************************
     * Storing the the angular descriptors for certain atom pair in term of spherical harmonics.
     *
     * The ordering of the array is central atom in first index. The second index is a
     * collapsed index over neighbor atom, l index, and m index as fastest.
     *******************************************************************************************/
    Vec2Real ylm;
    /*******************************************************************************************
     * Controlls size of the ylm array to store angular descriptors.
     *
     * Object is only used in in the case when algorithm is executed in std parallel execution.
     *******************************************************************************************/
    ArrayResizing2D ylmSize;
    /*******************************************************************************************
     * Storing the the angular descriptors for certain atom pair in term of spherical harmonics.
     *
     * The ordering of the array is central atom in first index. The second index is a
     * collapsed index over neighbor atom, l index, and m index as fastest.
     *******************************************************************************************/
    Vec2Real ylmDerivative;
    /*******************************************************************************************
     * Controlls size of the ylmDerivative array to store angular descriptors.
     *
     * Object is only used in in the case when algorithm is executed in std parallel execution.
     *******************************************************************************************/
    ArrayResizing2D ylmDerivativeSize;
    /*******************************************************************************************
     * stores shared pointer to the basis functions object which were used to compute the
     * DescriptorSHS2.
     *******************************************************************************************/
    std::shared_ptr<BasisFunctions> basisFunctions;
    /*******************************************************************************************
     *******************************************************************************************/
    ForcePreContractLookUp forcePreContractLookUp;
    /*******************************************************************************************
     * Look up table for the radial spline computation to flatten loops into single loop
     *
     * Flattens loop over central atoms, and indices of the radial basis functions, which
     * are l,nRoots. For more information check the documentation in the definition of the look up
     *table
     *******************************************************************************************/
    LookUpRadialSplineNN lookUpRadialSplineNN;
    /*******************************************************************************************
     * Object to check if  a resize of the lookUpRadialSplineNN is needed.
     *
     * Resize would be necessary if number of atoms, lmax or max number of roots changes
     *******************************************************************************************/
    ArrayResizing1D lookUpRadialSplineNNSize;
    /*******************************************************************************************
     * Stores look up table for flattening loop to compute spherical expansion coefficients.
     *
     *
     * Loop flattens central atom, l and magnetic index m. Further information can be found
     * in the definition of LookUpYlmNN.
     *******************************************************************************************/
    LookUpYlmNN lookUpYlmNN;
    /*******************************************************************************************
     * Object to check if a resize of the lookUpYlmNN is needed.
     *
     * Resize is needed if number of central atoms changes lmax changes.
     *******************************************************************************************/
    ArrayResizing1D lookUpYlmNNSize;
    /*******************************************************************************************
     * Stores look up table for flattening loop to compute pair descriptors clnm.
     *
     * Loop flattens central atom, l, nRoot and magnetic index m. Further information can be found
     * in the definition of LookUpYlmNN.
     *******************************************************************************************/
    LookUpclnmPairNN lookUpclnmPairNN;
    /*******************************************************************************************
     * Check if resize for lookUpclnmPairNN is needed.
     *
     * Resize is needed if number of central atoms changes lmax changes or nRoots changes.
     *******************************************************************************************/
    ArrayResizing1D lookUpclnmPairNNSize;
    /*******************************************************************************************
     * Look up table to flatten loop over central atom and xyz.
     *******************************************************************************************/
    LookUpForceTerms lookUpForceTerms;
    /*******************************************************************************************
     * Controll variable to determine if the clnmPair and clnmPair derivatives are computed
     * in neighbor major or basis function major order.
     *******************************************************************************************/
    ClnmPairOrder   clnmPairOrder;
    Vec1Size        centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
    // TODO: put those in a look up table
    Vec1Size        typeIndexArray;
    ArrayResizing1D typeIndexArraySize;

    Vec2Real dC00Head;
    Vec2Real dC00;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapDescriptorSHS2();

/*******************************************************************************************
 * @brief Copies all parameters from one DescriptorSHS2 object to another.
 *
 * This function copies the parameters of `descA` into `descB` by invoking
 * the `setParameters` method of `descB` with the values retrieved from `descA`.
 *
 * @param descA The source DescriptorSHS2 object from which parameters are copied.
 * @param descB The target DescriptorSHS2 object to which parameters are copied.
 *******************************************************************************************/
void copyParameters(const DescriptorSHS2& descA, DescriptorSHS2& descB);
/*******************************************************************************************
 * @brief Copies the neighbor list from one DescriptorSHS2 object to another.
 *
 * This function copies the neighbor list pointer from `descA` to `descB` using
 * the `set_neighborList` method of `descB`.
 *
 * @param descA The source DescriptorSHS2 object from which the neighbor list is copied.
 * @param descB The target DescriptorSHS2 object to which the neighbor list is copied.
 *******************************************************************************************/
void copyNeighborList(const DescriptorSHS2& descA, DescriptorSHS2& descB);
/*******************************************************************************************
 * @brief Copies basis functions from one DescriptorSHS2 object to another.
 *
 * This function copies the basis functions from `descA` to `descB` using
 * the `set_basisFunctions` method of `descB`.
 *
 * @param descA The source DescriptorSHS2 object from which basis functions are copied.
 * @param descB The target DescriptorSHS2 object to which basis functions are copied.
 *******************************************************************************************/
void copyBasisFunctions(const DescriptorSHS2& descA, DescriptorSHS2& descB);
/*******************************************************************************************
 * @brief Copies data arrays from one DescriptorSHS2 object to another.
 *
 * This function copies various data arrays from `descA` to `descB` using
 * the `setData` method of `descB`.
 *
 * @param descA The source DescriptorSHS2 object from which data arrays are copied.
 * @param descB The target DescriptorSHS2 object to which data arrays are copied.
 *******************************************************************************************/
void copyData(const DescriptorSHS2& descA, DescriptorSHS2& descB);
/*******************************************************************************************
 * @brief Extends data arrays in the target DescriptorSHS2 object using data from the source.
 *
 * This function extends the data arrays in `descB` by appending the data arrays
 * from `descA` using the `extendData` method of `descB`.
 *
 * @param descA The source DescriptorSHS2 object from which data arrays are appended.
 * @param descB The target DescriptorSHS2 object to which data arrays are extended.
 *******************************************************************************************/
void extendData(const DescriptorSHS2& descA, DescriptorSHS2& descB);
/*******************************************************************************************
 * @brief Adds element-specific data from the source DescriptorSHS2 object to the target.
 *
 * This function adds data for a specific element (identified by `atomIdx`) from `descA`
 * to `descB` using the `addElement` method of `descB`.
 *
 * @param descA The source DescriptorSHS2 object from which element-specific data is added.
 * @param descB The target DescriptorSHS2 object to which element-specific data is added.
 * @param atomIdx The index of the atom for which data is added.
 *******************************************************************************************/
void addElementData(const DescriptorSHS2& descA, DescriptorSHS2& descB, const Size atomIdx);

} // namespace vaspml

#endif
