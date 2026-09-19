#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "types.hpp"

namespace vaspml
{

class NearestNeighborNSquare;
class TypeMap;

enum class DescriptorType
{
    none = 0,
    bodyOrder2 = 1,
    bodyOrder3 = 2,
    bodyOrder3LinearElement = 3
};

class Descriptor
{
  public:
    Descriptor(ShRec record = nullptr);
    /*******************************************************************************************
     * @param weight weight factor of descriptor used during normalization
     * @param isNormalzied bool which defines if the descriptor is normalized in
     * DescriptorCollector
     *******************************************************************************************/
    Descriptor(const Real           weight,
               const bool           isNormalized,
               const DescriptorType descriptorType,
               ShRec                record = nullptr);
    /*******************************************************************************************
     *return the weight factor which determines the influence of the descriptor in a list of
     *descriptors
     *******************************************************************************************/
    Real get_weight() const;
    /**
     * get controll parameter of the descriptor will be normalized in the DescriptorCollector for
     * example
     */
    bool get_isNormalized() const;
    /*******************************************************************************************
     * @brief Retrieves the type of the descriptor.
     *
     * @return The type of the descriptor as a DescriptorType.
     *******************************************************************************************/
    DescriptorType get_descriptorType() const;
    /**
     *set the neighbor list from which the descriptor is computed
     */
    void set_neighborList(const std::shared_ptr<const NearestNeighborNSquare>& neighborList);
    /**
     * return neighbor list on wich the pair coefficients are based
     *
     * @return const reference to neighbor list
     */
    const NearestNeighborNSquare& get_neighborList() const;
    /**
     * return neighbor list on wich the pair coefficients are based
     *
     * @return syd::shared_ptr to neighbor list
     */
    const std::shared_ptr<const NearestNeighborNSquare>& get_neighborList_ptr() const;
    /** Get type index of central atom of some clnm
     *
     * @param atomIndex array entry of desired central atom
     */
    const Int& get_typeIndexCentral(const Size atomIndex) const;
    /*******************************************************************************************
     * @brief Retrieves the type index of the central atom for a given atom index.
     *
     * This function fetches the type index of the central atom associated with the
     * specified atom index from the neighbor list.
     *
     * @param atomIndex The index of the atom for which the central type index is requested.
     * @return A constant reference to the type index of the central atom.
     *******************************************************************************************/
    const Vec1Int& get_typeIndexCentral() const;
    /** Get type index of indxNeigh neighbor from central atom atomIndex
     *
     * @param atomIndex central atom
     * @param indxNeigh number of the neighbor in NearestNeighborNSquare
     */
    const Int& get_typeIndex(const Size atomIndex, const Size indxNeigh) const;
    /**
     * give number of atoms for which the expansion coefficients were computed
     *
     * @return number of atoms
     */
    Int get_nAtoms() const;
    /**
     *get number of atoms per type from nearest neighbor list
     */
    const Vec1Int& get_nAtomsType() const;
    /**
     * get descriptor array for certain central atom
     *@param centralAtom central atom index for which descriptor is retrieved
     */
    const Vec1Real& get_descriptor(const Size centralAtom) const;
    /**
     * get the size of the descriptor, get number of descriptor for central atom
     *@param centralAtom atom index for which the number of descriptors are returned
     *@NOTE
     *the number of descriptors is independent of the atom index. So all atoms
     *from for example the SHS2 descriptor for a given neighbor list will have
     *the same number of descriptors
     */
    Int get_sizeDescriptor(const Size centralAtom) const;
    Int get_totalNumberBasisFunctions() const;
    /**
     * rescale the descriptors for a certain central atom by a constant factor
     *
     *@param centralAtom atom index for which descriptors are rescaled
     *@param scaleFactor factor by which the elements of the descriptor are scaled
     */
    void rescale_descriptor(const Size centralAtom, const Real scaleFactor);
    /*******************************************************************************************
     * Computing the product of the kernel derivative times the derivative of the descriptor
     * vector with respect to the SHS2-2-body descriptors.
     *
     * @param derivativeMatrix stores the derivative of the kernel matrix with respect
     * to descriptor vector
     * @param forcePreContract stores the result given by the product of the kernel
     * matrix times the descriptor
     * @param typeMap used to convert the structure types to the machine learning force field types
     *
     * This equation shows the what product is computed by the routine.
     * @f[\sum_{B}w_{B}\frac{\partial K( \hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B} ) }
     * {\partial \hat{\mathbf{X}}_{i} } \frac{d\ \hat{\mathbf{X}}_{i}}{c^{iJ_{1}}_{lnm}}@f]
     * For a more information on computing this product check the child classes
     * DescriptorSHS2, DescriptorSHS3
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * compute size for the array where forcePreContract is saved to.
     *******************************************************************************************/
    Size get_forcePreContractSize(const Size atomIndex) const;
    /*******************************************************************************************
     * computing central and pair force terms. pair force terms are derivatives of kernel matrix
     * with respect to neighbor atom positions and central forces are those where derivative is
     *taken w.r.t to central atom
     *
     * @param forcePreContract stores the derivative of the kernel matrix times the 2-body
     *descriptors. This array has to be multiplied by 2-body descriptor derivatives with respect to
     *positions to get force terms.
     * @param pairForces on output stores the derivative of the kernel matrix with respect to the
     *positions of the neighbor atoms for a certain central atom. which is the force pair between
     *two atoms
     * @param centralForces on output stores the derivative of the kernel matrix wrt to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * @brief Resizes internal arrays based on the provided sizes, delegating the operation
     *        to the appropriate derived class based on the descriptor type.
     *
     * @param sizes A vector containing the new sizes for the arrays.
     *
     * This function checks the type of the descriptor (`descriptorType`) and delegates
     * the resizing operation to the corresponding derived class (`DescriptorSHS2` or
     *`DescriptorSHS3`). If the descriptor type is not recognized, no action is taken.
     *******************************************************************************************/
    void resizeArraysRefConf(const Vec1Size& sizes);
    /*******************************************************************************************
     * @brief Fills the descriptor reference configuration for SHS3.
     *
     * This function populates the descriptor reference configuration for a given atom
     * and spectrum data, specifically for descriptors of type `DescriptorSHS3`. If the
     * descriptor type is not `DescriptorType::bodyOrder3`, an exception is thrown.
     *
     * @param atom The index of the atom for which the descriptor is being filled.
     * @param spectrumVasp The spectrum data associated with the atom.
     *
     * @throws std::runtime_error If the descriptor type is not `DescriptorType::bodyOrder3`.
     *******************************************************************************************/
    void fillDescriptorRefConfSHS3(const Size& atom, const Vec1Real& spectrumVasp);
    /*******************************************************************************************
     * @brief Fills the descriptor reference configuration for SHS2.
     *
     * This function populates the descriptor reference configuration for a
     * second-order body descriptor (SHS2). It delegates the operation to the
     * `DescriptorSHS2` implementation if the descriptor type is `bodyOrder2`.
     * Otherwise, it throws an exception indicating that the operation is not
     * supported for other descriptor types.
     *
     * @param atom The index or identifier of the atom for which the descriptor is being filled.
     * @param clnmPair A vector of real values representing pairwise descriptor components.
     * @param clnmPairDerivative A vector of real values representing derivatives of the pairwise
     *descriptor components.
     * @param clnmVasp A vector of real values representing VASP-related descriptor components.
     * @param clnmDerivativeCentralVasp A vector of real values representing derivatives of the
     *central VASP descriptor components.
     *
     * @throws std::runtime_error If the descriptor type is not `DescriptorType::bodyOrder2`.
     *******************************************************************************************/
    void fillDescriptorRefConfSHS2(const Size&     atom,
                                   const Vec1Real& clnmPair,
                                   const Vec1Real& clnmPairDerivative,
                                   const Vec1Real& clnmVasp,
                                   const Vec1Real& clnmDerivativeCentralVasp);
    /*******************************************************************************************
     * @brief Creates a copy of the current Descriptor object with updated data.
     *
     * This method generates a new Descriptor object based on the current instance,
     * but with the provided `data` parameter. If the current descriptor type is
     * `DescriptorType::none`, a new Descriptor is created with the same weight,
     * normalization state, algorithm execution, and the provided data.
     *
     * @param data The ShRec object containing the data to be used in the new Descriptor.
     * @return A new Descriptor object with the updated data.
     *******************************************************************************************/
    Descriptor copy(const ShRec& data = nullptr) const;

    /*******************************************************************************************
     * @brief Precomputes the derivatives required for refitting the descriptor based on the
     *specified atom and its neighbors.
     *
     * This function determines the type of descriptor and delegates the computation to the
     *appropriate specialized implementation. If the descriptor type is not implemented or invalid,
     *an error message is displayed or logged.
     *
     * @param kAtom The index of the atom for which the derivatives are to be precomputed.
     * @param maxNeighborList The list of nearest neighbors for the specified atom.
     *
     * - If `descriptorType` is `DescriptorType::bodyOrder2`, the computation is delegated to
     *`DescriptorSHS2`.
     * - If `descriptorType` is `DescriptorType::bodyOrder3`, the computation is delegated to
     *`DescriptorSHS3`.
     * - If `descriptorType` is `DescriptorType::bodyOrder3LinearElement`, the function outputs a
     *"not implemented" message.
     * - For any other `descriptorType`, an error is logged using `global::tutor.bug`.
     *******************************************************************************************/
    void compute_RefitDescDerivatives(const Int                     kAtom,
                                      const NearestNeighborNSquare& maxNeighborList);

    const Vec2Real& get_refitHeadTerm() const;
    const ShRec&    get_dataRecord() const;
    ShRec&          get_dataRecord();

  protected:
    /*******************************************************************************************
     * Record to store data arrays which might be supplied from the outside world
     *******************************************************************************************/
    ShRec data;
    /// weight of the descriptor when combining several to supervector in DescriptorCollector class
    Real weight;
    /// define if the descriptor will be normalized in the Descriptor collector
    bool isNormalized;
    /// shared ptr to neighbor list from which the descriptor was computed
    std::shared_ptr<const NearestNeighborNSquare> neighborList;
    DescriptorType                                descriptorType;
    /// check if derivative sparse map is already computed. Has to be done once only
    bool derivativeMapComputed;
};

using DescriptorMap = std::map<String, std::shared_ptr<Descriptor>>;

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapDescriptor();

} // namespace vaspml

#endif
