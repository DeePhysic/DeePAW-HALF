#ifndef SPILLINGFACTOR_HPP
#define SPILLINGFACTOR_HPP

#include "types.hpp"

namespace vaspml
{

class MlMPI;
struct Record;
class TypeMap;
template<class T>
class ShmemArray2DVariableLen;

namespace setup
{
std::shared_ptr<ShmemArray2DVariableLen<Real>> makeInverseCovMatrixArray(
    const Record&                 inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn);
}

MapString dataMapSpillingFactor();

/*******************************************************************************************
 * Class can compute the spilling factor for a set of atoms.
 *
 * To compute a spilling factor over a set of atoms the Kernel matrix has to be supplied
 * to the class constructor. The kernel matrix here is the matrix between the local reference
 * configurations of the basis functions.
 * The spilling factor gives a measure for the quality of the force field on the current
 * structure. A spilling factor states that the prediction is trustworthy. A spilling factor
 * close to 1 states the force field is extrapolating and the prediction can not
 * be trusted.
 *******************************************************************************************/
class SpillingFactor
{
  public:
    /*******************************************************************************************
     * @brief Constructor for the SpillingFactor class.
     *
     * Initializes the SpillingFactor object by setting up data members based on input parameters,
     * MPI context, and shared record data. It assigns or creates the data record, extracts
     * the spilling factor vector, initializes local reference configurations, sets up the
     * inverse kernel matrix, and initializes computation and statistics flags and values.
     *
     * @param inputParameters A constant reference to a Record object containing various input parameters,
     *                        including "numLrc" for the number of local reference configurations
     *                        and parameters for inverse covariance matrix setup.
     * @param mpiIn A shared pointer to an MlMPI object, providing MPI context for distributed operations,
     *              specifically used in setting up the inverse kernel matrix.
     * @param dataIn A ShRec (shared record) object which is used to initialize the 'data' member.
     *               If 'dataIn' is null or empty, a new Record is created using 'dataMapSpillingFactor()'.
     *******************************************************************************************/
    SpillingFactor(const Record&                 inputParameters,
                   const std::shared_ptr<MlMPI>& mpiIn = nullptr,
                   const ShRec&                  dataIn = nullptr);
    /*******************************************************************************************
     * computes the spilling factor for the current configuration of atoms
     *
     * @param kernelMatrix stores the kernel matrix of the current configuration and should
     * be retrived from some kernel class as KernelPolynomial
     * @param nAtomsType determines the number of atoms per type
     * @param typeMap map which gives the relationship between force field and structure atom
     * types
     *******************************************************************************************/
    void computeSpillingFactor(const Vec2Real& kernelMatrix,
                               const Vec1Int&  nAtomsType,
                               const TypeMap&  typeMap);
    /*******************************************************************************************
     * computin statistics of the spilling factor
     *
     * statistics include min and max value of spilling factor per type. The average spilling
     * factor per atom type. Total average, minimum value and maximum value are computed
     * over the various atom types.
     *******************************************************************************************/
    void computeStatistics();
    /*******************************************************************************************
     * get spilling factor for certain atom of specified type.
     *
     * @param type determines atom type in current structure
     * @param atom determines number of atom within type of current structure
     *******************************************************************************************/
    const Real& get_spillingFactor(const Size type, const Size atom) const;
    /*******************************************************************************************
     * get vector of spilling factors for certain atom type.
     *
     * @param type determines atom type in current structure
     *******************************************************************************************/
    const Vec1Real& get_spillingFactor(const Size type) const;
    /*******************************************************************************************
     * get spilling factor array in format [type][atom]
     *******************************************************************************************/
    const Vec2Real& get_spillingFactor() const;
    /*******************************************************************************************
     * get the average spilling factor averaged over atom types and atoms within type
     *******************************************************************************************/
    const Real& get_averageSpillingFactor() const;
    /*******************************************************************************************
     * get the maximum spilling factor maximized over atom types and atoms within type
     *******************************************************************************************/
    const Real& get_maxSpillingFactor() const;
    /*******************************************************************************************
     * get the minimum spilling factor minimized over atom types and atoms within type
     *******************************************************************************************/
    const Real& get_minSpillingFactor() const;
    /*******************************************************************************************
     * get average spilling factor for various atom types. Average is computed over atoms
     * in type
     *******************************************************************************************/
    const Vec1Real& get_averageSpillingFactorType() const;
    /*******************************************************************************************
     * get max spilling factor for various atom types. Average is computed over atoms in type
     *******************************************************************************************/
    const Vec1Real& get_maxSpillingFactorType() const;
    /*******************************************************************************************
     * get min spilling factor for various atom types. Average is computed over atoms in type
     *******************************************************************************************/
    const Vec1Real& get_minSpillingFactorType() const;
    /*******************************************************************************************
     * get average spilling factor for various atom types. Average is computed over atoms
     * in type
     *
     * @param type for which the average spilling factor should be returned
     *******************************************************************************************/
    const Real& get_averageSpillingFactorType(const Size type) const;
    /*******************************************************************************************
     * get max spilling factor for various atom types. Average is computed over atoms in type
     *
     * @param type for which the maximum spilling factor should be returned
     *******************************************************************************************/
    const Real& get_maxSpillingFactorType(const Size type) const;
    /*******************************************************************************************
     * get min spilling factor for various atom types. Average is computed over atoms in type
     *
     * @param type for which the minimum spilling factor should be returned
     *******************************************************************************************/
    const Real& get_minSpillingFactorType(const Size type) const;
    /*******************************************************************************************
     * Write the computed spilling factor to screen.
     *
     * First column denotes atom number, second column is place holder third column is
     * atom type in integer code. And last column contains the spilling factor
     *******************************************************************************************/
    void writeToScreen() const;
  private:
    /*******************************************************************************************
     * Data Record used as data storage for arrays which might need communication between
     * different classes.
     *******************************************************************************************/
    ShRec data;
    /*******************************************************************************************
     * Stores the array of spilling factors for every atom.
     *
     * The first index accesses the atom type and the second index is the atom number
     * within an atom type
     *******************************************************************************************/
    Vec2Real& spillingFactor;
    /*******************************************************************************************
     * Number of local reference configurations per atom type.
     *
     * The total size of the kernel matrix is the numberLocalRefConfs squared for every
     * atom type.
     *******************************************************************************************/
    ShCVec1Int numberLocalRefConfs;
    /*******************************************************************************************
     * stores the inverse kernel matrix of the local reference configurations
     *******************************************************************************************/
    std::shared_ptr<ShmemArray2DVariableLen<Real>> inverseKernelMatrix;
    /*******************************************************************************************
     * Check if the spillingFactor array is computed
     *******************************************************************************************/
    bool isComputed;
    /*******************************************************************************************
     * Check if the statistics from the spillingFactor array are computed
     *******************************************************************************************/
    bool statisticsComputed;
    /*******************************************************************************************
     * Stores the average spillingFactor. Average is over all atoms, no type resolution.
     *******************************************************************************************/
    Real averageSpillingFactor;
    /*******************************************************************************************
     * Stores the maximum spillingFactor. Maximum is over all atoms, no type resolution.
     *******************************************************************************************/
    Real maxSpillingFactor;
    /*******************************************************************************************
     * Stores the minimum spillingFactor. Minimum is over all atoms, no type resolution.
     *******************************************************************************************/
    Real minSpillingFactor;
    /*******************************************************************************************
     * Stores the average spillingFactor with type resolution.
     *
     * Average is computed over the atoms within a single type.
     *******************************************************************************************/
    Vec1Real averageSpillingFactorType;
    /*******************************************************************************************
     * Stores the maximum spillingFactor with type resolution.
     *
     * Maximum is computed over the atoms within a single type.
     *******************************************************************************************/
    Vec1Real maxSpillingFactorType;
    /*******************************************************************************************
     * Stores the minimum spillingFactor with type resolution.
     *
     * Minimum is computed over the atoms within a single type.
     *******************************************************************************************/
    Vec1Real minSpillingFactorType;
};

} //namespace vaspml

#endif
