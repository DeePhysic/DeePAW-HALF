#ifndef FITMATRIX_HPP
#define FITMATRIX_HPP

#include "types.hpp"

namespace vaspml
{

class KernelPolynomial;
class MlMPI;
template<class T>
class ShmemArray2D;

class FitMatrix
{
  public:
    FitMatrix(ShRec record = nullptr);
    FitMatrix(const Vec1Int&         nAtoms,
              const Int              nRefConfs,
              std::shared_ptr<MlMPI> mpiIn = nullptr,
              ShRec                  record = nullptr);
    /*******************************************************************************************
     * @brief Initializes the design matrix mode for the FitMatrix object.
     *
     * This function sets up the shared memory array for the design matrix and computes
     * the force offsets based on the number of atoms per structure. The design matrix
     * is allocated with dimensions determined by the computed size of the design matrix
     * and the number of reference conformations. The force offsets are calculated using
     * a combination of partial sums and transformations to account for the structure and
     * atom-specific offsets.
     *
     * @note This function assumes that `nAtomsPerStructure`, `nRefConfs`, `nStrucs`, and `mpi`
     *       are properly initialized prior to calling this function.
     *******************************************************************************************/
    void initDesignMatrixMode();
    /*******************************************************************************************
     * @brief Initializes the normal equation mode for the FitMatrix class.
     *
     * This function sets up the necessary data structures and memory for performing
     * normal equation-based fitting. It initializes offsets, allocates shared memory
     * for the fit matrix, resizes the temporary design matrix, and prepares the target
     * data for the normal equation.
     *
     * Steps performed:
     * - Sets the first force offset to the number of structures.
     * - Allocates shared memory for the fit matrix with dimensions `nRefConfs x nRefConfs`.
     * - Determines the maximum number of atoms across all structures.
     * - Resizes the temporary design matrix based on the computed design matrix size.
     * - Initializes the target data for the normal equation with `nRefConfs` elements.
     *******************************************************************************************/
    void initNormalEqMode();
    /*******************************************************************************************
     * @brief Adds structural contributions to the design matrix.
     *
     * This function incorporates both energy and force contributions of a given structure,
     * specified by its index, into the design matrix using the provided kernel polynomial.
     *
     * @param kernelPoly The kernel polynomial used for the computation.
     * @param strucIdx The index of the structure to be added.
     *******************************************************************************************/
    void addStructure(const KernelPolynomial& kernelPoly, const Int strucIdx);
    /*******************************************************************************************
     * @brief Adds the design matrix components for a specific structure.
     *
     * This function updates the design matrix by incorporating energy and force contributions
     * for the specified structure index using the provided kernel polynomial.
     *
     * @param kernelPoly The kernel polynomial used for the design matrix computation.
     * @param strucIdx The index of the structure to be added to the design matrix.
     *******************************************************************************************/
    void addStructureDesignMatrix(const KernelPolynomial& kernelPoly, const Int strucIdx);
    /*******************************************************************************************
     * @brief Adds a structure's contribution to the normal equation for fitting a matrix.
     *
     * This function incorporates the energy, force, and stress targets of a structure into the
     * normal equation used for fitting. It computes the design matrix contributions and updates
     * the normal equation matrix and target vector accordingly.
     *
     * @param kernelPoly      The kernel polynomial object representing the structure's kernel data.
     * @param energyTarget    The target energy value for the structure.
     * @param forceTarget     The target force values for the structure, provided as a vector.
     * @param stressTarget    The target stress values for the structure, provided as a vector.
     *
     * The function performs the following steps:
     * - Adds the structure's energy and force contributions to the temporary design matrix.
     * - Computes the size of the design matrix based on the number of atom types and other
     *parameters.
     * - Updates the normal equation matrix by computing \( A^T \cdot A \), where \( A \) is the
     *design matrix.
     * - Constructs a contiguous target vector from the energy, force, and stress targets.
     * - Computes \( A^T \cdot y \), where \( y \) is the target vector, and adds it to the normal
     *equation target.
     *******************************************************************************************/
    void addStructureNormalEquation(const KernelPolynomial& kernelPoly,
                                    const Real&             energyTarget,
                                    const Vec1Real&         forceTarget,
                                    const Vec1Real&         stressTarget);
    /*******************************************************************************************
     * @brief Adds the structure energy contribution to the design matrix.
     *
     * This function computes and adds the contribution of a structure's energy
     * to the design matrix using the provided kernel polynomial. The energy
     * contribution is normalized by the total number of atoms in the structure.
     *
     * @param kernelPoly The kernel polynomial object containing the kernel matrix,
     *                   number of atoms per type, and number of reference configurations.
     * @param strucIdx   The index of the structure for which the energy contribution
     *                   is being added to the design matrix.
     *
     * The function iterates over all atom types, reference configurations, and atoms
     * to compute the normalized kernel matrix values and updates the design matrix
     * accordingly.
     *******************************************************************************************/
    void addStructureEnergy(const KernelPolynomial& kernelPoly,
                            const Int               strucIdx,
                            Real*                   fitMatrix);
    /*******************************************************************************************
     * @brief Adds structure force contributions to the fit matrix for a given structure index.
     *
     * This function computes the force contributions for a given structure index and updates the
     * fit matrix accordingly. The computation involves the use of kernel polynomial derivatives,
     * normalized descriptor vectors, and reference configurations. Both direct and indirect terms
     * are calculated for each atom and its neighbors.
     *
     * @param kernelPoly   The kernel polynomial object containing descriptors, derivatives, and
     *related data.
     * @param fitMatrix    Pointer to the fit matrix where the computed force contributions will be
     *added.
     * @param forceOffset  Offset in the fit matrix for the current structure's force contributions.
     *
     * The function performs the following steps:
     * - Extracts kernel polynomial derivatives and descriptor data.
     * - Iterates over all atoms and their neighbors to compute direct and indirect force terms.
     * - Computes outer products, dot products, and applies kernel polynomial derivatives.
     * - Updates the fit matrix with the computed force contributions.
     *
     * @note This function assumes that the input data (e.g., descriptors, neighbor lists) is
     *precomputed and consistent with the kernel polynomial object.
     *******************************************************************************************/
    void addStructureForce(const KernelPolynomial& kernelPoly,
                           Real*                   fitMatrix,
                           const Int               forceOffset = 0);
    /*******************************************************************************************
     * @brief Writes the design matrix to a file.
     *
     * This function saves the current state of the design matrix to the specified file.
     *
     * @param filename The name of the file to which the design matrix will be written.
     *******************************************************************************************/
    void write_fitMatrix(const String& filename) const;
    /*******************************************************************************************
     * @brief Writes the 2D design matrix to a file.
     *
     * This function outputs the contents of the 2D design matrix to a specified file.
     * Each row of the matrix corresponds to a structure, and each column corresponds
     * to a reference conformation. Values are separated by spaces, and rows are
     * separated by newlines.
     *
     * @param filename The name of the file to which the design matrix will be written.
     *
     * @note The file is opened in write mode, and any existing content will be overwritten.
     * @throws std::ios_base::failure if the file cannot be opened for writing.
     *******************************************************************************************/
    void write_fitMatrix2D(const String& filename) const;
    /*******************************************************************************************
     * @brief Retrieves a reference to the pointer of the design matrix data.
     *
     * This function provides access to the underlying pointer of the design matrix data,
     * allowing direct manipulation or inspection of the matrix contents.
     *
     * @return Real*& Reference to the pointer of the design matrix data.
     *******************************************************************************************/
    Real*& get_fitMatrixPtr();
    /*******************************************************************************************
     * @brief Retrieves a pointer to the design matrix data.
     *
     * This function provides access to the underlying data of the design matrix
     * as a constant pointer to `Real`. The returned pointer allows read-only
     * access to the matrix data.
     *
     * @return A constant pointer to the design matrix data.
     *******************************************************************************************/
    const Real* get_fitMatrixPtr() const;
    /*******************************************************************************************
     * @brief Retrieves the target vector for the normal equation.
     *
     * This function provides access to the "normal-equation-target" stored in the data object.
     *
     * @return A constant reference to the Vec1Real object representing the normal equation target.
     *******************************************************************************************/
    const Vec1Real& get_normalEquationTarget() const;
    Size            get_fitMatrixSize0() const;
    Size            get_fitMatrixSize1() const;

  private:
    /*******************************************************************************************
     * @brief Allocates and resizes the outer product matrix for a given key.
     *
     * This function initializes or resizes the outer product matrix associated with the specified
     *key. The matrix is stored as a 1D vector and its size is set to the square of the given size.
     *
     * @param key The unique identifier for the outer product matrix.
     * @param size The dimension of the square matrix (number of rows/columns).
     *******************************************************************************************/
    void allocate_outerProduct(const String& key, const Size size);
    /*******************************************************************************************
     * @brief Allocates and resizes a vector term for force part of the design matrix.
     *
     * This function creates or resizes a vector term identified by the given key
     * to the specified size
     *
     * @param key The unique identifier for the vector term.
     * @param size The desired size of the vector term.
     *******************************************************************************************/
    void allocate_vectorTerm(const String& key, const Size size);
    void allocate_bracketTermDerivative(const KernelPolynomial& kernelPolynomial);
    /*******************************************************************************************
     * @brief Computes and updates the outer product term for a given key and basis.
     *
     * This function calculates the outer product of the input vector `XhatAtom` with itself,
     * scales it by -1.0, and stores the result in the `outerProduct` map associated with the given
     *`key`. Additionally, it modifies the diagonal elements of the resulting matrix by adding 1.0
     *to each.
     *
     * @param key A unique identifier (string) used to store the computed outer product in the
     *`outerProduct` map.
     * @param XhatAtom Pointer to the input vector for which the outer product is computed.
     * @param nBasis The size of the basis (dimension of the input vector and resulting matrix).
     *******************************************************************************************/
    void computeOuterProductTerm(const String& key, const Real* XhatAtom, const Int nBasis);

    ShRec data;
    /*******************************************************************************************
     * @brief Number of structures in the dataset.
     *******************************************************************************************/
    Int nStrucs;
    /*******************************************************************************************
     * @brief Number of basis functions.
     *******************************************************************************************/
    Int nRefConfs;
    /*******************************************************************************************
     * @brief 2D shared memory array representing the design matrix.
     *
     * This matrix is typically used for storing data related to
     * the relationship between structures and reference conformations.
     *******************************************************************************************/
    std::shared_ptr<ShmemArray2D<Real>> fitMatrix;
    /*******************************************************************************************
     * @brief Shared pointer to an MPI wrapper object for parallel processing.
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mpi;
    /*******************************************************************************************
     * @brief Vector storing the number of atoms per structure.
     *******************************************************************************************/
    Vec1Int nAtomsPerStructure;
    /*******************************************************************************************
     * @brief Map storing outer product data, keyed by a string identifier.
     *******************************************************************************************/
    std::map<String, Vec1Real> outerProduct;
    /*******************************************************************************************
     * @brief Map storing contracted outer product term, keyed by a string identifier.
     *******************************************************************************************/
    std::map<String, Vec1Real> vectorTerm;
    /*******************************************************************************************
     * @brief Temporary storage for the design matrix as a 1D vector.
     *******************************************************************************************/
    Vec1Real tempDesignMatrix;
    /*******************************************************************************************
     * Variable to store where force entries for certain structure start in design matrix
     *******************************************************************************************/
    Vec1Int forceOffsets;

    std::map<String, Vec2Real> bracketTermDerivative;
};

/*******************************************************************************************
 * @brief Computes the size of the design matrix based on the number of atoms.
 *
 * This function calculates the total size of the design matrix by applying a
 * transformation and reduction operation on the input vector of atom counts.
 * For each atom count `n`, the transformation `3*n + 7` is applied, and the
 * results are summed up to produce the final size.
 *
 * @param nAtoms A vector containing the number of atoms for each entity.
 * @return The computed size of the design matrix as an integer.
 *******************************************************************************************/
Int computeDesignMatrixSize(const Vec1Int& nAtoms);

} //namespace vaspml

#endif
