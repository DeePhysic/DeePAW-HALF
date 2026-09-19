#ifndef INFERENCE_HPP
#define INFERENCE_HPP

#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "Predictor.hpp"
#include "types.hpp"

namespace vaspml
{

class BasisFunctions;
class MlMPI;
class Structure;
struct Record;

class Inference
{
  public:
    /*******************************************************************************************
     * @brief Constructor for the Inference class.
     *
     * Initializes the Inference object with the provided shared data structure.
     * Sets up internal components including frame, polyKernel, and predictor.
     *
     * @param data Shared data structure (ShRec) used to initialize the inference components.
     *******************************************************************************************/
    Inference(ShRec data = nullptr);
    /*******************************************************************************************
     * @brief Initializes the classic inference model using a machine-learned force field (MLFF)
     *file.
     *
     * This function reads the MLFF data, converts units to atomic units, and sets up the basis
     *functions for 2-body and 3-body interactions. It also initializes the frame, polynomial
     *kernel, and predictor components required for the inference process.
     *
     * @param mlffName The name (or path) of the MLFF file to be loaded.
     *
     * The initialization process includes:
     * - Converting units to atomic units.
     * - Configuring 2-body and 3-body basis functions using radial splines and angular functions.
     * - Initializing the `frame` object with the loaded data and basis functions.
     * - Setting up the polynomial kernel and predictor components.
     *******************************************************************************************/
    void init(const Record& mlff, const std::shared_ptr<MlMPI> mpiIn = nullptr);
    /*******************************************************************************************
     * @brief Performs a prediction based on the given POSCAR file.
     *
     * This function updates the internal frame using the provided POSCAR file name
     * and then executes the prediction process.
     *
     * @param poscarName The name of the POSCAR file to be used for updating the frame.
     *******************************************************************************************/
    void predict(const String& poscarName, const Real& distScale = 1.0);
    /*******************************************************************************************
     * @brief Predicts properties for a list of POSCAR files with a distance scaling factor.
     *
     * Iterates through the provided POSCAR file names, applies the distance scaling factor
     * to update each frame, and performs predictions for the updated frames.
     *
     * @param poscars Collection of POSCAR file names.
     * @param distScale Distance scaling factor for frame updates.
     *******************************************************************************************/
    void predict(const Vec1String& poscarName, const Real& distScale = 1.0);
    /*******************************************************************************************
     * @brief Performs prediction on the given structure.
     *
     * This function updates the internal frame with the provided structure
     * and then executes the prediction process.
     *
     * @param structure Reference to the Structure object to be used for prediction.
     *******************************************************************************************/
    void predict(Structure& structure);
    /*******************************************************************************************
     * @brief Performs predictions on a collection of structures.
     *
     * This function iterates over a vector of `Structure` objects, updates the
     * internal frame with each structure, and performs a prediction for each one.
     *
     * @param structure A vector of `Structure` objects to be processed.
     *******************************************************************************************/
    void predict(std::vector<Structure>& structure);
    /*******************************************************************************************
     * @brief Outputs the current simulation results to the screen, including energy, stress,
     *pressure, and atomic forces.
     *
     * This function writes the following information to the standard output:
     * - Total energy scaled by the provided energy scale.
     * - Stress tensor components scaled by the provided stress scale.
     * - Average pressure derived from the diagonal components of the stress tensor.
     * - Atomic forces for each atom, scaled by a force unit constant.
     *
     * @param energyScale Scaling factor for the total energy.
     * @param forceScale Scaling factor for the atomic forces.
     * @param stressScale Scaling factor for the stress tensor and pressure.
     *
     * @note The function assumes that the `predictor` and `polyKernel` objects are properly
     *initialized and provide access to the required simulation data, such as energy, stress tensor,
     *and atomic forces.
     * @note The atomic forces are scaled using a predefined constant `constants::FUNIT`.
     *******************************************************************************************/
    void  writeCurrentResultToScreen(const Real& energyScale = 1.0,
                                     const Real& forceScale = 1.0,
                                     const Real& stressScale = 1.0) const;
    void  writePredictionsToJson(const String& fname) const;
    ShRec getPredictions() const;

  private:
    /*******************************************************************************************
     * @brief Adds a new predicted data record to the predictions dataset.
     *
     * This function creates a new record in the "predictions" dataset and populates it with
     * various structural and prediction-related data. The data includes system metadata,
     * atomic structure information, and prediction results such as energy, forces, and stress.
     *
     * @param key A string identifier for the system being processed, stored in the "system" field
     *of the record.
     *
     * The following fields are added to the record:
     * - **system**: The identifier of the system (provided by the `key` parameter).
     * - **numTypes**: The number of distinct atom types in the structure.
     * - **numAtomsPerType**: A vector containing the number of atoms for each type.
     * - **numAtoms**: The total number of atoms in the structure.
     * - **typeNames**: A vector of strings representing the names of the atom types.
     * - **lattice**: A vector of real numbers representing the lattice components.
     * - **positions**: A vector of real numbers representing the atomic positions.
     * - **energy**: The total predicted energy of the system.
     * - **forces**: A vector of real numbers representing the predicted atomic forces.
     * - **stress**: A vector of real numbers representing the predicted stress tensor.
     *******************************************************************************************/
    void addPredictedData(const String& key);
    /*******************************************************************************************
     * @brief Perform prediction of atomic forces and stress tensor.
     *
     * This method updates the polynomial kernel with the current frame, retrieves the lattice
     *volume, and uses the predictor to compute atomic forces and the stress tensor based on the
     *updated kernel.
     *******************************************************************************************/
    void predict();

    ShRec                                             data;
    Frame                                             frame;
    std::map<String, std::shared_ptr<BasisFunctions>> basisFunctions;
    KernelPolynomial                                  polyKernel;
    Predictor                                         predictor;
    std::shared_ptr<MlMPI>                            mpi;
};

} // namespace vaspml

#endif
