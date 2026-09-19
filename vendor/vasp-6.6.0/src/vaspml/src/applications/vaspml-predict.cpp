#include "BasisFunctions.hpp"
#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "ParallelEnvironment.hpp"
#include "Predictor.hpp"
#include "Record.hpp"
#include "Timer.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "io_detail.hpp"
#include "nearest_neighbor.hpp"
#include "settings.hpp"
#include "setup.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace vaspml;

struct ToHeap
{
    BasisFunctionMap basisFunctions;
    Frame            singleStructure;
    KernelPolynomial polyKernel;
    Predictor        predictor;

    ToHeap(const Record& forceFieldData) : polyKernel(forceFieldData), predictor(forceFieldData) {}
};

int main(int argc, char** argv)
{
    bool convert = true;
    if (argc < 4)
    {
        std::cout << str("USAGE: %s <POSCAR> <ML_FF> <ML_PALGO>\n", argv[0]);
        throw std::runtime_error("ERROR: Wrong number of arguments.");
    }
    else if (argc >= 5)
    {
        std::cout << "WARNING: No unit conversion is applied to force field data (experimental).\n";
        convert = false;
    }

    /*============================================================================================+
     | Prepare dummy INCAR to create proper VASPml setup.
     +============================================================================================*/
    ShRec   data = std::make_shared<Record>();
    Record& setup = data->add("setup", "ShRec").dget<ShRec>("setup");
    Record& incar = setup.add("incar", "ShRec").dget<ShRec>("incar");

    String            dummyIncar = "ML_PALGO = " + String(argv[3]);
    std::stringstream strm(dummyIncar);
    io::readIncar(incar, strm);
    io::detail::addVaspInterfaceDataToSetup(setup, 1, 1, 0, false, true);
    setup.put<Vec1String>("TYPE", Vec1String{"?"});
    io::setupFromIncar(incar, setup);

    global::parallel.init(setup.cget<String>("ML_PALGO"));

    /*============================================================================================+
     | Prepare unit conversion factors.
     +============================================================================================*/
    Real energyScale = 1.0;
    Real forceScale = 1.0;
    Real stressScale = 1.0;
    Real distScale = 1.0;

    if (convert)
    {
        energyScale = constants::EUNIT;
        forceScale = constants::FUNIT;
        stressScale = constants::SUNIT;
        distScale = constants::AUTOA;
    }

    String poscar_name = argv[1];
    String forceFieldName = argv[2];

    /*============================================================================================+
     | Setting up force field.
     +============================================================================================*/
    Record forceFieldData;
    if (convert) io::readMlffAndConvertUnits(forceFieldData, forceFieldName);
    else
    {
        std::fstream file;

        io::open(file, forceFieldName, "rb");
        io::processMlff(forceFieldData, file, false);
        io::close(file);
    }

    settings::updateFromMlff(forceFieldData, setup);
    settings::checkFeatureSupport(setup);

    /*============================================================================================+
     | Write log file.
     +============================================================================================*/
    std::fstream out;
    io::open(out, "ML_LOGFILE", "w");
    io::writeMllogfile(setup, out);
    out.close();

    /*============================================================================================+
     | Setup and force field evaluation.
     +============================================================================================*/
    std::shared_ptr<ToHeap> onHeap = std::make_shared<ToHeap>(forceFieldData);
    BasisFunctionMap&       basisFunctions = onHeap->basisFunctions;
    Frame&                  singleStructure = onHeap->singleStructure;
    KernelPolynomial&       polyKernel = onHeap->polyKernel;
    Predictor&              predictor = onHeap->predictor;

    basisFunctions = setup::makeBasisFunctions(forceFieldData);

    singleStructure.init(forceFieldData, basisFunctions);
    singleStructure.update(poscar_name, 1.0 / distScale);

    VASPML_PROFILING_START("update");
    polyKernel.updatePolynomialKernel(singleStructure);
    VASPML_PROFILING_STOP("update");

    predictor.update(polyKernel);

    Real volume = polyKernel.get_descriptorCollection()
                      .getDescriptor("SHS2-2-body")
                      .get_neighborList_ptr()
                      ->get_latticeVolume();
    VASPML_PROFILING_START("Atomic forces");
    predictor.compute_atomicForces(polyKernel.get_descriptorCollection());
    VASPML_PROFILING_STOP("Atomic forces");
    predictor.compute_stressTensor(polyKernel.get_descriptorCollection(), volume);

    /*============================================================================================+
     | Output.
     +============================================================================================*/
    std::cout << str("# Output in units of Angstrom and eV, except stress and pressure in kbar.\n");
    Size numAtoms = polyKernel.get_descriptorCollection()
                        .getDescriptor("SHS2-2-body")
                        .get_neighborList_ptr()
                        ->get_nAtoms();
    std::cout << str("NATOM %zu\n", numAtoms);
    std::cout << str("VOLUME %24.16E\n", volume * std::pow(distScale, 3));
    std::cout << str("ENERGY %24.16E\n", predictor.get_totalEnergy() * energyScale);
    const Vec1Real fxyz = predictor.get_atomicForces();
    for (Size atom = 0; atom < numAtoms; atom++)
    {
        std::cout << str("FORCE %-6zu %24.16E %24.16E %24.16E\n",
                         atom + 1,
                         fxyz[3 * atom] * forceScale,
                         fxyz[3 * atom + 1] * forceScale,
                         fxyz[3 * atom + 2] * forceScale);
    }
    Real pressure = 0.0;
    for (Size i = 0; i < 3; i++)
    {
        std::cout << "STRESS ";
        for (Size j = 0; j < 3; j++)
        {
            std::cout << str(" %24.16E", predictor.get_totalStressTensor(i, j) * stressScale);
        }
        pressure += predictor.get_totalStressTensor(i, i);
        std::cout << "\n";
    }
    pressure *= stressScale / 3.0;
    std::cout << str("PRESSURE %24.16E\n", pressure);

    VASPML_PROFILING_WRITE();

    return 0;
}
