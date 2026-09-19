#include "Inference.hpp"
#include "Record.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <iostream>
#include <stdexcept>

using namespace vaspml;

int main(Int argc, char** argv)
{
    bool stdUnit = true;
    if (argc < 3)
    {
        std::cout << str("USAGE: <ML_FF> <POSCAR>\n", argv[0]);
        throw std::runtime_error("ERROR: Wrong number of arguments.");
    }

    String     forceFieldName = argv[1];
    String     poscar_name = argv[2];
    Vec1String poscars;
    if (argc > 3)
    {
        for (Int i = 2; i < argc; i++) { poscars.push_back(argv[i]); }
    }

    Record mlff;
    // TODO: This is only the version WITH unit conversion.
    io::readMlffAndConvertUnits(mlff, forceFieldName);

    Inference model;
    //model.initClassic(forceFieldName);
    model.init(mlff);
    std::cout << "Model initialized " << std::endl;
    Real distScale = 1.0;
    if (stdUnit) distScale = 1.0 / constants::AUTOA;
    if (argc <= 3) model.predict(poscar_name, distScale);
    else model.predict(poscars, distScale);
    std::cout << "Prediction done " << std::endl;

    Real energyScale = 1.0;
    Real forceScale = 1.0;
    Real stressScale = 1.0;
    if (stdUnit)
    {
        energyScale = constants::EUNIT;
        forceScale = constants::FUNIT;
        stressScale = constants::SUNIT;
    }
    if (argc <= 3) model.writeCurrentResultToScreen(energyScale, forceScale, stressScale);

    return 0;
}
