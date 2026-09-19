#include "GraceInterface.hpp"
#include "Lattice.hpp"
#include "Record.hpp"
#include "Structure.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

using namespace vaspml;

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
#ifdef VASPML_ENABLE_GRACE
    if (argc < 3)
    {
        std::cout << str("USAGE: %s <POSCAR> <MODEL> [<PADDING>]\n", argv[0]);
        throw std::runtime_error("ERROR: Wrong number of arguments.");
    }
    String poscarName = argv[1];
    String modelPath = argv[2];
    Real   padding = 1;
    if (argc >= 4) { padding = std::atof(argv[3]); }

    GraceInterface grace(modelPath, padding);
    Real           cutoff = grace.getMaximumCutoff() / constants::AUTOA;

    Structure structure;
    structure.readPoscar(poscarName);
    structure.convertUnits(1.0 / constants::AUTOA);
    if (structure.isDirect()) structure.directToCartesian();

    /*============================================================================================+
     | Set up a neighbor list instance.
     |
     | Here, we create a Record externally ("neighborListData") which will hold the actual neighbor
     | list data. In this way can directly access the memory and write a JSON file later on.
     | If the last constructor argument is ommitted the Record is set up internally and hence
     | private to the neighbor list instance.
     +============================================================================================*/
    ShRec                  neighborListData = std::make_shared<Record>();
    NearestNeighborNSquare neighborList(cutoff, true, false, neighborListData);

    /*============================================================================================+
     | Compute neighbor list for the given structure in POSCAR file.
     +============================================================================================*/
    neighborList.computeNearestNeighborsCartesianCoordinates(structure);
    grace.updateNeighborList(structure, neighborList);

    /*============================================================================================+
     | GRACE prediction.
     +============================================================================================*/
    grace.update();

    /*============================================================================================+
     | Output.
     +============================================================================================*/
    std::cout << str("# Output in units of Angstrom and eV, except stress and pressure in kbar.\n");
    Size nAtoms = neighborList.get_nAtoms();
    std::cout << str("NATOM %zu\n", nAtoms);
    const double volume = structure.get_lattice().get_volume() * std::pow(constants::AUTOA, 3);
    std::cout << str("VOLUME %24.16E\n", volume);
    std::cout << str("ENERGY %24.16E\n", grace.getEnergy());
    for (auto e : grace.getEnergyPerAtom()) { std::cout << str("EATOM %24.16E\n", e); }
    Vec1Real forces = grace.getForces();
    for (Size atom = 0; atom < nAtoms; ++atom)
    {
        std::cout << str("FORCE %-6zu %24.16E %24.16E %24.16E\n",
                         atom + 1,
                         forces[3 * atom],
                         forces[3 * atom + 1],
                         forces[3 * atom + 2]);
    }
    Vec1Real stress6 = grace.getStressTensor();
    Vec2Real stress3x3(3);
    for (auto& v : stress3x3) v.resize(3);
    const size_t XX = 0;
    const size_t YY = 1;
    const size_t ZZ = 2;
    const size_t XY = 3;
    const size_t YX = 3;
    const size_t XZ = 4;
    const size_t ZX = 4;
    const size_t YZ = 5;
    const size_t ZY = 5;
    const double factor = constants::EVTOJ * 1.0E22 / volume;
    stress3x3[0][0] = factor * stress6[XX];
    stress3x3[0][1] = factor * stress6[XY];
    stress3x3[0][2] = factor * stress6[XZ];
    stress3x3[1][0] = factor * stress6[YX];
    stress3x3[1][1] = factor * stress6[YY];
    stress3x3[1][2] = factor * stress6[YZ];
    stress3x3[2][0] = factor * stress6[ZX];
    stress3x3[2][1] = factor * stress6[ZY];
    stress3x3[2][2] = factor * stress6[ZZ];
    Real pressure = 0.0;
    for (Size i = 0; i < 3; i++)
    {
        std::cout << "STRESS ";
        for (Size j = 0; j < 3; j++) { std::cout << str(" %24.16E", stress3x3[i][j]); }
        pressure += stress3x3[i][i];
        std::cout << "\n";
    }
    pressure /= 3.0;
    std::cout << str("PRESSURE %24.16E\n", pressure);

    return 0;
#else
    std::cout << "ERROR: This application requires VASPml to be compiled with GRACE support "
                 "(-DVASPML_ENABLE_GRACE).\n";
    return 1;
#endif
}
