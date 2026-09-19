#include "Record.hpp"
#include "io.hpp"
#include "rec.hpp"
#include "types.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    Record       setup;
    Record&      incar = setup.add("incar", "ShRec").dget<ShRec>("incar");
    std::fstream in;

    io::open(in, "INCAR", "r");
    io::readIncar(incar, in);
    std::cout << incar.cget<String>("message");
    io::close(in);

    // Some settings which usually come from VASP.
    setup.put<bool>("LMLABEXIST", std::filesystem::exists("./ML_AB"));
    setup.put<bool>("LMLFFEXIST", std::filesystem::exists("./ML_FF"));
    setup.put<Int>("NSW", 1);
    setup.put<Int>("NTYP", 3);
    setup.put<Vec1String>("TYPE", {"Cs", "Pb", "Br", "F"});
    setup.put<Int>("IBRION", 0);
    // Create complete setup corresponding to ml_reader.F
    io::setupFromIncar(incar, setup);

    std::ofstream log;
    io::open(log, "ML_LOGFILE", "w");
    io::writeMllogfile(setup, log);
    log.close();

    std::cout << rec::printMemoryUsage(rec::computeMemory(setup));

    std::fstream out;
    io::open(out, "setup.json", "w");
    rec::toJson(setup, out);
    io::close(out);

    return 0;
}
