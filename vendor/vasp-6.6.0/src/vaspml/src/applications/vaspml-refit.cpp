#include "Record.hpp"
#include "Refit.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "io_detail.hpp"
#include "setup.hpp"
#include "types.hpp"

#include <fstream>
#include <iostream>

using namespace vaspml;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        global::tutor.error("NO INCAR or ML_AB file supplied. Please Supply an INCAR file "
                            "name and an ML_AB file name");
    }
    const String incarName = argv[1];
    const String mlabName = argv[2];

    ShRec data = std::make_shared<Record>();

    // Read data set from ML_AB.
    Record& mlab = data->add("mlab", "ShRec").dget<ShRec>("mlab");
    io::readMlab(mlab, mlabName);

    // Create setup from INCAR file.
    Record&       setup = data->add("setup", "ShRec").dget<ShRec>("setup");
    Record&       incar = setup.add("incar", "ShRec").dget<ShRec>("incar");
    std::ifstream incarFile;
    io::open(incarFile, incarName, "r");
    io::readIncar(incar, incarFile);
    io::close(incarFile);
    io::detail::addVaspInterfaceDataToSetup(setup, mlab.cget<Int>("maxTypes"), 1, 0, true);
    io::setupFromIncar(incar, setup);

    if (argc > 3)
    {
        std::cout << "CONVERTING TO ATOMIC UNITS" << std::endl;
        setup::rescaleIncarUnits(setup, 1.0 / constants::AUTOA);
        setup::rescaleMlabUnits(mlab,
                                1.0 / constants::EUNIT,
                                1.0 / constants::AUTOA,
                                1.0 / constants::FUNIT,
                                1.0 / constants::SUNIT);
    }

    std::cout << "START Refit construct" << std::endl;
    Refit fitting(data->get<ShRec>("setup"), data->get<ShRec>("mlab"));
    std::cout << "Refit construct done" << std::endl;

    std::cout << "Precompute data for fitting" << std::endl;
    fitting.prepareFit();
    std::cout << "Precompute data for fitting done" << std::endl;

    std::cout << "Starting SVD " << std::endl;
    fitting.computeFit();
    fitting.prepare_mlffRecord();
    fitting.write_mlffFile("ML_FF");
    //fitting.SVDFittingTest();
    std::cout << "SVD done " << std::endl;
    std::cout << "Check errors from fitting " << std::endl;

    //fitting.trainDataInference();
    std::cout << "Error checking done" << std::endl;
    fitting.designMatrixTest();

    return 0;
}
