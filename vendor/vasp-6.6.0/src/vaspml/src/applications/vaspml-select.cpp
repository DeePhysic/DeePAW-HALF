#include "MlMPI.hpp"

#include "Record.hpp"
#include "Selector.hpp"
#include "Timer.hpp"
#include "Tutor.hpp"
#include "io.hpp"
#include "io_detail.hpp"
#include "rec_mpi.hpp"
#include "types.hpp"

#include <fstream>

using namespace vaspml;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        global::tutor.error("No INCAR or ML_AB file supplied. Please Supply an INCAR file "
                            "name and an ML_AB file name");
    }

    const String incarName = argv[1];
    const String mlabName = argv[2];

    ShRec data = std::make_shared<Record>();

    std::shared_ptr<MlMPI> mlmpi = std::make_shared<MlMPI>();
    mlmpi->make_WorldComm();

    if (mlmpi->get_rank() == 0)
    {
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
    }
    rec::mpi::bcast(*data, *mlmpi);

    if (data->get<ShRec>("setup")->get<String>("ML_MODE") == "select")
    {
        data->add("selector", "ShRec");
        Selector selector(data->get<ShRec>("setup"), data->get<ShRec>("mlab"), nullptr, mlmpi);
        selector.select();
        //selector.SelectNumLocRefConfs();
        selector.writeMlabFile("ML_ABC");
    }
    if (mlmpi->get_rank() == 0) { VASPML_PROFILING_WRITE(); }
    MPI_Finalize();

    return 0;
}
