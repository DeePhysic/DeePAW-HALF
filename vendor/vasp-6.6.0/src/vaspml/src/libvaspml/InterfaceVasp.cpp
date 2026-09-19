#include "MlMPI.hpp"

#include "InterfaceVasp.hpp"

#include "ForceFieldGrace.hpp"
#include "ForceFieldKernel.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Tutor.hpp"
#include "io.hpp"
#include "io_detail.hpp"
#include "utils.hpp"

#include <fstream>
#include <sstream>

using namespace vaspml;

MapString vaspml::dataMapInterfaceVasp()
{
    return MapString{};
}

ShRec vaspml::createSetupRecord()
{
    ShRec   result = std::make_shared<Record>();
    Record& setup = result->add("setup", "ShRec").dget<ShRec>("setup");
    setup.add("incar", "ShRec");

    return result;
}

InterfaceVasp::InterfaceVasp() :
    data(assignOrMakeRecord(createSetupRecord(), dataMapInterfaceVasp())),
    setup(data->dget<ShRec>("setup")),
    incar(setup.dget<ShRec>("incar"))
{}

void InterfaceVasp::prepareSetup(const int*    ntyp,
                                 const int*    nsw,
                                 const int*    ibrion,
                                 const int8_t* lmlabexist,
                                 const int8_t* lmlffexist)
{
    // Implicit conversion from int_8_t to bool!
    io::detail::addVaspInterfaceDataToSetup(setup, *ntyp, *nsw, *ibrion, *lmlabexist, *lmlffexist);

    return;
}

void InterfaceVasp::readIncar(const char* incarAsChars)
{
    String            incarString(incarAsChars);
    std::stringstream strm(incarString);
    io::readIncar(incar, strm);

    return;
}

void InterfaceVasp::setupFromIncar(const char* typesAsChars)
{
    Vec1String types = string_tools::split(String(typesAsChars), ";");

    setup.put<Vec1String>("TYPE", types);
    io::setupFromIncar(incar, setup);

    global::parallel.init(setup.cget<String>("ML_PALGO"));

    return;
}

void InterfaceVasp::updateSetup(const int* outblockSync)
{
    Int& outblock = setup.get<Int>("ML_OUTBLOCK");

    if (outblock == *outblockSync) return;

    Int& stateOutblock = setup.dget<ShRec>("states").get<Int>("ML_OUTBLOCK");
    const Int stateIncar = static_cast<Int>(io::TagState::INCAR);
    const Int stateIncarAlt = static_cast<Int>(io::TagState::INCAR_ALT);
    if (stateOutblock == stateIncar || stateOutblock == stateIncarAlt)
    {
        global::tutor.warning(
            "The ML_OUTBLOCK value from the INCAR file (" + std::to_string(outblock)
            + ") was automatically reset to the minimum value (" + std::to_string(*outblockSync)
            + ") over all images used in this run (VCAIMAGES).");
        stateOutblock = static_cast<Int>(io::TagState::OVERRIDE);
    }
    outblock = *outblockSync;

    Int& estblock = setup.get<Int>("ML_ESTBLOCK");

    if (estblock == 0) return;

    Int& stateEstblock = setup.dget<ShRec>("states").get<Int>("ML_ESTBLOCK");

    if (stateEstblock != stateIncar && stateEstblock != stateIncarAlt) estblock = *outblockSync;

    return;
}

void InterfaceVasp::setupMpi(MPI_Fint* mpiComm)
{
    MPI_Comm comm = MPI_Comm_f2c(*mpiComm);
    if (comm == MPI_COMM_NULL)
    {
        global::tutor.warning(
            "Invalid MPI communicator received, continuing without proper MPI setup.");
        mpi = nullptr;
    }
    else
    {
        // TODO: Implement shared memory!
        mpi = std::make_shared<MlMPI>(comm, false);
    }

    if (setup.cget<String>("ML_TYPE") == "grace" && mpi != nullptr && mpi->get_numberRanks() > 1)
    {
        global::tutor.error("GRACE force fields cannot (yet) be combined with MPI parallelization. "
                            "Please run this job only on a single MPI rank.");
    }

    return;
}

void InterfaceVasp::setupForceField()
{
    String type = setup.cget<String>("ML_TYPE");
    if (type == "kernel") { forceField = std::make_unique<ForceFieldKernel>(mpi); }
    else if (type == "grace")
    {
        forceField = std::make_unique<ForceFieldGrace>(setup.cget<String>("ML_GRACE_MODEL"));
    }

    forceField->updateSetup(setup);

    if (mpi == nullptr || mpi->get_rank() == 0)
    {
        std::ofstream logFile;
        io::open(logFile, "ML_LOGFILE", "w");
        io::writeMllogfile(setup, logFile);
        io::close(logFile);
    }

    return;
}

Real InterfaceVasp::get_W1() const
{
    return forceField->get_W1();
}

Real InterfaceVasp::get_W2() const
{
    return forceField->get_W2();
}

Real InterfaceVasp::get_RCUT1() const
{
    return forceField->get_RCUT1();
}

Real InterfaceVasp::get_RCUT2() const
{
    return forceField->get_RCUT2();
}

void InterfaceVasp::set_typeMap(const char* typesAsChars)
{
    Vec1String types = string_tools::split(String(typesAsChars), ";");
    if (string_tools::trim(types.back()).empty()) types.pop_back();

    forceField->set_typeMap(types);

    return;
}

void InterfaceVasp::set_nAtomsType(const int*  ntypes,
                                   const int*  nAtomsType,
                                   const char* keyIn) const
{
    String key(keyIn);
    forceField->set_nAtomsType(*ntypes, nAtomsType, key);

    return;
}

void InterfaceVasp::resizeNeighborArrays(const int* nions, const char* keyIn)
{
    String key(keyIn);
    forceField->resizeNeighborArrays(*nions, key);

    return;
}

void InterfaceVasp::fillNeighborArraysSingleAtom(const int*    numberNeighbors,
                                                 const int*    atomNumber,
                                                 const int*    centralType,
                                                 const int*    neighborIndex,
                                                 const int*    neighborTypes,
                                                 const double* neighborDist,
                                                 const double* neighborConnect,
                                                 const char*   keyIn)
{
    String key(keyIn);
    forceField->fillNeighborArraysSingleAtom(*numberNeighbors,
                                             *atomNumber,
                                             *centralType,
                                             neighborIndex,
                                             neighborTypes,
                                             neighborDist,
                                             neighborConnect,
                                             key);

    return;
}

void InterfaceVasp::fillForceSingleAtom(const int* nions,
                                        const int* ionIndex,
                                        const int* localIndex,
                                        double*    forces) const
{
    forceField->fillForceSingleAtom(*nions, *ionIndex, *localIndex, forces);

    return;
}

void InterfaceVasp::computeForces(const int* nions, double* forces) const
{
    forceField->computeForces(*nions, forces);

    return;
}

void InterfaceVasp::makeGlobalIonIndex(const int* nions)
{
    forceField->makeGlobalIonIndex(*nions);

    return;
}

void InterfaceVasp::update(const double* volume)
{
    forceField->update(*volume);

    return;
}

void InterfaceVasp::getPotentialEnergy(double* energy) const
{
    forceField->getPotentialEnergy(energy);

    return;
}

void InterfaceVasp::getStressTensor(double* stress) const
{
    forceField->getStressTensor(stress);

    return;
}
