#include "FixtureMpi.hpp"

#include "boost_helpers.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

Int FixtureMpi::rank = -1;
Int FixtureMpi::numberRanks = -1;
MPI_Comm FixtureMpi::communicator;

void checkMpiReturnCode(const Int info, const String& funcName)
{
    REQUIRE_EQUAL_MESSAGE(info,
                          0,
                          "Attempt to construct/destruct \"MPI\" fixture failed (" + funcName
                              + " returned error code " + std::to_string(info) + ".)\n");
}

FixtureMpi::FixtureMpi() : Fixture( "MPI" )
{
    Int isInitialized = false;
    checkMpiReturnCode(MPI_Initialized(&isInitialized), "MPI_Initialized");
    REQUIRE_EQUAL_MESSAGE(
        isInitialized,
        false,
        "Attempt to construct \"MPI\" fixture failed (MPI_Initialized sets flag isInitialized to "
            + std::to_string(isInitialized) + ").\n");
    checkMpiReturnCode(MPI_Init(nullptr, nullptr), "MPI_Init");
    communicator = MPI_COMM_WORLD;
    checkMpiReturnCode(MPI_Comm_rank(communicator, &rank), "MPI_Comm_rank");
    checkMpiReturnCode(MPI_Comm_size(communicator, &numberRanks), "MPI_Comm_size");
}

FixtureMpi::~FixtureMpi()
{
    checkMpiReturnCode(MPI_Finalize(), "MPI_Finalize");
}
