#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE mpi_MlMPI

#include "MlMPI.hpp"

#include "boost_helpers.hpp"
#include "FixtureMpi.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_FIXTURE_TEST_CASE(MlMPITest, FixtureMpi )
{
    MlMPI mlmpi( communicator, false );
    REQUIRE_EQUAL_MESSAGE( rank, mlmpi.get_rank(), "MPI rank of MlMPI and FixtureMpi disagree\n" );
    REQUIRE_EQUAL_MESSAGE( numberRanks, mlmpi.get_numberRanks(), "MPI rank of MlMPI and FixtureMpi disagree\n" );
}

BOOST_AUTO_TEST_SUITE_END()
