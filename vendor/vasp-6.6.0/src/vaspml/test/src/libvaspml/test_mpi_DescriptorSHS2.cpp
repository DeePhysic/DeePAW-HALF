#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE mpi_DescriptorSHS2

#include "test_DescriptorSHS2.hpp" // IWYU pragma: associated

#include "boost_helpers.hpp"

#include "Descriptor.hpp"
#include "DescriptorSHS2.hpp"
#include "FixtureMpi.hpp"
#include "MlMPI.hpp"
#include "ParallelEnvironment.hpp"
#include "SampleNeighborList.hpp"
#include "TestCaseContainer.hpp"
#include "nearest_neighbor.hpp"
#include "rec_mpi.hpp"
#include "utils.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_DescriptorSHS2> container;

BOOST_TEST_GLOBAL_FIXTURE(FixtureMpi);

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeCLM_and_Derivatives, bdata::make(container.testCases), testCase)
{
    global::parallel.init("off");

    // extract an acronym for the neighbor list
    std::shared_ptr<SampleNeighborList> nList = testCase.nn_list;
    // extract an acronym for the descriptor
    std::shared_ptr<BasisFunctionsAngular> descriptor = testCase.descriptor;

    std::shared_ptr<MlMPI> mlmpi = std::make_shared<MlMPI>( FixtureMpi::communicator, false );

    Vec1Int roundRobin = mlmpi::getRoundRobinIndexes( nList->get_nAtoms(),
                                                      mlmpi->get_numberRanks(), mlmpi->get_rank() );

    std::shared_ptr<NearestNeighborNSquare> neighborDistList = 
                         std::make_shared<NearestNeighborNSquare>(neighbor_list::getNeighborElements(*nList,roundRobin));

    // set up compute clnm calculator
    DescriptorSHS2 coeffCalculator;
    coeffCalculator.set_basisFunctions(descriptor);

    coeffCalculator.updatePairCoefficients( neighborDistList);
    coeffCalculator.computeVaspCoefficientsFromPairCoefficients();

    ShRec clnmRecord = coeffCalculator.get_dataRecord();
    // collect data on root rank
    ShRec totalclnmRecord = rec::mpi::gatherV( *clnmRecord, *mlmpi );
    Vec1Int roundRobinAll;
    mlmpi->gatherV( roundRobin, roundRobinAll, 0 );
    if ( mlmpi -> get_rank() == 0 )
    {
        DescriptorSHS2 clnmTotal( 1.0, false, DescriptorType::none, totalclnmRecord );
        clnmTotal.set_neighborList( nList );
        Vec1Int newOrder = vector_tools::invertIndex( roundRobinAll );
        clnmTotal.reorderElements( newOrder );
        REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_pair,
                                  coeffCalculator.get_clnmPair(0),
                                  testCase.tolerance,
                                  "clnm_pair_values_mpi");
        REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_pair_derivative,
                                  coeffCalculator.get_clnmPairDerivative(0),
                                  testCase.tolerance,
                                  "clnm_pair_derivatives_values_mpi");
        REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_vasp,
                                  coeffCalculator.get_clnmVasp(0),
                                  testCase.tolerance,
                                  "clnm_vasp_values_mpi");
        REQUIRE_CLOSE_COLLECTIONS(testCase.derivative_clnm_central_vasp,
                              coeffCalculator.get_clnmDerivativeCentralVasp(0),
                              testCase.tolerance,
                              "clnm_vasp_values_mpi");
    }
}

BOOST_AUTO_TEST_SUITE_END()
