#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE mpi_NeighborList

#include "MlMPI.hpp"

#include "test_NeighborList.hpp" // IWYU pragma: associated

#include "boost_helpers.hpp"

#include "FixtureMpi.hpp"
#include "nearest_neighbor.hpp"
#include "rec_mpi.hpp"
#include "utils.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseNeighborList> container;

BOOST_TEST_GLOBAL_FIXTURE(FixtureMpi);

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeParallelNeighborList,
                     bdata::make(container.testCases),
                     testCase)
{
    // Set up neighbor list according to test case specification.
    NearestNeighborNSquare neighborList(testCase.cutoff, testCase.typeSort, testCase.distSort);

    std::shared_ptr<MlMPI> mlmpi = std::make_shared<MlMPI>( FixtureMpi::communicator, false );

    // Shortcut for test case structure.
    std::shared_ptr<SampleStructure> structure = testCase.structure;
    Vec1Int roundRobin = mlmpi::getRoundRobinIndexes( structure->get_numAtoms(), 
                                                      mlmpi->get_numberRanks(), mlmpi->get_rank() ); 
    Vec1Int blockDist  = mlmpi::getBlockDistributionIndexes( structure->get_numAtoms(),
                                                             mlmpi->get_numberRanks(), mlmpi->get_rank() ); 
    // Tests start with direct coordinates in structure. If necessary, convert now.
    if (!structure->isDirect()) structure->cartesianToDirect();

    Vec1Int roundRobinAll;
    Vec1Int blockDistAll;
    mlmpi->gatherV( roundRobin, roundRobinAll, 0 );
    mlmpi->gatherV( blockDist, blockDistAll, 0 );
    // Compute neighbor list from structure given in direct coordinates.
    BOOST_TEST_CONTEXT("MPI NeighborList from direct coordinates round Robin")
    {
        neighborList.computeNearestNeighborsDirectCoordinates(*structure,roundRobin);
        ShRec neighborData = neighborList.getNeighborRecord();

        ShRec totalData = rec::mpi::allGatherV( *neighborData, *mlmpi );
        //reorder neighbor list
        if ( mlmpi->get_rank() == 0 )
        {

            NearestNeighborNSquare nListTotal( testCase.cutoff, testCase.typeSort, testCase.distSort, totalData );
            Vec1Int newOrder = vector_tools::invertIndex( roundRobinAll );
            nListTotal.reorderNeighborList( newOrder );

            REQUIRE_EQUAL_MESSAGE( nListTotal.get_nAtoms(), testCase.nAtomsTot,
                       "Total number of atoms in MPI neighbor list does not agree.\n" );

            REQUIRE_EQUAL_COLLECTIONS( nListTotal.get_nAtomsType(), testCase.nAtomsPerType,
                           "Total neighbor atoms per type in MPI neighbor list." );

            REQUIRE_EQUAL_COLLECTIONS_2D( testCase.globalIndex,
                                          nListTotal.get_globalIndex(),
                                          "MPI neighborList position array" );

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                         nListTotal.get_distances(),
                                         testCase.tolerance,
                                         "MPI neighborList distances");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                         nListTotal.get_connectionVector(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVector");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                         nListTotal.get_connectionVectorNormalized(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVectorNormalized");
        } 
    }
    BOOST_TEST_CONTEXT("MPI NeighborList from direct coordinates block distributed")
    {
        neighborList.computeNearestNeighborsDirectCoordinates(*structure,blockDist);
        ShRec neighborData = neighborList.getNeighborRecord();

        ShRec totalData = rec::mpi::allGatherV( *neighborData, *mlmpi );
        //reorder neighbor list
        if ( mlmpi->get_rank() == 0 )
        {

            NearestNeighborNSquare nListTotal( testCase.cutoff, testCase.typeSort, testCase.distSort, totalData );
            Vec1Int newOrder = vector_tools::invertIndex( blockDistAll );
            nListTotal.reorderNeighborList( newOrder );

            REQUIRE_EQUAL_MESSAGE( nListTotal.get_nAtoms(), testCase.nAtomsTot,
                       "Total number of atoms in MPI neighbor list does not agree.\n" );

            REQUIRE_EQUAL_COLLECTIONS( nListTotal.get_nAtomsType(), testCase.nAtomsPerType,
                           "Total neighbor atoms per type in MPI neighbor list." );

            REQUIRE_EQUAL_COLLECTIONS_2D( testCase.globalIndex,
                                          nListTotal.get_globalIndex(),
                                          "MPI neighborList position array" );

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                         nListTotal.get_distances(),
                                         testCase.tolerance,
                                         "MPI neighborList distances");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                         nListTotal.get_connectionVector(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVector");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                         nListTotal.get_connectionVectorNormalized(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVectorNormalized");
        } 
    }
    BOOST_TEST_CONTEXT("MPI NeighborList from Cartesian coordinates round Robin")
    {
        structure->directToCartesian();
        neighborList.computeNearestNeighborsCartesianCoordinates(*structure,roundRobin);
        ShRec neighborData = neighborList.getNeighborRecord();

        ShRec totalData = rec::mpi::allGatherV( *neighborData, *mlmpi );
        //reorder neighbor list
        if ( mlmpi->get_rank() == 0 )
        {
            NearestNeighborNSquare nListTotal( testCase.cutoff, testCase.typeSort, testCase.distSort, totalData );
            Vec1Int newOrder = vector_tools::invertIndex( roundRobinAll );
            nListTotal.reorderNeighborList( newOrder );

            REQUIRE_EQUAL_MESSAGE( nListTotal.get_nAtoms(), testCase.nAtomsTot,
                       "Total number of atoms in MPI neighbor list does not agree.\n" );

            REQUIRE_EQUAL_COLLECTIONS( nListTotal.get_nAtomsType(), testCase.nAtomsPerType,
                           "Total neighbor atoms per type in MPI neighbor list." );

            REQUIRE_EQUAL_COLLECTIONS_2D( testCase.globalIndex,
                                          nListTotal.get_globalIndex(),
                                          "MPI neighborList position array" );

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                         nListTotal.get_distances(),
                                         testCase.tolerance,
                                         "MPI neighborList distances");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                         nListTotal.get_connectionVector(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVector");

            REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                         nListTotal.get_connectionVectorNormalized(),
                                         testCase.tolerance,
                                         "MPI neighborList connectionVectorNormalized");
        } 
    }
}

BOOST_AUTO_TEST_SUITE_END()
