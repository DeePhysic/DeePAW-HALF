#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE NeighborList

#include "test_NeighborList.hpp"

#include "boost_helpers.hpp"

#include "nearest_neighbor.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseNeighborList> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeNeighborList_CorrectNeighborIndicesAndDistances,
                     bdata::make(container.testCases),
                     testCase)
{
    // Set up neighbor list according to test case specification.
    NearestNeighborNSquare neighborList(testCase.cutoff, testCase.typeSort, testCase.distSort);

    // Shortcut for test case structure.
    std::shared_ptr<SampleStructure> structure = testCase.structure;

    // Tests start with direct coordinates in structure. If necessary, convert now.
    if (!structure->isDirect()) structure->cartesianToDirect();

    // Compute neighbor list from structure given in direct coordinates.
    BOOST_TEST_CONTEXT("NeighborList from direct coordinates")
    {


        neighborList.computeNearestNeighborsDirectCoordinates(*structure);

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     neighborList.get_globalIndex(),
                                     "NeighborList position array");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     neighborList.get_distances(),
                                     testCase.tolerance,
                                     "NeighborList distances");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     neighborList.get_connectionVector(),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     neighborList.get_connectionVectorNormalized(),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
    }

    BOOST_TEST_CONTEXT("NeighborList from Cartesian coordinates")
    {
        // Convert structure to contain positions in Cartesian coordinates.
        structure->directToCartesian();

        // until here it works print and compare
        neighborList.computeNearestNeighborsCartesianCoordinates(*structure);

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     neighborList.get_globalIndex(),
                                     "NeighborList position array");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     neighborList.get_distances(),
                                     testCase.tolerance,
                                     "NeighborList distances");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     neighborList.get_connectionVector(),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");

        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     neighborList.get_connectionVectorNormalized(),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
    }
}

BOOST_AUTO_TEST_SUITE_END()
