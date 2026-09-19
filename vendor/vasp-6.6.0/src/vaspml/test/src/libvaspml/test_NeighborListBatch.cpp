#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE NeighborListBatch

#include "test_NeighborListBatch.hpp"

#include "boost_helpers.hpp"

#include "Structure.hpp"
#include "BatchMap.hpp"
#include "nearest_neighbor.hpp"
#include "Record.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseNeighborListBatch> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeNeighborListBatch,
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
    BOOST_TEST_CONTEXT("NeighborListBatch from direct coordinates")
    {
        std::vector<Structure> strucVector;
        strucVector.push_back( *structure );
        strucVector.push_back( *structure );
        AtomBatchMap batchMap;
        Vec1String typeOrder = { "Pb", "Br", "Cs" };
        batchMap.makeMap( strucVector, typeOrder );

        neighborList.computeNeighborsDirectCoordinatesBatch( strucVector.cbegin(), strucVector.cend(), batchMap );
        REQUIRE_EQUAL_MESSAGE( neighborList.get_nAtoms(), testCase.nAtomsTot, 
                               "Total number of atoms in batched neighbor list does not agree.\n" );
        REQUIRE_EQUAL_COLLECTIONS( neighborList.get_nAtomsType(), testCase.nAtomsPerType, 
                                   "Total neighbor atoms per type in bacthed neighbor list." );
        ShRec list1 = neighbor_list::extractNeighborListFromBatchList( neighborList, batchMap, 0 ); 

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     list1->get<Vec2Int>( "globalIndex" ),
                                     "NeighborList position array");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     list1->get<Vec2Real>( "distances" ),
                                     testCase.tolerance,
                                     "NeighborList distances");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     list1->get<Vec2Real>( "connectionVector" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     list1->get<Vec2Real>( "connectionVectorNormalized" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
        
        // compare second part of list
        ShRec list2 = neighbor_list::extractNeighborListFromBatchList( neighborList, batchMap, 0 ); 

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     list2->get<Vec2Int>( "globalIndex" ),
                                     "NeighborList position array");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     list2->get<Vec2Real>( "distances" ),
                                     testCase.tolerance,
                                     "NeighborList distances");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     list2->get<Vec2Real>( "connectionVector" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     list2->get<Vec2Real>( "connectionVectorNormalized" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
    }

    BOOST_TEST_CONTEXT("NeighborList from Cartesian coordinates")
    {
        // Convert structure to contain positions in Cartesian coordinates.
        structure->directToCartesian();
        
        std::vector<Structure> strucVector;
        strucVector.push_back( *structure );
        strucVector.push_back( *structure );
        AtomBatchMap batchMap;
        Vec1String typeOrder = { "Pb", "Br", "Cs" };
        batchMap.makeMap( strucVector, typeOrder );

        // until here it works print and compare
        neighborList.computeNeighborsCartesianCoordinatesBatch( strucVector.cbegin(), strucVector.cend(), batchMap );
        
        REQUIRE_EQUAL_MESSAGE( neighborList.get_nAtoms(), testCase.nAtomsTot, 
                               "Total number of atoms in batched neighbor list does not agree\n" );
        REQUIRE_EQUAL_COLLECTIONS( neighborList.get_nAtomsType(), testCase.nAtomsPerType, 
                                   "Total neighbor atoms per type in bacthed neighbor list" );
        ShRec list1 = neighbor_list::extractNeighborListFromBatchList( neighborList, batchMap, 0 );

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     list1->get<Vec2Int>( "globalIndex" ),
                                     "NeighborList position array");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     list1->get<Vec2Real>( "distances" ),
                                     testCase.tolerance,
                                     "NeighborList distances");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     list1->get<Vec2Real>( "connectionVector" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     list1->get<Vec2Real>( "connectionVectorNormalized" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
        
        // compare second part of list
        ShRec list2 = neighbor_list::extractNeighborListFromBatchList( neighborList, batchMap, 0 ); 

        REQUIRE_EQUAL_COLLECTIONS_2D(testCase.globalIndex,
                                     list2->get<Vec2Int>( "globalIndex" ),
                                     "NeighborList position array");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.distances,
                                     list2->get<Vec2Real>( "distances" ),
                                     testCase.tolerance,
                                     "NeighborList distances");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVector,
                                     list2->get<Vec2Real>( "connectionVector" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVector");
        REQUIRE_CLOSE_COLLECTIONS_2D(testCase.connectionVectorNormalized,
                                     list2->get<Vec2Real>( "connectionVectorNormalized" ),
                                     testCase.tolerance,
                                     "NeighborList connectionVectorNormalized");
    }
}

BOOST_AUTO_TEST_SUITE_END()
