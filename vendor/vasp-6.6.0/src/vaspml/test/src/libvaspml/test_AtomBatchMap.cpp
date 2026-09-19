#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE AtomBatchMap

#include "test_AtomBatchMap.hpp"

#include "boost_helpers.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_AtomBatchMap>  container_AtomBatchMap;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(AtomBatchMapTest,
                     bdata::make(container_AtomBatchMap.testCases),
                     testCase)
{

    REQUIRE_EQUAL_COLLECTIONS( testCase.map.get_types(), testCase.initOrder, "original type order" );
    Vec1Int order;
    for ( Int i =0; i < (Int)testCase.map.get_typesTO().size(); i++) 
                   order.push_back( testCase.map.get_mapTO_OrigOrder(i) );
    REQUIRE_EQUAL_COLLECTIONS( order, testCase.typeOrderInt, 
                               "index reorder" );
    REQUIRE_EQUAL_COLLECTIONS( testCase.map.get_typesTO(), testCase.typeOrder, 
                               "types in type order" );


    Vec1String types;
    Vec1String typesTO;
    for ( Size type = 0; type < testCase.strucList.size(); type++ )
    {
        for ( Size locRef = 0; locRef < testCase.atomList[type].size(); locRef++ )
        {
            const Size& nStruc = testCase.strucList[type][locRef];
            const Size& atom   = testCase.atomList[type][locRef];
            const Size& idx    = testCase.map.get_mapStrucAtom_Batch( nStruc, atom );
            const Size& idxTO  = testCase.map.get_mapStrucAtom_BatchTO( nStruc, atom );
            types.push_back( testCase.map.get_types( idx ) );
            typesTO.push_back( testCase.map.get_typesTO( idxTO ) );
        }
    }
    REQUIRE_EQUAL_COLLECTIONS( types, testCase.typeBatchTest, 
                               "original order batch test" );
    REQUIRE_EQUAL_COLLECTIONS( typesTO, testCase.typeBatchTest, 
                               "type order order batch test" );

}


BOOST_AUTO_TEST_SUITE_END()
