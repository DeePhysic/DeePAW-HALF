#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Select

#include "test_Select.hpp"

#include "boost_helpers.hpp"
#include "FixtureRunInDirectory.hpp"

#include "Record.hpp"
#include "Selector.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "io_detail.hpp"
#include "setup.hpp"

#include "boost_unit_test_case.hpp"

#include <sstream>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseSelectOMPL2Norm> containerOMPL2Norm;
TestCaseContainer<TestCaseSelectOMPTreshL1Norm> containerOMPTreshL1Norm;

BOOST_AUTO_TEST_SUITE(UnitTests)

void makeSelection( const ShRec mlab, const String input,
                    Vec2Int& lrcStructure, Vec2Int& lrcAtom )
{
    ShRec setup = std::make_shared<Record>();
    Record&       incar = setup->add("incar", "ShRec").dget<ShRec>("incar");
    std::stringstream ss( input );
    io::readIncar( incar, ss );
    io::detail::addVaspInterfaceDataToSetup(*setup, mlab->cget<Int>("maxTypes"), 1, 0, true);
    io::setupFromIncar(incar, *setup);
    if (setup->get<String>("ML_MODE") == "select")
    {
        Selector selector(setup, mlab, nullptr, nullptr );
        selector.select();
        const ShCRec mlabOut = selector.get_mlab();
        lrcStructure  = mlabOut->cget<Vec2Int>( "lrcStructure" );
        lrcAtom       = mlabOut->cget<Vec2Int>( "lrcAtom" );
    }
}

BOOST_DATA_TEST_CASE_F(FixtureRunInDirectory,
                       SelectModeOMPL2Norm,
                       bdata::make(containerOMPL2Norm.testCases),
                       testCase)

{
    setup(testNameFull() + "_" + testCase.name);
    ShRec mlab = std::make_shared<Record>( testCase.sample->load() );
    setup::rescaleMlabUnits( *mlab,
                             1.0 / constants::EUNIT,
                             1.0 / constants::AUTOA,
                             1.0 / constants::FUNIT,
                             1.0 / constants::SUNIT );
    String input = "ML_LMLFF = .TRUE.; ML_MODE = Select; ML_SALGO = FPSN; ML_SMETRIC = l2norm; ML_SNCONF = 1 2 1;" 
                   "ML_RANDOM_SEED=248489752 0 0;ML_SPAR=openmp";
    Vec2Int lrcStructure;
    Vec2Int lrcAtom;
    makeSelection( mlab, input, lrcStructure, lrcAtom );
    REQUIRE_EQUAL_COLLECTIONS_2D( lrcStructure, testCase.lrcStructure, "TestCaseSelectOMPL2Norm lrcStructure" );
    REQUIRE_EQUAL_COLLECTIONS_2D( lrcAtom, testCase.lrcAtom, "TestCaseSelectOMPL2Norm lrcAtom" );
}

BOOST_DATA_TEST_CASE_F(FixtureRunInDirectory,
                       SelectModeMPIL1Norm,
                       bdata::make(containerOMPTreshL1Norm.testCases),
                       testCase)

{
    setup(testNameFull() + "_" + testCase.name);
    ShRec mlab = std::make_shared<Record>( testCase.sample->load() );
    setup::rescaleMlabUnits( *mlab,
                             1.0 / constants::EUNIT,
                             1.0 / constants::AUTOA,
                             1.0 / constants::FUNIT,
                             1.0 / constants::SUNIT );
    String input = "ML_LMLFF = .TRUE.; ML_MODE = Select; ML_SALGO = FPST; ML_SMETRIC = l1norm; ML_STRESH = 0.1 0.01 0.1;" 
                   "ML_RANDOM_SEED=248489752 0 0;ML_SPAR=openmp";
    Vec2Int lrcStructure;
    Vec2Int lrcAtom;
    makeSelection( mlab, input, lrcStructure, lrcAtom );
    REQUIRE_EQUAL_COLLECTIONS_2D( lrcStructure, testCase.lrcStructure, "TestCaseSelectOMPTreshL1Norm lrcStructure" );
    REQUIRE_EQUAL_COLLECTIONS_2D( lrcAtom, testCase.lrcAtom, "TestCaseSelectOMPTreshL1Norm lrcAtom" );
}

BOOST_AUTO_TEST_SUITE_END()
