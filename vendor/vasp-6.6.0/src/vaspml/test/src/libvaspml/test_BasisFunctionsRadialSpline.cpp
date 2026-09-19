#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE BasisFunctionsRadialSpline

#include "test_BasisFunctionsRadialSpline.hpp"

#include "boost_helpers.hpp"

#include "BasisFunctions.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "cutoff.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;

namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_BasisFunctionsRadialSpline> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(TestCase_BasisFunctionsRadialSpline,
                     bdata::make(container.testCases),
                     testCase)
{

    //DescriptorRadialSpline radial;
    //radial.set_weight( 0.1 );
    //radial.set_cutOffType( 1 );
    //radial.set_cutOff( 5.0 );
    //radial.set_widthBroadening( 0.5 );
    //radial.set_isNormed( true );
    //
    //radial.set_nGrid( 1000 );
    //radial.set_maxOrderBessel( 2 );
    //radial.set_nRootsBessel( 8 );

    //radial.update();
    //

    BasisFunctionsRadialSpline
        radial(math::CutoffType::BP, 5.0, 0.5, 1000, 2, 8, BasisFunctionType::bodyOrder2);

    Vec1Real values;
    Vec1Real derivative;

    Size nBasis = radial.get_totalNumberBasisFunctions();
    values.resize(nBasis);
    derivative.resize(nBasis);

    radial.interpolate(0.5, values, derivative);

    REQUIRE_CLOSE_COLLECTIONS(values, testCase.values, testCase.tolerance, "values");
    REQUIRE_CLOSE_COLLECTIONS(derivative, testCase.derivative, testCase.tolerance, "derviative");
}

BOOST_AUTO_TEST_SUITE_END()
