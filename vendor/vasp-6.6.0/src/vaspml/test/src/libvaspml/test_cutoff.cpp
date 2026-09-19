#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE cutoff

#include "test_cutoff.hpp"

#include "boost_helpers.hpp"

#include "SmartEnum.hpp"
#include "cutoff.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;

namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_cutoff> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeCutoff_CorrectValues, bdata::make(container.testCases), testCase)
{
    Vec1Real results(testCase.x.size(), 0.0);
    math::cutoff(testCase.x,
                 results,
                 testCase.parameters.at("rcut"),
                 SmartEnum<math::CutoffType>::toEnum(testCase.cutoffType));
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "cutoff function value");
}

BOOST_DATA_TEST_CASE(ComputeCutoffAndDerivatives_CorrectValues,
                     bdata::make(container.testCases),
                     testCase)
{
    Vec1Real results(testCase.x.size(), 0.0);
    Vec1Real dresults(testCase.x.size(), 0.0);
    math::cutoffAndDerivative(testCase.x,
                              results,
                              dresults,
                              testCase.parameters.at("rcut"),
                              SmartEnum<math::CutoffType>::toEnum(testCase.cutoffType));
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "cutoff function value");
    REQUIRE_CLOSE_COLLECTIONS(dresults, testCase.dfx, testCase.tolerance, "cutoff function deriv");
}

BOOST_AUTO_TEST_SUITE_END()
