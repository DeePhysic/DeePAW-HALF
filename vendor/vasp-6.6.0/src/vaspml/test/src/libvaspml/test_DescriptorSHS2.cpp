#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DescriptorSHS2

#include "test_DescriptorSHS2.hpp"

#include "boost_helpers.hpp"

#include "DescriptorSHS2.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_DescriptorSHS2> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ComputeCLM_and_Derivatives, bdata::make(container.testCases), testCase)
{
    // extract an acronym for the neighbor list
    std::shared_ptr<SampleNeighborList> nn_list = testCase.nn_list;
    // extract an acronym for the descriptor
    std::shared_ptr<BasisFunctionsAngular> descriptor = testCase.descriptor;

    // set up compute clnm calculator
    DescriptorSHS2 coeffCalculator;
    coeffCalculator.set_basisFunctions(descriptor);

    coeffCalculator.updatePairCoefficients(nn_list);

    REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_pair,
                              coeffCalculator.get_clnmPair(0),
                              testCase.tolerance,
                              "clnm_pair_values");

    REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_pair_derivative,
                              coeffCalculator.get_clnmPairDerivative(0),
                              testCase.tolerance,
                              "clnm_pair_derivatives_values");

    coeffCalculator.computeVaspCoefficientsFromPairCoefficients();

    REQUIRE_CLOSE_COLLECTIONS(testCase.clnm_vasp,
                              coeffCalculator.get_clnmVasp(0),
                              testCase.tolerance,
                              "clnm_vasp_values");

    REQUIRE_CLOSE_COLLECTIONS(testCase.derivative_clnm_central_vasp,
                              coeffCalculator.get_clnmDerivativeCentralVasp(0),
                              testCase.tolerance,
                              "clnm_vasp_values");
}

BOOST_AUTO_TEST_SUITE_END()
