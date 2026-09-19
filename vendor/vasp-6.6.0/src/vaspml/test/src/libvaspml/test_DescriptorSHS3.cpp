#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DescriptorSHS3

#include "test_DescriptorSHS3.hpp"

#include "boost_helpers.hpp"

#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_DescriptorSHS3_CsPbBr3> containerCsPbBr3;
TestCaseContainer<TestCase_DescriptorSHS3_CaO>     containerCaO;
TestCaseContainer<TestCase_DescriptorSHS3_HCN>     containerHCN;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(Compute_DescriptorSHS3_CsPbBr3,
                     bdata::make(containerCsPbBr3.testCases),
                     testCase)
{
    // extract an acronym for the neighbor list
    std::shared_ptr<SampleNeighborList> nn_list = testCase.nn_list;
    // extract an acronym for the descriptor
    std::shared_ptr<BasisFunctionsAngular> descriptor = testCase.descriptor;

    // set up compute clnm calculator
    std::shared_ptr<DescriptorSHS2> coeffCalculator = std::make_shared<DescriptorSHS2>();
    coeffCalculator->set_basisFunctions(descriptor);

    coeffCalculator->updatePairCoefficients(nn_list);
    coeffCalculator->computeVaspCoefficientsFromPairCoefficients();

    DescriptorSHS3 descriptor3body(testCase.descriptorListMLFF,
                                   true,
                                   2,
                                   0.002,
                                   coeffCalculator->get_nRoots(),
                                   coeffCalculator->get_maxOrder());

    descriptor3body.computeSHS3(coeffCalculator, testCase.typeMap);

    REQUIRE_CLOSE_COLLECTIONS(descriptor3body.get_SHS3Atom_vasp(0),
                              testCase.targetSHS3,
                              testCase.tolerance,
                              "SHSPbBrCs_atom0");
}

BOOST_DATA_TEST_CASE(Compute_DescriptorSHS3_CaO, bdata::make(containerCaO.testCases), testCase)
{
    // extract an acronym for the neighbor list
    std::shared_ptr<SampleNeighborList> nn_list = testCase.nn_list;
    // extract an acronym for the descriptor
    std::shared_ptr<BasisFunctionsAngular> descriptor = testCase.descriptor;

    // set up compute clnm calculator
    std::shared_ptr<DescriptorSHS2> coeffCalculator = std::make_shared<DescriptorSHS2>();
    coeffCalculator->set_basisFunctions(descriptor);

    coeffCalculator->updatePairCoefficients(nn_list);
    coeffCalculator->computeVaspCoefficientsFromPairCoefficients();

    DescriptorSHS3 descriptor3body(testCase.descriptorListMLFF,
                                   true,
                                   2,
                                   0.002,
                                   coeffCalculator->get_nRoots(),
                                   coeffCalculator->get_maxOrder());

    descriptor3body.computeSHS3(coeffCalculator, testCase.typeMap);

    REQUIRE_CLOSE_COLLECTIONS(descriptor3body.get_SHS3Atom_vasp(0),
                              testCase.targetSHS3,
                              testCase.tolerance,
                              "CaO_atom0");
}

BOOST_DATA_TEST_CASE(Compute_DescriptorSHS3_HCN, bdata::make(containerHCN.testCases), testCase)
{
    // extract an acronym for the neighbor list
    std::shared_ptr<SampleNeighborList> nn_list = testCase.nn_list;
    // extract an acronym for the descriptor
    std::shared_ptr<BasisFunctionsAngular> descriptor = testCase.descriptor;

    // set up compute clnm calculator
    std::shared_ptr<DescriptorSHS2> coeffCalculator = std::make_shared<DescriptorSHS2>();
    coeffCalculator->set_basisFunctions(descriptor);

    coeffCalculator->updatePairCoefficients(nn_list);
    coeffCalculator->computeVaspCoefficientsFromPairCoefficients();

    DescriptorSHS3 descriptor3body(testCase.descriptorListMLFF,
                                   true,
                                   2,
                                   0.002,
                                   coeffCalculator->get_nRoots(),
                                   coeffCalculator->get_maxOrder());

    descriptor3body.computeSHS3(coeffCalculator, testCase.typeMap);

    REQUIRE_CLOSE_COLLECTIONS(descriptor3body.get_SHS3Atom_vasp(0),
                              testCase.targetSHS3,
                              testCase.tolerance,
                              "HCN_atom0");
}

BOOST_AUTO_TEST_SUITE_END()
