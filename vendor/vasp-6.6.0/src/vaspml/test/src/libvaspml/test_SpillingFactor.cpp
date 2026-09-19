#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE SpillingFactor

#include "test_SpillingFactor.hpp"

#include "boost_helpers.hpp"

#include "BasisFunctions.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "SpillingFactor.hpp"
#include "Structure.hpp"
#include "constants.hpp"
#include "setup.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_SpillingFactor> container;

void spillingFactor(const TestCase_SpillingFactor& testCase,
                    Vec1Real&                      sfAtom,
                    Real&                          min,
                    Real&                          max,
                    Real&                          mean)
{
    const Record& ff = testCase.ff;

    BasisFunctionMap basisFunctions = setup::makeBasisFunctions(ff);

    Frame frame;
    frame.init(ff, basisFunctions);

    ShRec     strucRecord = std::make_shared<Record>();
    Structure struc(strucRecord, false);
    struc.readPoscar(testCase.structure);
    strucRecord->get<Vec1Real>("positions")[0] = testCase.pos0;
    struc.convertUnits(1.0 / constants::AUTOA);
    frame.update(struc);

    KernelPolynomial kernel(ff);
    kernel.updatePolynomialKernel(frame);

    SpillingFactor sf(ff);
    sf.computeSpillingFactor(kernel.get_kernelMatrix(),
                             kernel.get_nAtomsType(),
                             *frame.get_typeMap());
    sf.computeStatistics();
    for (const auto& x : sf.get_spillingFactor())
    {
        for (const auto& y : x) sfAtom.push_back(y);
    }
    min = sf.get_minSpillingFactor();
    max = sf.get_maxSpillingFactor();
    mean = sf.get_averageSpillingFactor();
    return;
}

BOOST_AUTO_TEST_SUITE(IntegrationTests)

BOOST_DATA_TEST_CASE(SpillingFactor, bdata::make(container.testCases), testCase)
{
    Vec1Real sfAtom;
    Real     min;
    Real     max;
    Real     mean;

    spillingFactor(testCase, sfAtom, min, max, mean);

    BOOST_REQUIRE_SMALL(min - testCase.min, testCase.tolerance);
    BOOST_REQUIRE_SMALL(max - testCase.max, testCase.tolerance);
    BOOST_REQUIRE_SMALL(mean - testCase.mean, testCase.tolerance);
    REQUIRE_CLOSE_COLLECTIONS(sfAtom,
                              testCase.sfAtom,
                              testCase.tolerance,
                              "spilling factor per atom");
}

BOOST_AUTO_TEST_SUITE_END()
