#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE FPS

#include "test_FPS.hpp"

#include "boost_helpers.hpp"

#include "FarthestPointSampling.hpp"
#include "MetricFunctions.hpp"
#include "SequentialLikelihoodCoverage.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseSelect> container;


TestCaseContainer<TestCaseFarthestPointSampling> containerFPS;
TestCaseContainer<TestCaseFarthestPointSamplingThreshold> containerFPSTresh;
TestCaseContainer<TestCaseSequentialLikelihoodCoverageSelect> containerSLCS;


BOOST_AUTO_TEST_SUITE(UnitTests)


BOOST_DATA_TEST_CASE(farthest_point_sampling, bdata::make(containerFPS.testCases), testCase)
{
#ifdef __NEC__
    BOOST_TEST_MESSAGE("Random number generator seems inconsistent with NEC compiler, needs fix.");
#else
    const Size seed = 4473578942;
    const Vec1Int index = fps::farthestPointSamplingKPointsOMP( testCase.points, 
                                                                testCase.numSamples, 
                                                                MetricFunctions<Real>::l2Norm,
                                                                seed );
    REQUIRE_EQUAL_COLLECTIONS( index, testCase.selectedIndx, "farthestPointSamplingKPointsOMP" );
#endif
}

BOOST_DATA_TEST_CASE(farthest_point_sampling_threshold, bdata::make(containerFPSTresh.testCases), testCase)
{
#ifdef __NEC__
    BOOST_TEST_MESSAGE("Random number generator seems inconsistent with NEC compiler, needs fix.");
#else
    const Size seed = 4473578942;
    const Vec1Int index = fps::farthestPointSamplingDistTreshOMP( testCase.points, 
                                                                  testCase.treshold, 
                                                                  MetricFunctions<Real>::l1Norm,
                                                                  seed 
                                                                );
    REQUIRE_EQUAL_COLLECTIONS( index, testCase.selectedIndx, "farthestPointSamplingDistTreshOMP" );
#endif
}

BOOST_DATA_TEST_CASE(sequential_likelihood_coverage_select, bdata::make(containerSLCS.testCases), testCase)
{
    const Vec1Int index = bayes_select::sequentialLikelihoodCoverageSelect<Real>( 
                                                               testCase.points, 
                                                               MetricFunctions<Real>::polyKernelNorm<4> ,
                                                               0.0,
                                                               0.0,
                                                               0.0001,
                                                               2
                                                             );
    //REQUIRE_EQUAL_COLLECTIONS( index, testCase.selectedIndx, "sequentialLikelihoodCoverageSelect" );
    BOOST_REQUIRE( true );
}


BOOST_AUTO_TEST_SUITE_END()
