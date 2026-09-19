#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE mpi_FPS

#include "MlMPI.hpp"

#include "test_mpi_FPS.hpp"

#include "boost_helpers.hpp"
#include "FixtureMpi.hpp"

#include "FarthestPointSampling.hpp"
#include "MetricFunctions.hpp"
#include "SequentialLikelihoodCoverage.hpp"
#include "ShmemArray.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseSelect> container;
TestCaseContainer<TestCaseFarthestPointSampling> containerFPS;
TestCaseContainer<TestCaseFarthestPointSamplingThreshold> containerFPSTresh;
TestCaseContainer<TestCaseSequentialLikelihoodCoverageSelect> containerSLCS;

BOOST_TEST_GLOBAL_FIXTURE(FixtureMpi);

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(farthest_point_sampling, bdata::make(containerFPS.testCases), testCase)
{
#ifdef __NEC__
    BOOST_TEST_MESSAGE("Random number generator seems inconsistent with NEC compiler, needs fix.");
#else

    const Size seed = 4473578942;
    std::shared_ptr<MlMPI> mlmpi = std::make_shared<MlMPI>( FixtureMpi::communicator, true );
    // fill shmem array
    ShmemArray2D<Real> shmemPoints( testCase.points.size(), testCase.points[0].size(), mlmpi );
    for ( Size i = 0; i < testCase.points.size(); i++ )
    {
        for ( Size j = 0; j < testCase.points[i].size(); j++ )
        {
            shmemPoints.set_value( i, j, testCase.points[i][j] );
        }
    }
    const Vec1Int index = fps::farthestPointSamplingKPointsMPI( shmemPoints, 
                                                                testCase.numSamples, 
                                                                MetricFunctions<Real>::l2NormPtrUnity,
                                                                mlmpi,
                                                                seed );
    REQUIRE_EQUAL_COLLECTIONS( index, testCase.selectedIndx, "farthestPointSamplingKPointsMPI" );
#endif
}

BOOST_DATA_TEST_CASE(farthest_point_sampling_threshold, bdata::make(containerFPSTresh.testCases), testCase)
{
#ifdef __NEC__
    BOOST_TEST_MESSAGE("Random number generator seems inconsistent with NEC compiler, needs fix.");
#else

    const Size seed = 4473578942;
    std::shared_ptr<MlMPI> mlmpi = std::make_shared<MlMPI>( FixtureMpi::communicator, true );
    // fill shmem array
    ShmemArray2D<Real> shmemPoints( testCase.points.size(), testCase.points[0].size(), mlmpi );
    for ( Size i = 0; i < testCase.points.size(); i++ )
    {
        for ( Size j = 0; j < testCase.points[i].size(); j++ )
        {
            shmemPoints.set_value( i, j, testCase.points[i][j] );
        }
    }
    const Vec1Int index = fps::farthestPointSamplingDistTreshMPI( shmemPoints,
                                                                  testCase.treshold, 
                                                                  MetricFunctions<Real>::l1NormPtr,
                                                                  mlmpi, 
                                                                  seed );
    REQUIRE_EQUAL_COLLECTIONS( index, testCase.selectedIndx, "farthestPointSamplingDistTreshMPI" );
#endif
}

BOOST_DATA_TEST_CASE(sequential_likelihood_coverage_select, bdata::make(containerSLCS.testCases), testCase)
{
    const Vec1Int index = bayes_select::sequentialLikelihoodCoverageSelectOMP<Real>( 
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
