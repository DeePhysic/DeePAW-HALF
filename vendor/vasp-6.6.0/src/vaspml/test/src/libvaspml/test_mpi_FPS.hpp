#ifndef TEST_MPI_FPS_HPP
#define TEST_MPI_FPS_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "SampleTrainingData.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{


Real const defaultTolerance = 100.0 * std::numeric_limits<Real>::epsilon();

struct TestCaseFarthestPointSampling : public TestCase
{
    Vec2Real points;
    Int numSamples;
    Vec1Int selectedIndx;
    TestCaseFarthestPointSampling(String name) :
        TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseFarthestPointSampling>::setup()
{
    
    TestCaseFarthestPointSampling* tc = nullptr;
    testCases.push_back(TestCaseFarthestPointSampling("TestCaseFarthestPointSamplingMPI"));
    tc = &(testCases.back());

    tc->points.push_back({0.0, 0.0});
    tc->points.push_back({1.0, 1.0});
    tc->points.push_back({2.0, 2.0});
    tc->points.push_back({3.0, 3.0});
    tc->numSamples = 2;

    tc->selectedIndx.push_back( 3 );
    tc->selectedIndx.push_back( 0 );

    return;
}

struct TestCaseFarthestPointSamplingThreshold : public TestCase
{
    Vec2Real points;
    Real treshold;
    Vec1Int selectedIndx;
    TestCaseFarthestPointSamplingThreshold( String name ) :
        TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseFarthestPointSamplingThreshold>::setup()
{
    TestCaseFarthestPointSamplingThreshold* tc = nullptr;
    testCases.push_back(TestCaseFarthestPointSamplingThreshold("TestCaseFarthestPointSamplingThreshold"));
    tc = &(testCases.back());

    tc->points.push_back({0.0, 0.0});
    tc->points.push_back({1.0, 1.0});
    tc->points.push_back({2.0, 2.0});
    tc->points.push_back({3.0, 3.0});
    tc->treshold = 5.0;

    tc->selectedIndx.push_back( 3 );
    tc->selectedIndx.push_back( 0 );
    return;
}

struct TestCaseSequentialLikelihoodCoverageSelect : public TestCase
{
    Vec2Real points;
    Vec1Int selectedIndx;
    Int maxSamples;
    TestCaseSequentialLikelihoodCoverageSelect( String name ) :
        TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseSequentialLikelihoodCoverageSelect>::setup()
{
    TestCaseSequentialLikelihoodCoverageSelect* tc = nullptr;
    testCases.push_back(TestCaseSequentialLikelihoodCoverageSelect("TestCaseSequentialLikelihoodCoverageSelect"));
    tc = &(testCases.back());

    tc->points.push_back({0.0, 0.0});
    tc->points.push_back({1.0, 1.0});
    tc->points.push_back({2.0, 2.0});
    tc->points.push_back({3.0, 3.0});

    tc->maxSamples = 2;

    tc->selectedIndx.push_back( 2 );
    tc->selectedIndx.push_back( 3 );
    return;
}

struct TestCaseSelect : public TestCase
{
    TestCaseSelect(String name) :
        TestCase(name),
        sample(std::make_shared<SampleTrainingData>(name))
    {}

    std::shared_ptr<SampleTrainingData> sample;
};

template<>
inline void TestCaseContainer<TestCaseSelect>::setup()
{
    String relativePath;

    Vec1String versions{"6.5.1"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_AB.CsPbBr3";
        testCases.push_back(TestCaseSelect(relativePath));
    }

    return;
}

} //namespace vaspml

#endif
