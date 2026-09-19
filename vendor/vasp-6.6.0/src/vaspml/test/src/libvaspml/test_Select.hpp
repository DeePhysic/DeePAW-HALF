#ifndef TEST_SELECT_HPP
#define TEST_SELECT_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "SampleTrainingData.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{

Real const defaultTolerance = 100.0 * std::numeric_limits<Real>::epsilon();

struct TestCaseSelectOMPL2Norm : public TestCase
{
    TestCaseSelectOMPL2Norm(String name) :
        TestCase(name),
        sample(std::make_shared<SampleTrainingData>(name))
    {}
    std::shared_ptr<SampleTrainingData> sample;
    Vec2Int lrcStructure;
    Vec2Int lrcAtom;
};

template<>
inline void TestCaseContainer<TestCaseSelectOMPL2Norm>::setup()
{
    String relativePath;

    Vec1String versions{"6.5.1"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_AB.CsPbBr3";
        testCases.push_back(TestCaseSelectOMPL2Norm(relativePath));
        TestCaseSelectOMPL2Norm* tc = &(testCases.back());
        tc->lrcStructure = {{1},{1,1},{1}};
        tc->lrcAtom = {{1},{3,2},{5}};
    }

    return;
}

struct TestCaseSelectOMPTreshL1Norm : public TestCase
{
    TestCaseSelectOMPTreshL1Norm(String name) :
        TestCase(name),
        sample(std::make_shared<SampleTrainingData>(name))
    {}
    std::shared_ptr<SampleTrainingData> sample;
    Vec2Int lrcStructure;
    Vec2Int lrcAtom;
};

template<>
inline void TestCaseContainer<TestCaseSelectOMPTreshL1Norm>::setup()
{
    String relativePath;

    Vec1String versions{"6.5.1"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_AB.CsPbBr3";
        testCases.push_back(TestCaseSelectOMPTreshL1Norm(relativePath));
        TestCaseSelectOMPTreshL1Norm* tc = &(testCases.back());
        tc->lrcStructure = {{1},{1,1},{1}};
        tc->lrcAtom = {{1},{3,2},{5}};
    }

    return;
}

} //namespace vaspml

#endif
