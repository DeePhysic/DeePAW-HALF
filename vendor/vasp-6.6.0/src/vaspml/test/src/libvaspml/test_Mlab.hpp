#ifndef TEST_MLAB_HPP
#define TEST_MLAB_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "SampleTrainingData.hpp"

#include "types.hpp"

namespace vaspml
{

struct TestCaseMlab : public TestCase
{
    TestCaseMlab(String name) :
        TestCase(name),
        sample(std::make_shared<SampleTrainingData>(name))
    {}

    std::shared_ptr<SampleTrainingData> sample;
};

template<>
inline void TestCaseContainer<TestCaseMlab>::setup()
{
    String relativePath;

    Vec1String versions{"6.5.1"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_AB.MAPbI3";
        testCases.push_back(TestCaseMlab(relativePath));
    }

    return;
}

} //namespace vaspml

#endif
