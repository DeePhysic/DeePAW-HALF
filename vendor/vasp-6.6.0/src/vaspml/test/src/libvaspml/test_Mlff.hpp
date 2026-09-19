#ifndef TEST_MLFF_HPP
#define TEST_MLFF_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "SampleForceField.hpp"

#include "types.hpp"

namespace vaspml
{

struct TestCaseMlff : public TestCase
{
    TestCaseMlff(String name) :
        TestCase(name),
        sample(std::make_shared<SampleForceField>(name))
    {}

    std::shared_ptr<SampleForceField> sample;
    bool hasHeader = true;
};

template<>
inline void TestCaseContainer<TestCaseMlff>::setup()
{
    String relativePath;
    TestCaseMlff* tc = nullptr;

    /*============================================================================================+
     | Default INCAR settings
     +============================================================================================*/
    relativePath = "6.3.2/train/ML_FF.MAPbI3";
    testCases.push_back(TestCaseMlff(relativePath));
    tc = &(testCases.back());
    tc->hasHeader = false;

    Vec1String versions{"6.4.0", "6.4.1", "6.4.3", "6.5.1", "6.6.0"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_FF.MAPbI3";
        testCases.push_back(TestCaseMlff(relativePath));
        tc = &(testCases.back());

        relativePath = v + "/refit/ML_FF.MAPbI3";
        testCases.push_back(TestCaseMlff(relativePath));
        tc = &(testCases.back());
    }

    versions = {"6.5.1", "6.6.0"};
    for (auto const& v : versions)
    {
        relativePath = v + "/refit/ML_FF.CsPbBr3";
        testCases.push_back(TestCaseMlff(relativePath));
        tc = &(testCases.back());
    }

    /*============================================================================================+
     | ML_DESC_TYPE = 1
     +============================================================================================*/
    relativePath = "6.4.3/refit/ML_FF.MAPbI3.DESC_TYPE_1";
    testCases.push_back(TestCaseMlff(relativePath));
    tc = &(testCases.back());

    versions = {"6.5.1", "6.6.0"};
    for (auto const& v : versions)
    {
        relativePath = v + "/train/ML_FF.MAPbI3.DESC_TYPE_1";
        testCases.push_back(TestCaseMlff(relativePath));
        tc = &(testCases.back());

        relativePath = v + "/refit/ML_FF.MAPbI3.DESC_TYPE_1";
        testCases.push_back(TestCaseMlff(relativePath));
        tc = &(testCases.back());
    }

    return;
}

} //namespace vaspml

#endif
