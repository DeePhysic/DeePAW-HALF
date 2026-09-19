#ifndef TEST_STRINGALTERNATIVE_HPP
#define TEST_STRINGALTERNATIVE_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "types.hpp"

namespace vaspml
{

struct TestCaseStringAlternative : public TestCase
{
    Vec2String alternatives;

    TestCaseStringAlternative(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseStringAlternative>::setup()
{
    const Vec1String alt1 = {"main", "alternative1", "alternative2", "alternative3"};
    const Vec1String alt2 = {"MAIN_TAG", "ALT_TAG1", "ALT_TAG2"};
    const Vec1String alt3 = {"Fe", "Fe_2", "Fe_+", "Fe_3", "Fe_--", "Fe_~"};

    String                     testName = "";
    TestCaseStringAlternative* tc = nullptr;

    testName = "Single entry 1";
    testCases.push_back(TestCaseStringAlternative(testName));
    tc = &(testCases.back());
    tc->alternatives.push_back(alt1);

    testName = "Single entry 2";
    testCases.push_back(TestCaseStringAlternative(testName));
    tc = &(testCases.back());
    tc->alternatives.push_back(alt2);

    testName = "Single entry 3";
    testCases.push_back(TestCaseStringAlternative(testName));
    tc = &(testCases.back());
    tc->alternatives.push_back(alt3);

    testName = "Two entries";
    testCases.push_back(TestCaseStringAlternative(testName));
    tc = &(testCases.back());
    tc->alternatives.push_back(alt1);
    tc->alternatives.push_back(alt2);

    testName = "Three entries";
    testCases.push_back(TestCaseStringAlternative(testName));
    tc = &(testCases.back());
    tc->alternatives.push_back(alt1);
    tc->alternatives.push_back(alt2);
    tc->alternatives.push_back(alt3);

    return;
}

} //namespace vaspml

#endif
