#ifndef TEST_RECORD_HPP
#define TEST_RECORD_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "record_helpers.hpp"

#include "Record.hpp"
#include "types.hpp"

namespace vaspml
{

struct TestCaseRecord : public TestCase
{
    Record record;

    TestCaseRecord(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseRecord>::setup()
{
    String          testName = "";
    TestCaseRecord* tc = nullptr;

    testName = "Nested 1";
    testCases.push_back(TestCaseRecord(testName));
    tc = &(testCases.back());
    tc->record = generateSampleRecord(testName);

    return;
}

} //namespace vaspml

#endif
