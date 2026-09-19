#ifndef TEST_EXTXYZ_HPP
#define TEST_EXTXYZ_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "types.hpp"

namespace vaspml
{

struct TestCaseExtXyzCommentLine : public TestCase
{
    String line;

    TestCaseExtXyzCommentLine(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseExtXyzCommentLine>::setup()
{
    String          testName = "";
    TestCaseExtXyzCommentLine* tc = nullptr;

    testName = "Test 1";
    testCases.push_back(TestCaseExtXyzCommentLine(testName));
    tc = &(testCases.back());

    return;
}

struct TestCaseExtXyz : public TestCase
{
    TestCaseExtXyz(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseExtXyz>::setup()
{
    String          testName = "";
    TestCaseExtXyz* tc = nullptr;

    testName = "Test 1";
    testCases.push_back(TestCaseExtXyz(testName));
    tc = &(testCases.back());

    return;
}

} //namespace vaspml

#endif
