#ifndef TESTCASECONTAINER_HPP
#define TESTCASECONTAINER_HPP

#include <vector>

namespace vaspml
{

template<class T>
struct TestCaseContainer
{
    TestCaseContainer() { setup(); }

    void setup() { return; }

    std::vector<T> testCases;
};

} //namespace vaspml

#endif
