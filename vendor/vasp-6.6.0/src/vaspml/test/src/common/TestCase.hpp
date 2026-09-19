#ifndef TESTCASE_HPP
#define TESTCASE_HPP

#include "types.hpp"

#include <ostream>

namespace vaspml
{

struct TestCase
{
    String name;

    TestCase(String name);
};

std::ostream& operator<<(std::ostream& os, vaspml::TestCase const& testCase);

} //namespace vaspml

#endif
