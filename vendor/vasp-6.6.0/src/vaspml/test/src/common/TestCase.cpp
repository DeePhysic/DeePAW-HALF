#include "TestCase.hpp"

using namespace vaspml;

TestCase::TestCase(String name) : name(name)
{}

std::ostream& vaspml::operator<<(std::ostream& os, TestCase const& testCase)
{
    os << "Test case name: \"" << testCase.name << "\"";

    return os;
}
