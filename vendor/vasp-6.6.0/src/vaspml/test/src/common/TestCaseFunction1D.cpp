#include "TestCaseFunction1D.hpp"

using namespace vaspml;

TestCaseFunction1D::TestCaseFunction1D(String name) :
    TestCase(name),
    tolerance(defaultTolerance),
    dtolerance(defaultTolerance)
{}
