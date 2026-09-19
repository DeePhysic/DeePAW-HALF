#ifndef TESTCASEFUNCTION1D_HPP
#define TESTCASEFUNCTION1D_HPP

#include "TestCase.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{

Real const defaultTolerance = 10.0 * std::numeric_limits<Real>::epsilon();

/// Generic test case for 1-dimensional function and its derivative.
struct TestCaseFunction1D : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    /// Tolerance for derivative values.
    Real dtolerance;
    /// x-values for which the function should be computed.
    Vec1Real x;
    /// Expected function values for given x-values.
    Vec1Real fx;
    /// Expected derivatives of function for given x-values.
    Vec1Real dfx;
    /// Additional parameters for this function.
    std::map<String, Real> parameters;

    TestCaseFunction1D(String name);
};

} //namespace vaspml

#endif
