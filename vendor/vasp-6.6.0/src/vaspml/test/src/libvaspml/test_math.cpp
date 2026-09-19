#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE math

#include "test_math.hpp"

#include "boost_helpers.hpp"

#include "math.hpp"

#include "boost_unit_test_case.hpp"

#include <algorithm>
#include <cmath>

using namespace vaspml;

namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_math_clebschGordan>           containerClebschGordan;
TestCaseContainer<TestCase_math_clebschGordanM0>         containerClebschGordanM0;
TestCaseContainer<TestCase_math_gaussian>                containerGaussian;
TestCaseContainer<TestCase_math_sphericalBessel>         containerSphericalBessel;
TestCaseContainer<TestCase_math_sphericalBesselRoots>    containerSphericalBesselRoots;
TestCaseContainer<TestCase_math_modifiedSphericalBessel> containerModifiedSphericalBessel;
TestCaseContainer<TestCase_math_expModSphBessel>         containerExpModSphBessel;
TestCaseContainer<TestCase_math_calculateMean>           containerMean;
TestCaseContainer<TestCase_math_calculateVariance>       containerVariance;
TestCaseContainer<TestCase_math_calculateRMSE>           containerRMSE;
TestCaseContainer<TestCase_math_calculateMAE>            containerMAE;
TestCaseContainer<TestCase_math_pearsonCorr>             containerPearson;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(ComputeFactorials_CorrectVales)
{
    // Compute factorials up to 10!!, i.e. size needs to be (10 + 1) = 11.
    Vec1Real f(11);
    math::factorial(f);
    // clang-format off
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 0])), 1);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 1])), 1);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 2])), 2);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 3])), 6);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 4])), 24);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 5])), 120);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 6])), 720);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 7])), 5040);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 8])), 40320);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 9])), 362880);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[10])), 3628800);
    // clang-format on
}

BOOST_AUTO_TEST_CASE(ComputeDoubleFactorialsOdd_CorrectVales)
{
    // Compute Real factorials up to 13!!, i.e. size needs to be (13 + 1) / 2 = 7.
    Vec1Real f(7);
    math::doubleFactorialOdd(f);
    // clang-format off
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 0])), 1);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 1])), 3);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 2])), 15);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 3])), 105);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 4])), 945);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 5])), 10395);
    BOOST_REQUIRE_EQUAL(static_cast<UInt>(std::round(f[ 6])), 135135);
    // clang-format on
}

BOOST_DATA_TEST_CASE(ComputeClebschGordan_CorrectCoefficients,
                     bdata::make(containerClebschGordan.testCases),
                     testCase)
{
    Vec1Real results;
    Vec1Real ref;
    Vec1Real factorials(4 * testCase.jmax + 2);
    math::factorial(factorials);
    for (auto const& cg : testCase.cgs)
    {
        ref.push_back(cg.coeff);
        results.push_back(math::clebschGordan(cg.j1, cg.j2, cg.J, cg.m1, cg.m2, cg.M, factorials));
    }
    REQUIRE_CLOSE_COLLECTIONS(results, ref, testCase.tolerance, "Clebsch-Gordan coefficients");
}

BOOST_DATA_TEST_CASE(ComputeClebschGordanM0_CorrectCoefficients,
                     bdata::make(containerClebschGordanM0.testCases),
                     testCase)
{
    Vec1Real results;
    Vec1Real ref;
    Vec1Real factorials(4 * testCase.jmax + 2);
    math::factorial(factorials);
    for (auto const& cg : testCase.cgs)
    {
        ref.push_back(cg.coeff);
        results.push_back(math::clebschGordanM0(cg.j1, cg.j2, cg.J, factorials));
    }
    REQUIRE_CLOSE_COLLECTIONS(results, ref, testCase.tolerance, "Clebsch-Gordan coefficients (M0)");
}

BOOST_DATA_TEST_CASE(ComputeGaussians_CorrectValues,
                     bdata::make(containerGaussian.testCases),
                     testCase)
{
    Vec1Real results(testCase.x.size(), 0.0);
    math::gaussian(testCase.x, results, testCase.parameters.at("width"));
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "Gaussian");
}

BOOST_DATA_TEST_CASE(ComputeGaussiansAndDerivatives_CorrectValues,
                     bdata::make(containerGaussian.testCases),
                     testCase)
{
    Size     n = testCase.x.size();
    Vec1Real results(n, 0.0);
    Vec1Real dresults(n, 0.0);
    math::gaussianAndDerivative(testCase.x, results, dresults, testCase.parameters.at("width"));
    REQUIRE_CLOSE_COLLECTIONS(dresults, testCase.dfx, testCase.dtolerance, "Gaussian derivatives");
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "Gaussian");
}

BOOST_DATA_TEST_CASE(ComputeSphericalBessel_CorrectValues,
                     bdata::make(containerSphericalBessel.testCases),
                     testCase)
{
    Vec1Real results(testCase.x.size(), 0.0);
    math::sphericalBessel(testCase.x, results, testCase.parameters.at("nu"));
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "Spherical Bessel");
}

BOOST_DATA_TEST_CASE(ComputeSphericalBesselAndDerivative_CorrectValues,
                     bdata::make(containerSphericalBessel.testCases),
                     testCase)
{
    Size     n = testCase.x.size();
    Vec1Real results(n, 0.0);
    Vec1Real dresults(n, 0.0);
    math::sphericalBesselAndDerivative(testCase.x, results, dresults, testCase.parameters.at("nu"));
    REQUIRE_CLOSE_COLLECTIONS(results, testCase.fx, testCase.tolerance, "Spherical Bessel");
    REQUIRE_CLOSE_COLLECTIONS(dresults,
                              testCase.dfx,
                              testCase.dtolerance,
                              "Spherical Bessel derivatives");
}

BOOST_DATA_TEST_CASE(ComputeSphericalBesselRoots_CorrectRoots,
                     bdata::make(containerSphericalBesselRoots.testCases),
                     testCase)
{
    Vec1Real results(testCase.roots.size(), 0.0);
    math::sphericalBesselRoots(results, testCase.order);
    REQUIRE_CLOSE_COLLECTIONS(results,
                              testCase.roots,
                              testCase.tolerance,
                              "Spherical Bessel roots");
}

BOOST_DATA_TEST_CASE(ComputeModifiedSphericalBesselAndDerivative_CorrectValues,
                     bdata::make(containerModifiedSphericalBessel.testCases),
                     testCase)
{
    Vec2Real results(testCase.fx.size());
    for (auto& ri : results) ri.resize(testCase.x.size());
    Vec2Real dresults(testCase.fx.size());
    for (auto& dri : dresults) dri.resize(testCase.x.size());
    math::modifiedSphericalBesselAndDerivative(testCase.x, results, dresults);
    REQUIRE_CLOSE_COLLECTIONS_2D(results,
                                 testCase.fx,
                                 testCase.tolerance,
                                 "Modified Spherical Bessel");
    REQUIRE_CLOSE_COLLECTIONS_2D(dresults,
                                 testCase.dfx,
                                 testCase.tolerance,
                                 "Modified Spherical Bessel derivatives");
}

BOOST_DATA_TEST_CASE(ComputeExpModSphBessel_CorrectValues,
                     bdata::make(containerExpModSphBessel.testCases),
                     testCase)
{
    Vec2Real results(testCase.f.size());
    for (auto& ri : results) ri.resize(testCase.rp.size());
    const Real a2 = testCase.r * testCase.r / (2.0 * testCase.sigma * testCase.sigma);
    Vec1Real   ab(testCase.rp.size());
    std::transform(testCase.rp.begin(),
                   testCase.rp.end(),
                   ab.begin(),
                   [&](const Real& rpi)
                   { return testCase.r * rpi / (testCase.sigma * testCase.sigma); });
    Vec1Real b2(testCase.rp.size());
    std::transform(testCase.rp.begin(),
                   testCase.rp.end(),
                   b2.begin(),
                   [&](const Real& rpi)
                   { return rpi * rpi / (2.0 * testCase.sigma * testCase.sigma); });
    math::expModSphBessel(a2, ab, b2, results);
    REQUIRE_CLOSE_COLLECTIONS_2D(results, testCase.f, testCase.tolerance, "expModSphBessel");
}

BOOST_DATA_TEST_CASE(ComputeMean,
                     bdata::make(containerMean.testCases),
                     testCase)
{
    Real computed = math::stats::calculateMean( testCase.data );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
}

BOOST_DATA_TEST_CASE(ComputeVariance,
                     bdata::make(containerVariance.testCases),
                     testCase)
{
    Real mean     = math::stats::calculateMean( testCase.data );
    Real computed = math::stats::calculateVariance( testCase.data, mean );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
    computed = math::stats::calculateVariance( testCase.data );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
}

BOOST_DATA_TEST_CASE(ComputeRMSE,
                     bdata::make(containerRMSE.testCases),
                     testCase)
{
    Real computed = math::stats::calculateRMSE( testCase.dataA, testCase.dataB );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
}

BOOST_DATA_TEST_CASE(ComputeMAE,
                     bdata::make(containerMAE.testCases),
                     testCase)
{
    Real computed = math::stats::calculateMAE( testCase.dataA, testCase.dataB );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
}

BOOST_DATA_TEST_CASE(PearsonCorr,
                     bdata::make(containerPearson.testCases),
                     testCase)
{
    Real computed = math::stats::pearsonCorr( testCase.dataA, testCase.dataB );
    BOOST_REQUIRE_CLOSE(computed, testCase.result, testCase.tolerance );
}


BOOST_AUTO_TEST_SUITE_END()
