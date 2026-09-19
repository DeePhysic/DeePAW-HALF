#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Optimize

#include "test_Optimize.hpp"

#include "boost_helpers.hpp"

#include "BFGS.hpp"
#include "Bracket.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseOptimize_BFGS>    container_BFGS;
TestCaseContainer<TestCaseOptimize_Bracket> container_Bracket;

BOOST_AUTO_TEST_SUITE(UnitTests)
BOOST_DATA_TEST_CASE( TestCaseOptimize_BFGS, 
                      bdata::make(container_BFGS.testCases), testCase)
{
    Vec1Real pos = testCase.startPos;
    BFGS<Real> bfgs;
    Vec1Real funcValue(1);
    Vec1Size nIter(1);
    Function<Real> func;
    bfgs.minimize( pos, testCase.gradTol, nIter[0], funcValue[0], func );
    REQUIRE_CLOSE_COLLECTIONS( pos, testCase.xPosFinal, testCase.tolerance, "BFGS::minimum_pos");
    Vec1Real desiredValue(1);
    desiredValue[0] = testCase.funcFinal;
    REQUIRE_CLOSE_COLLECTIONS( funcValue, desiredValue, testCase.tolerance, "BFGS::function_value" );
    Vec1Size nIterDesired(1);
    nIterDesired[0] = testCase.nIterFinal;
    REQUIRE_EQUAL_COLLECTIONS( nIter, nIterDesired, "BFGS::number_iterations" );
}

BOOST_DATA_TEST_CASE( TestCaseOptimize_Bracket, 
                      bdata::make(container_Bracket.testCases), testCase)
{
    BracketMinimum<double> bracket;
    Function<Real> func;
    bracket.bracket( testCase.x00, testCase.x01, func );
    Vec1Real resultX;
    resultX.push_back( bracket.get_ax() );
    resultX.push_back( bracket.get_bx() );
    resultX.push_back( bracket.get_cx() );
    Vec1Real resultFunc;
    resultFunc.push_back( bracket.get_fa() );
    resultFunc.push_back( bracket.get_fb() );
    resultFunc.push_back( bracket.get_fc() );

    REQUIRE_CLOSE_COLLECTIONS( resultX, testCase.xBracket, testCase.tolerance, "BracketS::x_values" );
    REQUIRE_CLOSE_COLLECTIONS( resultFunc, testCase.functionBracket, testCase.tolerance, "BracketS::function_values" );
}

BOOST_AUTO_TEST_SUITE_END()
