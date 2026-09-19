#ifndef TEST_OPTIMIZE_HPP
#define TEST_OPTIMIZE_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "types.hpp"

#include <cmath>
#include <limits>
#include <type_traits>

namespace vaspml
{

Real const defaultTolerance = 100.0 * std::numeric_limits<Real>::epsilon();


template<typename T>
class Function
{
    static_assert(std::is_floating_point<T>::value, "Wrong argument type in CubicSpline\n");
    public:
    Function( void );
    T operator()( const T& x );
    T operator()( const std::vector<T>& x );
    void df( const std::vector<T>& x, std::vector<T>& grad );
    private:
    std::vector<T> params;
};


template<typename T>
Function<T>::Function( void )
{
    params.push_back( 1 );
    params.push_back( -1 );
    params.push_back( -10 );
    params.push_back( 1 );
    params.push_back( 5 );
}


template<typename T>
T Function<T>::operator()( const T& x )
{
    T x2 = x*x;
    T x3 = x2*x;
    T x4 = x2*x2;
    return params[0] + params[1] * x + params[2] * x2 + params[3] * x3 + params[4] * x4;
}

template<typename T>
T Function<T>::operator()( const std::vector<T>& x )
{
    if ( x[0]*x[0] + x[1]*x[1] < 5 )
    {
        return x[0] * x[0] + x[1] * x[1] - 3.0;
    }
    else
    {
        return ( x[0] * x[0] + x[1] * x[1] ) * std::cos( x[0] ) * std::cos( x[0] )
                                             * std::sin( x[1] ) * std::sin( x[1] );
    }
}


template<typename T>
void Function<T>::df( const std::vector<T>& x, std::vector<T>& grad )
{
    if ( x[0]*x[0] + x[1]*x[1] < 5  )
    {
        grad[0] = (T) 2 * x[0];
        grad[1] = (T) 2 * x[1];
    }
    else
    {
        grad[0] = (T)2 * x[0] * std::cos(x[0])*std::cos( x[0] ) * std::sin(x[1])*std::sin(x[1])
                    + (x[0]*x[0] + x[1]*x[1] ) * (-(T)2 * std::sin(x[0]) * std::cos(x[0])) * std::sin(x[1]) * std::sin (x[1]);
        grad[1] = (T)2 * x[1] * std::cos(x[0])*std::cos( x[0] ) * std::sin(x[1])*std::sin(x[1])
                    + (x[0]*x[0] + x[1]*x[1]) * std::cos(x[0]*std::cos(x[0])) * ((T)2 * std::sin(x[1])*std::cos(x[1]));
    }
}
    

struct TestCaseOptimize_BFGS : public TestCase
{
    Real     tolerance;
    Vec1Real xPosFinal;
    Vec1Real startPos;
    Real funcFinal;
    Real gradTol;
    Size nIterFinal;
    TestCaseOptimize_BFGS(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
inline void TestCaseContainer<TestCaseOptimize_BFGS>::setup()
{
    TestCaseOptimize_BFGS* tc = nullptr;
    testCases.push_back(TestCaseOptimize_BFGS("TestCaseOptimize_BFGS"));
    tc = &(testCases.back());

    tc-> xPosFinal.resize(2);
    tc->startPos.resize(2);
    tc->xPosFinal[0] = 0.0;
    tc->xPosFinal[1] = 0.0;
    
    tc->startPos[0] = 2.4;
    tc->startPos[1] = 2.1;

    tc->gradTol = 1e-14;
    tc->nIterFinal = 7;
    tc->funcFinal = -3.0;
}

struct TestCaseOptimize_Bracket : public TestCase
{   
    Real     tolerance;
    Vec1Real xBracket;
    Vec1Real functionBracket;
    Real x00;
    Real x01;
    TestCaseOptimize_Bracket(String name) : TestCase(name), tolerance(defaultTolerance) {}
};  

template<>
inline void TestCaseContainer<TestCaseOptimize_Bracket>::setup()
{
    TestCaseOptimize_Bracket* tc = nullptr;
    testCases.push_back(TestCaseOptimize_Bracket("TestCaseOptimize_Bracket"));
    tc = &(testCases.back());

    tc->xBracket.resize(3);
    tc->functionBracket.resize(3);
    tc->xBracket[0]  =  -5.000000000000000;
    tc->xBracket[1]  =   1.000000000000000;
    tc->xBracket[2]  =  10.708204000000000;
    tc->functionBracket[0]  =  2756.000000000000000;
    tc->functionBracket[1]  =  -4.000000000000000;
    tc->functionBracket[2]  =  65812.535304113727761;
    tc->x00 =  1.0;
    tc->x01 = -5.0;
    tc->tolerance = 2.0E-11;
}

}//namespace vaspml
#endif
