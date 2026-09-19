#ifndef TEST_BASISFUNCTIONSRADIALSPLINE_HPP
#define TEST_BASISFUNCTIONSRADIALSPLINE_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{

// test of accuracy relative 1e-7
Real const defaultTolerance = 10000000000.0 * std::numeric_limits<Real>::epsilon();

struct TestCase_BasisFunctionsRadialSpline : public TestCase
{

    Real     tolerance;
    Vec1Real values;
    Vec1Real derivative;

    TestCase_BasisFunctionsRadialSpline(String name) : TestCase(name), tolerance(defaultTolerance)
    {}
};

template<>
void TestCaseContainer<TestCase_BasisFunctionsRadialSpline>::setup()
{
    TestCase_BasisFunctionsRadialSpline* tc = nullptr;

    testCases.push_back(TestCase_BasisFunctionsRadialSpline("BasisFunctionsRadialSpline.update()"));
    tc = &(testCases.back());

    tc->values.resize(22);

    tc->values[0] = 0.36295320603219249911;
    tc->values[1] = 0.5953772563943472429;
    tc->values[2] = 0.64028501345248456378;
    tc->values[3] = 0.53284514088583068148;
    tc->values[4] = 0.35934411018769513158;
    tc->values[5] = 0.1985955707494624356;
    tc->values[6] = 0.08894272563471666837;
    tc->values[7] = 0.030824631099712073573;

    tc->values[8] = 0.07536847958467078612;
    tc->values[9] = 0.17293242109392251149;
    tc->values[10] = 0.24012861106443025849;
    tc->values[11] = 0.247441757065722856;
    tc->values[12] = 0.20336143747240431612;
    tc->values[13] = 0.13776474607734762179;
    tc->values[14] = 0.078183685399780794567;

    tc->values[15] = 0.013651800692651289;
    tc->values[16] = 0.039280395147663197;
    tc->values[17] = 0.065199470067950349;
    tc->values[18] = 0.077860321721855641;
    tc->values[19] = 0.072751154583277161;
    tc->values[20] = 0.055436636313064197;
    tc->values[21] = 0.035257578581798045;

    tc->derivative.resize(22);

    tc->derivative[0] = -0.060159595646645206;
    tc->derivative[1] = -0.22023203762450785;
    tc->derivative[2] = -0.46741795959427745;
    tc->derivative[3] = -0.68358893009813415;
    tc->derivative[4] = -0.75444861869686597;
    tc->derivative[5] = -0.66021786344476374;
    tc->derivative[6] = -0.47095317008294091;
    tc->derivative[7] = -0.27797543001640518;

    tc->derivative[8] = 0.13711418766420119;
    tc->derivative[9] = 0.28665009020259735;
    tc->derivative[10] = 0.33805840136155729;
    tc->derivative[11] = 0.26227935326912544;
    tc->derivative[12] = 0.12167150342811399;
    tc->derivative[13] = -0.00023269310019207094;
    tc->derivative[14] = -0.060923966079636044;

    tc->derivative[15] = 0.051946088933082439;
    tc->derivative[16] = 0.14380442693779893;
    tc->derivative[17] = 0.22531135611118677;
    tc->derivative[18] = 0.24796339506618859;
    tc->derivative[19] = 0.20673210273102602;
    tc->derivative[20] = 0.13402424949439823;
    tc->derivative[21] = 0.066997074791968980;

    return;
}

} //namespace vaspml

#endif
