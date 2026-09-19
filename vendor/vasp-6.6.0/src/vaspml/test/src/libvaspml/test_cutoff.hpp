#ifndef TEST_CUTOFF_HPP
#define TEST_CUTOFF_HPP

#include "TestCaseContainer.hpp"
#include "TestCaseFunction1D.hpp"

#include "types.hpp"

namespace vaspml
{

struct TestCase_cutoff : public TestCaseFunction1D
{
    /// Cutoff function type as string.
    String cutoffType;

    TestCase_cutoff(String name) : TestCaseFunction1D(name) {};
};

template<>
inline void TestCaseContainer<TestCase_cutoff>::setup()
{
    TestCase_cutoff* tc = nullptr;

    Real   rcut = 6.0;
    String type = "BP";
    testCases.push_back(TestCase_cutoff("type = " + type + ", rcut = " + std::to_string(rcut)));
    tc = &(testCases.back());
    tc->parameters["rcut"] = rcut;
    tc->cutoffType = type;
    // regular grid x-values: xmin = -0.5, xmax = 6.5, n = 11
    // function parameters =  {'rcut': 6.0}
    // clang-format off
    //tc->x.push_back(-5.0000000000000000E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 1.9999999999999996E-01); tc->fx.push_back( 9.9726094768413664E-01); tc->dfx.push_back(-2.7365487691057549E-02);
    tc->x.push_back( 8.9999999999999991E-01); tc->fx.push_back( 9.4550326209418389E-01); tc->dfx.push_back(-1.1885443489844326E-01);
    tc->x.push_back( 1.5999999999999996E+00); tc->fx.push_back( 8.3456530317942912E-01); tc->dfx.push_back(-1.9455486035608752E-01);
    tc->x.push_back( 2.2999999999999998E+00); tc->fx.push_back( 6.7918397477265025E-01); tc->dfx.push_back(-2.4441078411823625E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 5.0000000000000000E-01); tc->dfx.push_back(-2.6179938779914941E-01);
    tc->x.push_back( 3.6999999999999993E+00); tc->fx.push_back( 3.2081602522734998E-01); tc->dfx.push_back(-2.4441078411823627E-01);
    tc->x.push_back( 4.3999999999999995E+00); tc->fx.push_back( 1.6543469682057105E-01); tc->dfx.push_back(-1.9455486035608760E-01);
    tc->x.push_back( 5.0999999999999996E+00); tc->fx.push_back( 5.4496737905816106E-02); tc->dfx.push_back(-1.1885443489844327E-01);
    tc->x.push_back( 5.7999999999999998E+00); tc->fx.push_back( 2.7390523158633551E-03); tc->dfx.push_back(-2.7365487691057625E-02);
    tc->x.push_back( 6.5000000000000000E+00); tc->fx.push_back( 0.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    // clang-format on

    rcut = 6.0;
    type = "MIWA";
    testCases.push_back(TestCase_cutoff("type = " + type + ", rcut = " + std::to_string(rcut)));
    tc = &(testCases.back());
    tc->parameters["rcut"] = rcut;
    tc->cutoffType = type;
    // regular grid x-values: xmin = -0.5, xmax = 6.5, n = 11
    // function parameters =  {'rcut': 6.0}
    // clang-format off
    //tc->x.push_back(-5.0000000000000000E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 1.9999999999999996E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 8.9999999999999991E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 1.5999999999999996E+00); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 2.2999999999999998E+00); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 3.6999999999999993E+00); tc->fx.push_back( 8.6207407407407421E-01); tc->dfx.push_back(-3.5777777777777758E-01);
    tc->x.push_back( 4.3999999999999995E+00); tc->fx.push_back( 5.4992592592592604E-01); tc->dfx.push_back(-4.9777777777777782E-01);
    tc->x.push_back( 5.0999999999999996E+00); tc->fx.push_back( 2.1600000000000008E-01); tc->dfx.push_back(-4.2000000000000010E-01);
    tc->x.push_back( 5.7999999999999998E+00); tc->fx.push_back( 1.2740740740740719E-02); tc->dfx.push_back(-1.2444444444444443E-01);
    tc->x.push_back( 6.5000000000000000E+00); tc->fx.push_back( 0.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    // clang-format on

    rcut = 6.0;
    type = "POLY2";
    testCases.push_back(TestCase_cutoff("type = " + type + ", rcut = " + std::to_string(rcut)));
    tc = &(testCases.back());
    tc->parameters["rcut"] = rcut;
    tc->cutoffType = type;
    // regular grid x-values: xmin = -0.5, xmax = 6.5, n = 11
    // function parameters =  {'rcut': 6.0}
    // clang-format off
    //tc->x.push_back(-5.0000000000000000E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 1.9999999999999996E-01); tc->fx.push_back( 9.9964790123456793E-01); tc->dfx.push_back(-5.1913580246913558E-03);
    tc->x.push_back( 8.9999999999999991E-01); tc->fx.push_back( 9.7338812500000005E-01); tc->dfx.push_back(-8.1281250000000013E-02);
    tc->x.push_back( 1.5999999999999996E+00); tc->fx.push_back( 8.7813135802469144E-01); tc->dfx.push_back(-1.9120987654320987E-01);
    tc->x.push_back( 2.2999999999999998E+00); tc->fx.push_back( 7.1093986882716043E-01); tc->dfx.push_back(-2.7939853395061737E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 5.0000000000000000E-01); tc->dfx.push_back(-3.1250000000000000E-01);
    tc->x.push_back( 3.6999999999999993E+00); tc->fx.push_back( 2.8906013117283968E-01); tc->dfx.push_back(-2.7939853395061731E-01);
    tc->x.push_back( 4.3999999999999995E+00); tc->fx.push_back( 1.2186864197530889E-01); tc->dfx.push_back(-1.9120987654320976E-01);
    tc->x.push_back( 5.0999999999999996E+00); tc->fx.push_back( 2.6611875000000729E-02); tc->dfx.push_back(-8.1281249999999527E-02);
    tc->x.push_back( 5.7999999999999998E+00); tc->fx.push_back( 3.5209876543151886E-04); tc->dfx.push_back(-5.1913580246917652E-03);
    tc->x.push_back( 6.5000000000000000E+00); tc->fx.push_back( 0.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    // clang-format on

    rcut = 6.0;
    type = "BUMP";
    testCases.push_back(TestCase_cutoff("type = " + type + ", rcut = " + std::to_string(rcut)));
    tc = &(testCases.back());
    tc->parameters["rcut"] = rcut;
    tc->cutoffType = type;
    // regular grid x-values: xmin = -0.5, xmax = 6.5, n = 11
    // function parameters =  {'rcut': 6.0}
    // clang-format off
    //tc->x.push_back(-5.0000000000000000E-01); tc->fx.push_back( 1.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    tc->x.push_back( 1.9999999999999996E-01); tc->fx.push_back( 9.9888827137637848E-01); tc->dfx.push_back(-1.1123463646280323E-02);
    tc->x.push_back( 8.9999999999999991E-01); tc->fx.push_back( 9.7724498818209626E-01); tc->dfx.push_back(-5.1137550810478542E-02);
    tc->x.push_back( 1.5999999999999996E+00); tc->fx.push_back( 9.2630194411928546E-01); tc->dfx.push_back(-9.5427273838437396E-02);
    tc->x.push_back( 2.2999999999999998E+00); tc->fx.push_back( 8.4176315412801961E-01); tc->dfx.push_back(-1.4780550558904373E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 7.1653131057378927E-01); tc->dfx.push_back(-2.1230557350334497E-01);
    tc->x.push_back( 3.6999999999999993E+00); tc->fx.push_back( 5.4138415085500935E-01); tc->dfx.push_back(-2.8976145750283022E-01);
    tc->x.push_back( 4.3999999999999995E+00); tc->fx.push_back( 3.1240291237977874E-01); tc->dfx.push_back(-3.5743214133688617E-01);
    tc->x.push_back( 5.0999999999999996E+00); tc->fx.push_back( 7.4006407358280485E-02); tc->dfx.push_back(-2.7229584721819494E-01);
    tc->x.push_back( 5.7999999999999998E+00); tc->fx.push_back( 6.4485697550942261E-07); tc->dfx.push_back(-4.8350379376029631E-05);
    tc->x.push_back( 6.5000000000000000E+00); tc->fx.push_back( 0.0000000000000000E+00); tc->dfx.push_back( 0.0000000000000000E+00);
    // clang-format on
}

} //namespace vaspml

#endif
