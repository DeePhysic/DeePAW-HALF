#ifndef TEST_SPILLINGFACTOR_HPP
#define TEST_SPILLINGFACTOR_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "Record.hpp"
#include "io.hpp"
#include "types.hpp"

#include <limits>

namespace vaspml
{

Real const defaultTolerance = 1000000000.0 * std::numeric_limits<Real>::epsilon();

struct TestCase_SpillingFactor : public TestCase
{
    Real     tolerance;
    String   structure;
    Record   ff;
    Vec1Real sfAtom;
    Real     min;
    Real     max;
    Real     mean;
    Real     pos0;
    TestCase_SpillingFactor(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCase_SpillingFactor>::setup()
{
    String              structureName = "";
    TestCase_SpillingFactor* tc = nullptr;

    structureName = "MAPbI3.cubic";
    testCases.push_back(TestCase_SpillingFactor("TestCase_SpillingFactor_" + structureName));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/6.5.1/refit/ML_FF.MAPbI3");
    tc->sfAtom.push_back(0.43166408096669784    ); 
    tc->sfAtom.push_back(4.2555204978953043E-002);
    tc->sfAtom.push_back(3.5816833735907494E-002);
    tc->sfAtom.push_back(0.47656826643026307    ); 
    tc->sfAtom.push_back(4.3745026236619022E-004);
    tc->sfAtom.push_back(7.9243402148931885E-005);
    tc->sfAtom.push_back(2.8193500279585226E-003);
    tc->sfAtom.push_back(4.8581335719602992E-003);
    tc->sfAtom.push_back(4.8581386228089274E-003);
    tc->sfAtom.push_back(2.7844214880654583E-003);
    tc->sfAtom.push_back(2.7844205365976649E-003);
    tc->sfAtom.push_back(1.4378487493627023E-004);
    tc->tolerance = defaultTolerance;
    tc->pos0      = 0.6400224700078113;
    tc->min       = 7.92434021e-05;
    tc->max       = 4.76568266e-01;
    tc->mean      = 8.37807774e-02;

    return;
}

} //namespace vaspml

#endif
