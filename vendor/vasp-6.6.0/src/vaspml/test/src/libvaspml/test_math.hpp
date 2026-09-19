#ifndef TEST_MATH_HPP
#define TEST_MATH_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "TestCaseFunction1D.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{

struct TestCase_math_clebschGordan : public TestCase
{
    struct SingleClebschGordan
    {
        Int  j1;
        Int  m1;
        Int  j2;
        Int  m2;
        Int  J;
        Int  M;
        Real coeff;
    };

    Int                              jmax;
    Real                             tolerance = 10 * std::numeric_limits<Real>::epsilon();
    std::vector<SingleClebschGordan> cgs;

    TestCase_math_clebschGordan(String name) : TestCase(name) {};
};

template<>
inline void TestCaseContainer<TestCase_math_clebschGordan>::setup()
{
    TestCase_math_clebschGordan* tc = nullptr;
    using SingleClebschGordan = TestCase_math_clebschGordan::SingleClebschGordan;

    Int jmax = 3;
    testCases.push_back(TestCase_math_clebschGordan("j1max = j2max = " + std::to_string(jmax)));
    tc = &(testCases.back());
    tc->jmax = jmax;
    //// j1-index from 0 to jmax, m1-index from -j1 to j1.
    //Int n = (jmax + 1) * (2 * jmax + 1);
    //// Same for j2, m2.
    //n *= n;
    //// J-index from 0 to jmax * jmax, M-index from -J to J.
    //n *= (2 * jmax + 1) * (4 * jmax + 3);
    for (Int j1 = 0; j1 <= jmax; ++j1)
    {
        for (Int m1 = -j1; m1 <= j1; ++m1)
        {
            for (Int j2 = 0; j2 <= jmax; ++j2)
            {
                for (Int m2 = -j2; m2 <= j2; ++m2)
                {
                    for (Int J = 0; J <= j1 + j2; ++J)
                    {
                        for (Int M = -J; M <= J; ++M)
                        {
                            tc->cgs.push_back(SingleClebschGordan({j1, m1, j2, m2, J, M, 0.0}));
                        }
                    }
                }
            }
        }
    }

    // clang-format off
    tc->cgs[    0].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   0   0   0   0 1
    tc->cgs[    2].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   1  -1   1  -1 1
    tc->cgs[    7].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   1   0   1   0 1
    tc->cgs[   12].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   1   1   1   1 1
    tc->cgs[   17].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   2  -2   2  -2 1
    tc->cgs[   27].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   2  -1   2  -1 1
    tc->cgs[   37].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   2   0   2   0 1
    tc->cgs[   47].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   2   1   2   1 1
    tc->cgs[   57].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   2   2   2   2 1
    tc->cgs[   67].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3  -3   3  -3 1
    tc->cgs[   84].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3  -2   3  -2 1
    tc->cgs[  101].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3  -1   3  -1 1
    tc->cgs[  118].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3   0   3   0 1
    tc->cgs[  135].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3   1   3   1 1
    tc->cgs[  152].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3   2   3   2 1
    tc->cgs[  169].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   0   0   3   3   3   3 1
    tc->cgs[  171].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1  -1   0   0   1  -1 1
    tc->cgs[  178].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1  -1   1  -1   2  -2 1
    tc->cgs[  184].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   1   0   1  -1 -sqrt(2)/2
    tc->cgs[  188].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   1   0   2  -1 sqrt(2)/2
    tc->cgs[  192].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   1   1   0   0 sqrt(3)/3
    tc->cgs[  194].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   1   1   1   0 -sqrt(2)/2
    tc->cgs[  198].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   1  -1   1   1   2   0 sqrt(6)/6
    tc->cgs[  210].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1  -1   2  -2   3  -3 1
    tc->cgs[  221].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2  -1   2  -2 -sqrt(3)/3
    tc->cgs[  227].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2  -1   3  -2 sqrt(6)/3
    tc->cgs[  234].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   0   1  -1 sqrt(10)/10
    tc->cgs[  238].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   0   2  -1 -sqrt(2)/2
    tc->cgs[  244].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   0   3  -1 sqrt(10)/5
    tc->cgs[  251].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   1   1   0 sqrt(30)/10
    tc->cgs[  255].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   1   2   0 -sqrt(2)/2
    tc->cgs[  261].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   1   3   0 sqrt(5)/5
    tc->cgs[  268].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   2   1   1 sqrt(15)/5
    tc->cgs[  272].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   2   2   1 -sqrt(3)/3
    tc->cgs[  278].coeff =   2.5819888974716100E-1; // j1 m1 j2 m2 J M coeff:   1  -1   2   2   3   1 sqrt(15)/15
    tc->cgs[  297].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1  -1   3  -3   4  -4 1
    tc->cgs[  315].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3  -2   3  -3 -1/2
    tc->cgs[  323].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3  -2   4  -3 sqrt(3)/2
    tc->cgs[  335].coeff =   2.1821789023599200E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3  -1   2  -2 sqrt(21)/21
    tc->cgs[  341].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3  -1   3  -2 -sqrt(15)/6
    tc->cgs[  349].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3  -1   4  -2 sqrt(105)/14
    tc->cgs[  361].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   0   2  -1 sqrt(7)/7
    tc->cgs[  367].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   0   3  -1 -sqrt(2)/2
    tc->cgs[  375].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   0   4  -1 sqrt(70)/14
    tc->cgs[  387].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   1   2   0 sqrt(14)/7
    tc->cgs[  393].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   1   3   0 -sqrt(2)/2
    tc->cgs[  401].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   1   4   0 sqrt(42)/14
    tc->cgs[  413].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   2   2   1 sqrt(210)/21
    tc->cgs[  419].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   2   3   1 -sqrt(15)/6
    tc->cgs[  427].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   2   4   1 sqrt(21)/14
    tc->cgs[  439].coeff =   8.4515425472851700E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   3   2   2 sqrt(35)/7
    tc->cgs[  445].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   3   3   2 -1/2
    tc->cgs[  453].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   1  -1   3   3   4   2 sqrt(7)/14
    tc->cgs[  458].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1   0   0   0   1   0 1
    tc->cgs[  461].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   0   1  -1   1  -1 sqrt(2)/2
    tc->cgs[  465].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   0   1  -1   2  -1 sqrt(2)/2
    tc->cgs[  469].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   0   1   0   0   0 -sqrt(3)/3
    tc->cgs[  475].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   1   0   1   0   2   0 sqrt(6)/3
    tc->cgs[  481].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   0   1   1   1   1 -sqrt(2)/2
    tc->cgs[  485].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   0   1   1   2   1 sqrt(2)/2
    tc->cgs[  491].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2  -2   2  -2 sqrt(6)/3
    tc->cgs[  497].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2  -2   3  -2 sqrt(3)/3
    tc->cgs[  504].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2  -1   1  -1 -sqrt(30)/10
    tc->cgs[  508].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   1   0   2  -1   2  -1 sqrt(6)/6
    tc->cgs[  514].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   1   0   2  -1   3  -1 2*sqrt(30)/15
    tc->cgs[  521].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   0   1   0 -sqrt(10)/5
    tc->cgs[  531].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   0   3   0 sqrt(15)/5
    tc->cgs[  538].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   1   1   1 -sqrt(30)/10
    tc->cgs[  542].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   1   2   1 -sqrt(6)/6
    tc->cgs[  548].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   1   3   1 2*sqrt(30)/15
    tc->cgs[  559].coeff =  -8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   2   2   2 -sqrt(6)/3
    tc->cgs[  565].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   0   2   2   3   2 sqrt(3)/3
    tc->cgs[  576].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -3   3  -3 sqrt(3)/2
    tc->cgs[  584].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -3   4  -3 1/2
    tc->cgs[  596].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -2   2  -2 -sqrt(105)/21
    tc->cgs[  602].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -2   3  -2 sqrt(3)/3
    tc->cgs[  610].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -2   4  -2 sqrt(21)/7
    tc->cgs[  622].coeff =  -6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -1   2  -1 -2*sqrt(42)/21
    tc->cgs[  628].coeff =   2.8867513459481300E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -1   3  -1 sqrt(3)/6
    tc->cgs[  636].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   1   0   3  -1   4  -1 sqrt(105)/14
    tc->cgs[  648].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   0   2   0 -sqrt(21)/7
    tc->cgs[  662].coeff =   7.5592894601845500E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   0   4   0 2*sqrt(7)/7
    tc->cgs[  674].coeff =  -6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   1   2   1 -2*sqrt(42)/21
    tc->cgs[  680].coeff =  -2.8867513459481300E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   1   3   1 -sqrt(3)/6
    tc->cgs[  688].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   1   4   1 sqrt(105)/14
    tc->cgs[  700].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   2   2   2 -sqrt(105)/21
    tc->cgs[  706].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   2   3   2 -sqrt(3)/3
    tc->cgs[  714].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   2   4   2 sqrt(21)/7
    tc->cgs[  732].coeff =  -8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   3   3   3 -sqrt(3)/2
    tc->cgs[  740].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1   0   3   3   4   3 1/2
    tc->cgs[  745].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1   1   0   0   1   1 1
    tc->cgs[  746].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   1   1  -1   0   0 sqrt(3)/3
    tc->cgs[  748].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   1  -1   1   0 sqrt(2)/2
    tc->cgs[  752].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   1   1   1  -1   2   0 sqrt(6)/6
    tc->cgs[  758].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   1   0   1   1 sqrt(2)/2
    tc->cgs[  762].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   1   0   2   1 sqrt(2)/2
    tc->cgs[  772].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1   1   1   1   2   2 1
    tc->cgs[  774].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -2   1  -1 sqrt(15)/5
    tc->cgs[  778].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -2   2  -1 sqrt(3)/3
    tc->cgs[  784].coeff =   2.5819888974716100E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -2   3  -1 sqrt(15)/15
    tc->cgs[  791].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -1   1   0 sqrt(30)/10
    tc->cgs[  795].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -1   2   0 sqrt(2)/2
    tc->cgs[  801].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   1   1   2  -1   3   0 sqrt(5)/5
    tc->cgs[  808].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   1   1   2   0   1   1 sqrt(10)/10
    tc->cgs[  812].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   2   0   2   1 sqrt(2)/2
    tc->cgs[  818].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   1   1   2   0   3   1 sqrt(10)/5
    tc->cgs[  829].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   1   1   2   1   2   2 sqrt(3)/3
    tc->cgs[  835].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   1   1   2   1   3   2 sqrt(6)/3
    tc->cgs[  852].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1   1   2   2   3   3 1
    tc->cgs[  857].coeff =   8.4515425472851700E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -3   2  -2 sqrt(35)/7
    tc->cgs[  863].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -3   3  -2 1/2
    tc->cgs[  871].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -3   4  -2 sqrt(7)/14
    tc->cgs[  883].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -2   2  -1 sqrt(210)/21
    tc->cgs[  889].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -2   3  -1 sqrt(15)/6
    tc->cgs[  897].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -2   4  -1 sqrt(21)/14
    tc->cgs[  909].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -1   2   0 sqrt(14)/7
    tc->cgs[  915].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -1   3   0 sqrt(2)/2
    tc->cgs[  923].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   1   1   3  -1   4   0 sqrt(42)/14
    tc->cgs[  935].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   0   2   1 sqrt(7)/7
    tc->cgs[  941].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   0   3   1 sqrt(2)/2
    tc->cgs[  949].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   0   4   1 sqrt(70)/14
    tc->cgs[  961].coeff =   2.1821789023599200E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   1   2   2 sqrt(21)/21
    tc->cgs[  967].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   1   3   2 sqrt(15)/6
    tc->cgs[  975].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   1   4   2 sqrt(105)/14
    tc->cgs[  993].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   2   3   3 1/2
    tc->cgs[ 1001].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   1   1   3   2   4   3 sqrt(3)/2
    tc->cgs[ 1027].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   1   1   3   3   4   4 1
    tc->cgs[ 1032].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2  -2   0   0   2  -2 1
    tc->cgs[ 1046].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2  -2   1  -1   3  -3 1
    tc->cgs[ 1057].coeff =  -8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   1   0   2  -2 -sqrt(6)/3
    tc->cgs[ 1063].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   1   0   3  -2 sqrt(3)/3
    tc->cgs[ 1070].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   1   1   1  -1 sqrt(15)/5
    tc->cgs[ 1074].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   1   1   2  -1 -sqrt(3)/3
    tc->cgs[ 1080].coeff =   2.5819888974716100E-1; // j1 m1 j2 m2 J M coeff:   2  -2   1   1   3  -1 sqrt(15)/15
    tc->cgs[ 1101].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2  -2   2  -2   4  -4 1
    tc->cgs[ 1119].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2  -1   3  -3 -sqrt(2)/2
    tc->cgs[ 1127].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2  -1   4  -3 sqrt(2)/2
    tc->cgs[ 1139].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   0   2  -2 sqrt(14)/7
    tc->cgs[ 1145].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   0   3  -2 -sqrt(2)/2
    tc->cgs[ 1153].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   0   4  -2 sqrt(42)/14
    tc->cgs[ 1161].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   1   1  -1 -sqrt(5)/5
    tc->cgs[ 1165].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   1   2  -1 sqrt(21)/7
    tc->cgs[ 1171].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   1   3  -1 -sqrt(30)/10
    tc->cgs[ 1179].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   1   4  -1 sqrt(14)/14
    tc->cgs[ 1185].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   2   0   0 sqrt(5)/5
    tc->cgs[ 1187].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   2   1   0 -sqrt(10)/5
    tc->cgs[ 1191].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   2   2   0 sqrt(14)/7
    tc->cgs[ 1197].coeff =  -3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   2   3   0 -sqrt(10)/10
    tc->cgs[ 1205].coeff =   1.1952286093343900E-1; // j1 m1 j2 m2 J M coeff:   2  -2   2   2   4   0 sqrt(70)/70
    tc->cgs[ 1235].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2  -2   3  -3   5  -5 1
    tc->cgs[ 1262].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3  -2   4  -4 -sqrt(10)/5
    tc->cgs[ 1272].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3  -2   5  -4 sqrt(15)/5
    tc->cgs[ 1291].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3  -1   3  -3 sqrt(6)/6
    tc->cgs[ 1299].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3  -1   4  -3 -sqrt(2)/2
    tc->cgs[ 1309].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3  -1   5  -3 sqrt(3)/3
    tc->cgs[ 1322].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   0   2  -2 -sqrt(14)/14
    tc->cgs[ 1328].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   0   3  -2 sqrt(3)/3
    tc->cgs[ 1336].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   0   4  -2 -sqrt(21)/7
    tc->cgs[ 1346].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   0   5  -2 sqrt(6)/6
    tc->cgs[ 1355].coeff =   1.6903085094570300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   1   1  -1 sqrt(35)/35
    tc->cgs[ 1359].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   1   2  -1 -sqrt(42)/14
    tc->cgs[ 1365].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   1   3  -1 sqrt(10)/5
    tc->cgs[ 1373].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   1   4  -1 -sqrt(14)/7
    tc->cgs[ 1383].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   1   5  -1 sqrt(14)/14
    tc->cgs[ 1392].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   2   1   0 sqrt(7)/7
    tc->cgs[ 1396].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   2   2   0 -sqrt(70)/14
    tc->cgs[ 1402].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   2   3   0 sqrt(3)/3
    tc->cgs[ 1410].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   2   4   0 -sqrt(7)/7
    tc->cgs[ 1420].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   2   5   0 sqrt(42)/42
    tc->cgs[ 1429].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   3   1   1 sqrt(21)/7
    tc->cgs[ 1433].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   3   2   1 -sqrt(70)/14
    tc->cgs[ 1439].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   3   3   1 sqrt(6)/6
    tc->cgs[ 1447].coeff =  -2.0701966780270600E-1; // j1 m1 j2 m2 J M coeff:   2  -2   3   3   4   1 -sqrt(210)/70
    tc->cgs[ 1457].coeff =   6.9006555934235400E-2; // j1 m1 j2 m2 J M coeff:   2  -2   3   3   5   1 sqrt(210)/210
    tc->cgs[ 1467].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2  -1   0   0   2  -1 1
    tc->cgs[ 1475].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1  -1   2  -2 sqrt(3)/3
    tc->cgs[ 1481].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1  -1   3  -2 sqrt(6)/3
    tc->cgs[ 1488].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   0   1  -1 -sqrt(30)/10
    tc->cgs[ 1492].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   0   2  -1 -sqrt(6)/6
    tc->cgs[ 1498].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   0   3  -1 2*sqrt(30)/15
    tc->cgs[ 1505].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   1   1   0 sqrt(30)/10
    tc->cgs[ 1509].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   1   2   0 -sqrt(2)/2
    tc->cgs[ 1515].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   1   1   3   0 sqrt(5)/5
    tc->cgs[ 1528].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2  -2   3  -3 sqrt(2)/2
    tc->cgs[ 1536].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2  -2   4  -3 sqrt(2)/2
    tc->cgs[ 1548].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2  -1   2  -2 -sqrt(21)/7
    tc->cgs[ 1562].coeff =   7.5592894601845500E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2  -1   4  -2 2*sqrt(7)/7
    tc->cgs[ 1570].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   0   1  -1 sqrt(30)/10
    tc->cgs[ 1574].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   0   2  -1 -sqrt(14)/14
    tc->cgs[ 1580].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   0   3  -1 -sqrt(5)/5
    tc->cgs[ 1588].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   0   4  -1 sqrt(21)/7
    tc->cgs[ 1594].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   1   0   0 -sqrt(5)/5
    tc->cgs[ 1596].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   1   1   0 sqrt(10)/10
    tc->cgs[ 1600].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   1   2   0 sqrt(14)/14
    tc->cgs[ 1606].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   1   3   0 -sqrt(10)/5
    tc->cgs[ 1614].coeff =   4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   1   4   0 2*sqrt(70)/35
    tc->cgs[ 1622].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   2   1   1 -sqrt(5)/5
    tc->cgs[ 1626].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   2   2   1 sqrt(21)/7
    tc->cgs[ 1632].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   2   3   1 -sqrt(30)/10
    tc->cgs[ 1640].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2  -1   2   2   4   1 sqrt(14)/14
    tc->cgs[ 1660].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -3   4  -4 sqrt(15)/5
    tc->cgs[ 1670].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -3   5  -4 sqrt(10)/5
    tc->cgs[ 1689].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -2   3  -3 -sqrt(15)/6
    tc->cgs[ 1697].coeff =   2.2360679774997900E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -2   4  -3 sqrt(5)/10
    tc->cgs[ 1707].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -2   5  -3 2*sqrt(30)/15
    tc->cgs[ 1720].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -1   2  -2 sqrt(42)/14
    tc->cgs[ 1726].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -1   3  -2 -1/2
    tc->cgs[ 1734].coeff =  -1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -1   4  -2 -sqrt(7)/14
    tc->cgs[ 1744].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3  -1   5  -2 sqrt(2)/2
    tc->cgs[ 1753].coeff =  -2.9277002188456000E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   0   1  -1 -sqrt(105)/35
    tc->cgs[ 1757].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   0   2  -1 sqrt(14)/7
    tc->cgs[ 1763].coeff =  -1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   0   3  -1 -sqrt(30)/30
    tc->cgs[ 1771].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   0   4  -1 -sqrt(42)/14
    tc->cgs[ 1781].coeff =   6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   0   5  -1 2*sqrt(42)/21
    tc->cgs[ 1790].coeff =  -4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   1   1   0 -2*sqrt(70)/35
    tc->cgs[ 1794].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   1   2   0 sqrt(7)/7
    tc->cgs[ 1800].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   1   3   0 sqrt(30)/30
    tc->cgs[ 1808].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   1   4   0 -sqrt(70)/14
    tc->cgs[ 1818].coeff =   4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   1   5   0 sqrt(105)/21
    tc->cgs[ 1827].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   2   1   1 -sqrt(14)/7
    tc->cgs[ 1837].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   2   3   1 1/2
    tc->cgs[ 1845].coeff =  -5.9160797830996200E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   2   4   1 -sqrt(35)/10
    tc->cgs[ 1855].coeff =   3.3806170189140700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   2   5   1 2*sqrt(35)/35
    tc->cgs[ 1868].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   3   2   2 -sqrt(70)/14
    tc->cgs[ 1874].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   3   3   2 sqrt(15)/6
    tc->cgs[ 1882].coeff =  -4.3915503282684000E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   3   4   2 -3*sqrt(105)/70
    tc->cgs[ 1892].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2  -1   3   3   5   2 sqrt(30)/30
    tc->cgs[ 1902].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   0   0   0   2   0 1
    tc->cgs[ 1906].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2   0   1  -1   1  -1 sqrt(10)/10
    tc->cgs[ 1910].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   0   1  -1   2  -1 sqrt(2)/2
    tc->cgs[ 1916].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   0   1  -1   3  -1 sqrt(10)/5
    tc->cgs[ 1923].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   0   1   0   1   0 -sqrt(10)/5
    tc->cgs[ 1933].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2   0   1   0   3   0 sqrt(15)/5
    tc->cgs[ 1940].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2   0   1   1   1   1 sqrt(10)/10
    tc->cgs[ 1944].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   0   1   1   2   1 -sqrt(2)/2
    tc->cgs[ 1950].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   0   1   1   3   1 sqrt(10)/5
    tc->cgs[ 1957].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -2   2  -2 sqrt(14)/7
    tc->cgs[ 1963].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -2   3  -2 sqrt(2)/2
    tc->cgs[ 1971].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -2   4  -2 sqrt(42)/14
    tc->cgs[ 1979].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -1   1  -1 -sqrt(30)/10
    tc->cgs[ 1983].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -1   2  -1 -sqrt(14)/14
    tc->cgs[ 1989].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -1   3  -1 sqrt(5)/5
    tc->cgs[ 1997].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   0   2  -1   4  -1 sqrt(21)/7
    tc->cgs[ 2003].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   0   0   0 sqrt(5)/5
    tc->cgs[ 2009].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   0   2   0 -sqrt(14)/7
    tc->cgs[ 2023].coeff =   7.1713716560063600E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   0   4   0 3*sqrt(70)/35
    tc->cgs[ 2031].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   1   1   1 sqrt(30)/10
    tc->cgs[ 2035].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   1   2   1 -sqrt(14)/14
    tc->cgs[ 2041].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   1   3   1 -sqrt(5)/5
    tc->cgs[ 2049].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   1   4   1 sqrt(21)/7
    tc->cgs[ 2061].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   2   2   2 sqrt(14)/7
    tc->cgs[ 2067].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   2   3   2 -sqrt(2)/2
    tc->cgs[ 2075].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   0   2   2   4   2 sqrt(42)/14
    tc->cgs[ 2087].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -3   3  -3 sqrt(15)/6
    tc->cgs[ 2095].coeff =   6.7082039324993700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -3   4  -3 3*sqrt(5)/10
    tc->cgs[ 2105].coeff =   3.6514837167011100E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -3   5  -3 sqrt(30)/15
    tc->cgs[ 2118].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -2   2  -2 -sqrt(70)/14
    tc->cgs[ 2132].coeff =   5.8554004376912000E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -2   4  -2 2*sqrt(105)/35
    tc->cgs[ 2142].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -2   5  -2 sqrt(30)/10
    tc->cgs[ 2151].coeff =   4.1403933560541300E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -1   1  -1 sqrt(210)/35
    tc->cgs[ 2155].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -1   2  -1 -sqrt(7)/7
    tc->cgs[ 2161].coeff =  -3.8729833462074200E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -1   3  -1 -sqrt(15)/10
    tc->cgs[ 2169].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -1   4  -1 sqrt(21)/14
    tc->cgs[ 2179].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3  -1   5  -1 sqrt(21)/7
    tc->cgs[ 2188].coeff =   5.0709255283711000E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   0   1   0 3*sqrt(35)/35
    tc->cgs[ 2198].coeff =  -5.1639777949432200E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   0   3   0 -2*sqrt(15)/15
    tc->cgs[ 2216].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   0   5   0 sqrt(210)/21
    tc->cgs[ 2225].coeff =   4.1403933560541300E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   1   1   1 sqrt(210)/35
    tc->cgs[ 2229].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   1   2   1 sqrt(7)/7
    tc->cgs[ 2235].coeff =  -3.8729833462074200E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   1   3   1 -sqrt(15)/10
    tc->cgs[ 2243].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   1   4   1 -sqrt(21)/14
    tc->cgs[ 2253].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   1   5   1 sqrt(21)/7
    tc->cgs[ 2266].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   2   2   2 sqrt(70)/14
    tc->cgs[ 2280].coeff =  -5.8554004376912000E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   2   4   2 -2*sqrt(105)/35
    tc->cgs[ 2290].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   2   5   2 sqrt(30)/10
    tc->cgs[ 2309].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   3   3   3 sqrt(15)/6
    tc->cgs[ 2317].coeff =  -6.7082039324993700E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   3   4   3 -3*sqrt(5)/10
    tc->cgs[ 2327].coeff =   3.6514837167011100E-1; // j1 m1 j2 m2 J M coeff:   2   0   3   3   5   3 sqrt(30)/15
    tc->cgs[ 2337].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   1   0   0   2   1 1
    tc->cgs[ 2341].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   1   1  -1   1   0 sqrt(30)/10
    tc->cgs[ 2345].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   1   1  -1   2   0 sqrt(2)/2
    tc->cgs[ 2351].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   1   1  -1   3   0 sqrt(5)/5
    tc->cgs[ 2358].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   1   1   0   1   1 -sqrt(30)/10
    tc->cgs[ 2362].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2   1   1   0   2   1 sqrt(6)/6
    tc->cgs[ 2368].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   2   1   1   0   3   1 2*sqrt(30)/15
    tc->cgs[ 2379].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   1   1   1   2   2 -sqrt(3)/3
    tc->cgs[ 2385].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   2   1   1   1   3   2 sqrt(6)/3
    tc->cgs[ 2388].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -2   1  -1 sqrt(5)/5
    tc->cgs[ 2392].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -2   2  -1 sqrt(21)/7
    tc->cgs[ 2398].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -2   3  -1 sqrt(30)/10
    tc->cgs[ 2406].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -2   4  -1 sqrt(14)/14
    tc->cgs[ 2412].coeff =  -4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -1   0   0 -sqrt(5)/5
    tc->cgs[ 2414].coeff =  -3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -1   1   0 -sqrt(10)/10
    tc->cgs[ 2418].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -1   2   0 sqrt(14)/14
    tc->cgs[ 2424].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -1   3   0 sqrt(10)/5
    tc->cgs[ 2432].coeff =   4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   2   1   2  -1   4   0 2*sqrt(70)/35
    tc->cgs[ 2440].coeff =  -5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   0   1   1 -sqrt(30)/10
    tc->cgs[ 2444].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   0   2   1 -sqrt(14)/14
    tc->cgs[ 2450].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   0   3   1 sqrt(5)/5
    tc->cgs[ 2458].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   0   4   1 sqrt(21)/7
    tc->cgs[ 2470].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   1   2   2 -sqrt(21)/7
    tc->cgs[ 2484].coeff =   7.5592894601845500E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   1   4   2 2*sqrt(7)/7
    tc->cgs[ 2502].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   2   3   3 -sqrt(2)/2
    tc->cgs[ 2510].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   1   2   2   4   3 sqrt(2)/2
    tc->cgs[ 2516].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -3   2  -2 sqrt(70)/14
    tc->cgs[ 2522].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -3   3  -2 sqrt(15)/6
    tc->cgs[ 2530].coeff =   4.3915503282684000E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -3   4  -2 3*sqrt(105)/70
    tc->cgs[ 2540].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -3   5  -2 sqrt(30)/30
    tc->cgs[ 2549].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -2   1  -1 -sqrt(14)/7
    tc->cgs[ 2559].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -2   3  -1 1/2
    tc->cgs[ 2567].coeff =   5.9160797830996200E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -2   4  -1 sqrt(35)/10
    tc->cgs[ 2577].coeff =   3.3806170189140700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -2   5  -1 2*sqrt(35)/35
    tc->cgs[ 2586].coeff =  -4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -1   1   0 -2*sqrt(70)/35
    tc->cgs[ 2590].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -1   2   0 -sqrt(7)/7
    tc->cgs[ 2596].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -1   3   0 sqrt(30)/30
    tc->cgs[ 2604].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -1   4   0 sqrt(70)/14
    tc->cgs[ 2614].coeff =   4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   2   1   3  -1   5   0 sqrt(105)/21
    tc->cgs[ 2623].coeff =  -2.9277002188456000E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   0   1   1 -sqrt(105)/35
    tc->cgs[ 2627].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   0   2   1 -sqrt(14)/7
    tc->cgs[ 2633].coeff =  -1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   0   3   1 -sqrt(30)/30
    tc->cgs[ 2641].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   0   4   1 sqrt(42)/14
    tc->cgs[ 2651].coeff =   6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   0   5   1 2*sqrt(42)/21
    tc->cgs[ 2664].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   1   2   2 -sqrt(42)/14
    tc->cgs[ 2670].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   1   3   2 -1/2
    tc->cgs[ 2678].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   1   4   2 sqrt(7)/14
    tc->cgs[ 2688].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   1   5   2 sqrt(2)/2
    tc->cgs[ 2707].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   2   3   3 -sqrt(15)/6
    tc->cgs[ 2715].coeff =  -2.2360679774997900E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   2   4   3 -sqrt(5)/10
    tc->cgs[ 2725].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   2   5   3 2*sqrt(30)/15
    tc->cgs[ 2752].coeff =  -7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   3   4   4 -sqrt(15)/5
    tc->cgs[ 2762].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   1   3   3   5   4 sqrt(10)/5
    tc->cgs[ 2772].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   2   0   0   2   2 1
    tc->cgs[ 2776].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2   2   1  -1   1   1 sqrt(15)/5
    tc->cgs[ 2780].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   2   1  -1   2   1 sqrt(3)/3
    tc->cgs[ 2786].coeff =   2.5819888974716100E-1; // j1 m1 j2 m2 J M coeff:   2   2   1  -1   3   1 sqrt(15)/15
    tc->cgs[ 2797].coeff =   8.1649658092772600E-1; // j1 m1 j2 m2 J M coeff:   2   2   1   0   2   2 sqrt(6)/3
    tc->cgs[ 2803].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   2   1   0   3   2 sqrt(3)/3
    tc->cgs[ 2820].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   2   1   1   3   3 1
    tc->cgs[ 2821].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -2   0   0 sqrt(5)/5
    tc->cgs[ 2823].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -2   1   0 sqrt(10)/5
    tc->cgs[ 2827].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -2   2   0 sqrt(14)/7
    tc->cgs[ 2833].coeff =   3.1622776601683800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -2   3   0 sqrt(10)/10
    tc->cgs[ 2841].coeff =   1.1952286093343900E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -2   4   0 sqrt(70)/70
    tc->cgs[ 2849].coeff =   4.4721359549995800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -1   1   1 sqrt(5)/5
    tc->cgs[ 2853].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -1   2   1 sqrt(21)/7
    tc->cgs[ 2859].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -1   3   1 sqrt(30)/10
    tc->cgs[ 2867].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   2   2  -1   4   1 sqrt(14)/14
    tc->cgs[ 2879].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   2   2   0   2   2 sqrt(14)/7
    tc->cgs[ 2885].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2   0   3   2 sqrt(2)/2
    tc->cgs[ 2893].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   2   2   0   4   2 sqrt(42)/14
    tc->cgs[ 2911].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2   1   3   3 sqrt(2)/2
    tc->cgs[ 2919].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   2   2   1   4   3 sqrt(2)/2
    tc->cgs[ 2945].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   2   2   2   4   4 1
    tc->cgs[ 2947].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -3   1  -1 sqrt(21)/7
    tc->cgs[ 2951].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -3   2  -1 sqrt(70)/14
    tc->cgs[ 2957].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -3   3  -1 sqrt(6)/6
    tc->cgs[ 2965].coeff =   2.0701966780270600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -3   4  -1 sqrt(210)/70
    tc->cgs[ 2975].coeff =   6.9006555934235400E-2; // j1 m1 j2 m2 J M coeff:   2   2   3  -3   5  -1 sqrt(210)/210
    tc->cgs[ 2984].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -2   1   0 sqrt(7)/7
    tc->cgs[ 2988].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -2   2   0 sqrt(70)/14
    tc->cgs[ 2994].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -2   3   0 sqrt(3)/3
    tc->cgs[ 3002].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -2   4   0 sqrt(7)/7
    tc->cgs[ 3012].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -2   5   0 sqrt(42)/42
    tc->cgs[ 3021].coeff =   1.6903085094570300E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -1   1   1 sqrt(35)/35
    tc->cgs[ 3025].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -1   2   1 sqrt(42)/14
    tc->cgs[ 3031].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -1   3   1 sqrt(10)/5
    tc->cgs[ 3039].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -1   4   1 sqrt(14)/7
    tc->cgs[ 3049].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   2   3  -1   5   1 sqrt(14)/14
    tc->cgs[ 3062].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   0   2   2 sqrt(14)/14
    tc->cgs[ 3068].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   0   3   2 sqrt(3)/3
    tc->cgs[ 3076].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   0   4   2 sqrt(21)/7
    tc->cgs[ 3086].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   0   5   2 sqrt(6)/6
    tc->cgs[ 3105].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   1   3   3 sqrt(6)/6
    tc->cgs[ 3113].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   1   4   3 sqrt(2)/2
    tc->cgs[ 3123].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   1   5   3 sqrt(3)/3
    tc->cgs[ 3150].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   2   4   4 sqrt(10)/5
    tc->cgs[ 3160].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   2   2   3   2   5   4 sqrt(15)/5
    tc->cgs[ 3197].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   2   2   3   3   5   5 1
    tc->cgs[ 3207].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -3   0   0   3  -3 1
    tc->cgs[ 3230].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -3   1  -1   4  -4 1
    tc->cgs[ 3248].coeff =  -8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   3  -3   1   0   3  -3 -sqrt(3)/2
    tc->cgs[ 3256].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3  -3   1   0   4  -3 1/2
    tc->cgs[ 3268].coeff =   8.4515425472851700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   1   1   2  -2 sqrt(35)/7
    tc->cgs[ 3274].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3  -3   1   1   3  -2 -1/2
    tc->cgs[ 3282].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3  -3   1   1   4  -2 sqrt(7)/14
    tc->cgs[ 3314].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -3   2  -2   5  -5 1
    tc->cgs[ 3341].coeff =  -7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2  -1   4  -4 -sqrt(15)/5
    tc->cgs[ 3351].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2  -1   5  -4 sqrt(10)/5
    tc->cgs[ 3370].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   0   3  -3 sqrt(15)/6
    tc->cgs[ 3378].coeff =  -6.7082039324993700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   0   4  -3 -3*sqrt(5)/10
    tc->cgs[ 3388].coeff =   3.6514837167011100E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   0   5  -3 sqrt(30)/15
    tc->cgs[ 3401].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   1   2  -2 -sqrt(70)/14
    tc->cgs[ 3407].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   1   3  -2 sqrt(15)/6
    tc->cgs[ 3415].coeff =  -4.3915503282684000E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   1   4  -2 -3*sqrt(105)/70
    tc->cgs[ 3425].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   1   5  -2 sqrt(30)/30
    tc->cgs[ 3434].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   2   1  -1 sqrt(21)/7
    tc->cgs[ 3438].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   2   2  -1 -sqrt(70)/14
    tc->cgs[ 3444].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   2   3  -1 sqrt(6)/6
    tc->cgs[ 3452].coeff =  -2.0701966780270600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   2   2   4  -1 -sqrt(210)/70
    tc->cgs[ 3462].coeff =   6.9006555934235400E-2; // j1 m1 j2 m2 J M coeff:   3  -3   2   2   5  -1 sqrt(210)/210
    tc->cgs[ 3505].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -3   3  -3   6  -6 1
    tc->cgs[ 3543].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3  -2   5  -5 -sqrt(2)/2
    tc->cgs[ 3555].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3  -2   6  -5 sqrt(2)/2
    tc->cgs[ 3583].coeff =   5.2223296786709400E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3  -1   4  -4 sqrt(33)/11
    tc->cgs[ 3593].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3  -1   5  -4 -sqrt(2)/2
    tc->cgs[ 3605].coeff =   4.7673129462279600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3  -1   6  -4 sqrt(110)/22
    tc->cgs[ 3625].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   0   3  -3 -sqrt(6)/6
    tc->cgs[ 3633].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   0   4  -3 3*sqrt(22)/22
    tc->cgs[ 3643].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   0   5  -3 -sqrt(3)/3
    tc->cgs[ 3655].coeff =   3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   0   6  -3 sqrt(11)/11
    tc->cgs[ 3669].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   1   2  -2 sqrt(210)/42
    tc->cgs[ 3675].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   1   3  -2 -sqrt(3)/3
    tc->cgs[ 3683].coeff =   5.9215652546379200E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   1   4  -2 3*sqrt(231)/77
    tc->cgs[ 3693].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   1   5  -2 -sqrt(6)/6
    tc->cgs[ 3705].coeff =   1.7407765595569800E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   1   6  -2 sqrt(33)/33
    tc->cgs[ 3715].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   1  -1 -sqrt(21)/14
    tc->cgs[ 3719].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   2  -1 5*sqrt(21)/42
    tc->cgs[ 3725].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   3  -1 -sqrt(3)/3
    tc->cgs[ 3733].coeff =   4.4136741475237500E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   4  -1 sqrt(1155)/77
    tc->cgs[ 3743].coeff =  -2.4397501823713300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   5  -1 -sqrt(105)/42
    tc->cgs[ 3755].coeff =   8.7038827977848900E-2; // j1 m1 j2 m2 J M coeff:   3  -3   3   2   6  -1 sqrt(33)/66
    tc->cgs[ 3763].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   0   0 sqrt(7)/7
    tc->cgs[ 3765].coeff =  -5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   1   0 -3*sqrt(7)/14
    tc->cgs[ 3769].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   2   0 5*sqrt(21)/42
    tc->cgs[ 3775].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   3   0 -sqrt(6)/6
    tc->cgs[ 3783].coeff =   2.4174688920761400E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   4   0 3*sqrt(154)/154
    tc->cgs[ 3793].coeff =  -1.0910894511799600E-1; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   5   0 -sqrt(21)/42
    tc->cgs[ 3805].coeff =   3.2897584747988400E-2; // j1 m1 j2 m2 J M coeff:   3  -3   3   3   6   0 sqrt(231)/462
    tc->cgs[ 3822].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -2   0   0   3  -2 1
    tc->cgs[ 3837].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1  -1   3  -3 1/2
    tc->cgs[ 3845].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1  -1   4  -3 sqrt(3)/2
    tc->cgs[ 3857].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   0   2  -2 -sqrt(105)/21
    tc->cgs[ 3863].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   0   3  -2 -sqrt(3)/3
    tc->cgs[ 3871].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   0   4  -2 sqrt(21)/7
    tc->cgs[ 3883].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   1   2  -1 sqrt(210)/21
    tc->cgs[ 3889].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   1   3  -1 -sqrt(15)/6
    tc->cgs[ 3897].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3  -2   1   1   4  -1 sqrt(21)/14
    tc->cgs[ 3919].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2  -2   4  -4 sqrt(10)/5
    tc->cgs[ 3929].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2  -2   5  -4 sqrt(15)/5
    tc->cgs[ 3948].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2  -1   3  -3 -sqrt(15)/6
    tc->cgs[ 3956].coeff =  -2.2360679774997900E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2  -1   4  -3 -sqrt(5)/10
    tc->cgs[ 3966].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2  -1   5  -3 2*sqrt(30)/15
    tc->cgs[ 3979].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   0   2  -2 sqrt(70)/14
    tc->cgs[ 3993].coeff =  -5.8554004376912000E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   0   4  -2 -2*sqrt(105)/35
    tc->cgs[ 4003].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   0   5  -2 sqrt(30)/10
    tc->cgs[ 4012].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   1   1  -1 -sqrt(14)/7
    tc->cgs[ 4022].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   1   3  -1 1/2
    tc->cgs[ 4030].coeff =  -5.9160797830996200E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   1   4  -1 -sqrt(35)/10
    tc->cgs[ 4040].coeff =   3.3806170189140700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   1   5  -1 2*sqrt(35)/35
    tc->cgs[ 4049].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   2   1   0 sqrt(7)/7
    tc->cgs[ 4053].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   2   2   0 -sqrt(70)/14
    tc->cgs[ 4059].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   2   3   0 sqrt(3)/3
    tc->cgs[ 4067].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   2   4   0 -sqrt(7)/7
    tc->cgs[ 4077].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3  -2   2   2   5   0 sqrt(42)/42
    tc->cgs[ 4108].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -3   5  -5 sqrt(2)/2
    tc->cgs[ 4120].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -3   6  -5 sqrt(2)/2
    tc->cgs[ 4148].coeff =  -6.7419986246324200E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -2   4  -4 -sqrt(55)/11
    tc->cgs[ 4170].coeff =   7.3854894587599600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -2   6  -4 sqrt(66)/11
    tc->cgs[ 4190].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -1   3  -3 sqrt(3)/3
    tc->cgs[ 4198].coeff =  -3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -1   4  -3 -sqrt(11)/11
    tc->cgs[ 4208].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -1   5  -3 -sqrt(6)/6
    tc->cgs[ 4220].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3  -1   6  -3 3*sqrt(22)/22
    tc->cgs[ 4234].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   0   2  -2 -sqrt(105)/21
    tc->cgs[ 4240].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   0   3  -2 sqrt(6)/6
    tc->cgs[ 4248].coeff =   1.3957263155977100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   0   4  -2 sqrt(462)/154
    tc->cgs[ 4258].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   0   5  -2 -sqrt(3)/3
    tc->cgs[ 4270].coeff =   4.9236596391733100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   0   6  -2 2*sqrt(66)/33
    tc->cgs[ 4280].coeff =   4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   1   1  -1 sqrt(35)/14
    tc->cgs[ 4284].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   1   2  -1 -sqrt(35)/14
    tc->cgs[ 4298].coeff =   4.5584230583855200E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   1   4  -1 4*sqrt(77)/77
    tc->cgs[ 4308].coeff =  -5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   1   5  -1 -3*sqrt(7)/14
    tc->cgs[ 4320].coeff =   3.3709993123162100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   1   6  -1 sqrt(55)/22
    tc->cgs[ 4328].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   0   0 -sqrt(7)/7
    tc->cgs[ 4330].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   1   0 sqrt(7)/7
    tc->cgs[ 4340].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   3   0 -sqrt(6)/6
    tc->cgs[ 4348].coeff =   5.6407607481776600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   4   0 sqrt(154)/22
    tc->cgs[ 4358].coeff =  -4.3643578047198500E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   5   0 -2*sqrt(21)/21
    tc->cgs[ 4370].coeff =   1.9738550848793100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   2   6   0 sqrt(231)/77
    tc->cgs[ 4380].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   1   1 -sqrt(21)/14
    tc->cgs[ 4384].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   2   1 5*sqrt(21)/42
    tc->cgs[ 4390].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   3   1 -sqrt(3)/3
    tc->cgs[ 4398].coeff =   4.4136741475237500E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   4   1 sqrt(1155)/77
    tc->cgs[ 4408].coeff =  -2.4397501823713300E-1; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   5   1 -sqrt(105)/42
    tc->cgs[ 4420].coeff =   8.7038827977848900E-2; // j1 m1 j2 m2 J M coeff:   3  -2   3   3   6   1 sqrt(33)/66
    tc->cgs[ 4437].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3  -1   0   0   3  -1 1
    tc->cgs[ 4446].coeff =   2.1821789023599200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1  -1   2  -2 sqrt(21)/21
    tc->cgs[ 4452].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1  -1   3  -2 sqrt(15)/6
    tc->cgs[ 4460].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1  -1   4  -2 sqrt(105)/14
    tc->cgs[ 4472].coeff =  -6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   0   2  -1 -2*sqrt(42)/21
    tc->cgs[ 4478].coeff =  -2.8867513459481300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   0   3  -1 -sqrt(3)/6
    tc->cgs[ 4486].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   0   4  -1 sqrt(105)/14
    tc->cgs[ 4498].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   1   2   0 sqrt(14)/7
    tc->cgs[ 4504].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   1   3   0 -sqrt(2)/2
    tc->cgs[ 4512].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   1   1   4   0 sqrt(42)/14
    tc->cgs[ 4526].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -2   3  -3 sqrt(6)/6
    tc->cgs[ 4534].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -2   4  -3 sqrt(2)/2
    tc->cgs[ 4544].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -2   5  -3 sqrt(3)/3
    tc->cgs[ 4557].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -1   2  -2 -sqrt(42)/14
    tc->cgs[ 4563].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -1   3  -2 -1/2
    tc->cgs[ 4571].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -1   4  -2 sqrt(7)/14
    tc->cgs[ 4581].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2  -1   5  -2 sqrt(2)/2
    tc->cgs[ 4590].coeff =   4.1403933560541300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   0   1  -1 sqrt(210)/35
    tc->cgs[ 4594].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   0   2  -1 sqrt(7)/7
    tc->cgs[ 4600].coeff =  -3.8729833462074200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   0   3  -1 -sqrt(15)/10
    tc->cgs[ 4608].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   0   4  -1 -sqrt(21)/14
    tc->cgs[ 4618].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   0   5  -1 sqrt(21)/7
    tc->cgs[ 4627].coeff =  -4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   1   1   0 -2*sqrt(70)/35
    tc->cgs[ 4631].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   1   2   0 sqrt(7)/7
    tc->cgs[ 4637].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   1   3   0 sqrt(30)/30
    tc->cgs[ 4645].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   1   4   0 -sqrt(70)/14
    tc->cgs[ 4655].coeff =   4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   1   5   0 sqrt(105)/21
    tc->cgs[ 4664].coeff =   1.6903085094570300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   2   1   1 sqrt(35)/35
    tc->cgs[ 4668].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   2   2   1 -sqrt(42)/14
    tc->cgs[ 4674].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   2   3   1 sqrt(10)/5
    tc->cgs[ 4682].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   2   4   1 -sqrt(14)/7
    tc->cgs[ 4692].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   2   2   5   1 sqrt(14)/14
    tc->cgs[ 4713].coeff =   5.2223296786709400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -3   4  -4 sqrt(33)/11
    tc->cgs[ 4723].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -3   5  -4 sqrt(2)/2
    tc->cgs[ 4735].coeff =   4.7673129462279600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -3   6  -4 sqrt(110)/22
    tc->cgs[ 4755].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -2   3  -3 -sqrt(3)/3
    tc->cgs[ 4763].coeff =  -3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -2   4  -3 -sqrt(11)/11
    tc->cgs[ 4773].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -2   5  -3 sqrt(6)/6
    tc->cgs[ 4785].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -2   6  -3 3*sqrt(22)/22
    tc->cgs[ 4799].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -1   2  -2 sqrt(14)/7
    tc->cgs[ 4813].coeff =  -5.0964719143762500E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -1   4  -2 -2*sqrt(385)/77
    tc->cgs[ 4835].coeff =   6.7419986246324200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3  -1   6  -2 sqrt(55)/11
    tc->cgs[ 4845].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   1  -1 -sqrt(42)/14
    tc->cgs[ 4849].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   2  -1 sqrt(42)/42
    tc->cgs[ 4855].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   3  -1 sqrt(6)/6
    tc->cgs[ 4863].coeff =  -3.1209389196618000E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   4  -1 -sqrt(2310)/154
    tc->cgs[ 4873].coeff =  -3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   5  -1 -sqrt(210)/42
    tc->cgs[ 4885].coeff =   6.1545745489666400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   0   6  -1 5*sqrt(66)/66
    tc->cgs[ 4893].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   0   0 sqrt(7)/7
    tc->cgs[ 4895].coeff =  -1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   1   0 -sqrt(7)/14
    tc->cgs[ 4899].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   2   0 -sqrt(21)/14
    tc->cgs[ 4905].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   3   0 sqrt(6)/6
    tc->cgs[ 4913].coeff =   8.0582296402538000E-2; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   4   0 sqrt(154)/154
    tc->cgs[ 4923].coeff =  -5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   5   0 -5*sqrt(21)/42
    tc->cgs[ 4935].coeff =   4.9346377121982700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   1   6   0 5*sqrt(231)/154
    tc->cgs[ 4945].coeff =   4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   2   1   1 sqrt(35)/14
    tc->cgs[ 4949].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   2   2   1 -sqrt(35)/14
    tc->cgs[ 4963].coeff =   4.5584230583855200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   2   4   1 4*sqrt(77)/77
    tc->cgs[ 4973].coeff =  -5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   2   5   1 -3*sqrt(7)/14
    tc->cgs[ 4985].coeff =   3.3709993123162100E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   2   6   1 sqrt(55)/22
    tc->cgs[ 4999].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   3   2   2 sqrt(210)/42
    tc->cgs[ 5005].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   3   3   2 -sqrt(3)/3
    tc->cgs[ 5013].coeff =   5.9215652546379200E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   3   4   2 3*sqrt(231)/77
    tc->cgs[ 5023].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   3   5   2 -sqrt(6)/6
    tc->cgs[ 5035].coeff =   1.7407765595569800E-1; // j1 m1 j2 m2 J M coeff:   3  -1   3   3   6   2 sqrt(33)/33
    tc->cgs[ 5052].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   0   0   0   3   0 1
    tc->cgs[ 5061].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   0   1  -1   2  -1 sqrt(7)/7
    tc->cgs[ 5067].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   0   1  -1   3  -1 sqrt(2)/2
    tc->cgs[ 5075].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   0   1  -1   4  -1 sqrt(70)/14
    tc->cgs[ 5087].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   0   1   0   2   0 -sqrt(21)/7
    tc->cgs[ 5101].coeff =   7.5592894601845500E-1; // j1 m1 j2 m2 J M coeff:   3   0   1   0   4   0 2*sqrt(7)/7
    tc->cgs[ 5113].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   0   1   1   2   1 sqrt(7)/7
    tc->cgs[ 5119].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   0   1   1   3   1 -sqrt(2)/2
    tc->cgs[ 5127].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   0   1   1   4   1 sqrt(70)/14
    tc->cgs[ 5135].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -2   2  -2 sqrt(14)/14
    tc->cgs[ 5141].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -2   3  -2 sqrt(3)/3
    tc->cgs[ 5149].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -2   4  -2 sqrt(21)/7
    tc->cgs[ 5159].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -2   5  -2 sqrt(6)/6
    tc->cgs[ 5168].coeff =  -2.9277002188456000E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -1   1  -1 -sqrt(105)/35
    tc->cgs[ 5172].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -1   2  -1 -sqrt(14)/7
    tc->cgs[ 5178].coeff =  -1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -1   3  -1 -sqrt(30)/30
    tc->cgs[ 5186].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -1   4  -1 sqrt(42)/14
    tc->cgs[ 5196].coeff =   6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   3   0   2  -1   5  -1 2*sqrt(42)/21
    tc->cgs[ 5205].coeff =   5.0709255283711000E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   0   1   0 3*sqrt(35)/35
    tc->cgs[ 5215].coeff =  -5.1639777949432200E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   0   3   0 -2*sqrt(15)/15
    tc->cgs[ 5233].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   0   5   0 sqrt(210)/21
    tc->cgs[ 5242].coeff =  -2.9277002188456000E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   1   1   1 -sqrt(105)/35
    tc->cgs[ 5246].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   1   2   1 sqrt(14)/7
    tc->cgs[ 5252].coeff =  -1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   1   3   1 -sqrt(30)/30
    tc->cgs[ 5260].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   1   4   1 -sqrt(42)/14
    tc->cgs[ 5270].coeff =   6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   1   5   1 2*sqrt(42)/21
    tc->cgs[ 5283].coeff =  -2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   2   2   2 -sqrt(14)/14
    tc->cgs[ 5289].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   2   3   2 sqrt(3)/3
    tc->cgs[ 5297].coeff =  -6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   2   4   2 -sqrt(21)/7
    tc->cgs[ 5307].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   2   2   5   2 sqrt(6)/6
    tc->cgs[ 5320].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -3   3  -3 sqrt(6)/6
    tc->cgs[ 5328].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -3   4  -3 3*sqrt(22)/22
    tc->cgs[ 5338].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -3   5  -3 sqrt(3)/3
    tc->cgs[ 5350].coeff =   3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -3   6  -3 sqrt(11)/11
    tc->cgs[ 5364].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -2   2  -2 -sqrt(105)/21
    tc->cgs[ 5370].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -2   3  -2 -sqrt(6)/6
    tc->cgs[ 5378].coeff =   1.3957263155977100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -2   4  -2 sqrt(462)/154
    tc->cgs[ 5388].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -2   5  -2 sqrt(3)/3
    tc->cgs[ 5400].coeff =   4.9236596391733100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -2   6  -2 2*sqrt(66)/33
    tc->cgs[ 5410].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   1  -1 sqrt(42)/14
    tc->cgs[ 5414].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   2  -1 sqrt(42)/42
    tc->cgs[ 5420].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   3  -1 -sqrt(6)/6
    tc->cgs[ 5428].coeff =  -3.1209389196618000E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   4  -1 -sqrt(2310)/154
    tc->cgs[ 5438].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   5  -1 sqrt(210)/42
    tc->cgs[ 5450].coeff =   6.1545745489666400E-1; // j1 m1 j2 m2 J M coeff:   3   0   3  -1   6  -1 5*sqrt(66)/66
    tc->cgs[ 5458].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   0   0   0 -sqrt(7)/7
    tc->cgs[ 5464].coeff =   4.3643578047198500E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   0   2   0 2*sqrt(21)/21
    tc->cgs[ 5478].coeff =  -4.8349377841522800E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   0   4   0 -3*sqrt(154)/77
    tc->cgs[ 5500].coeff =   6.5795169495976900E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   0   6   0 10*sqrt(231)/231
    tc->cgs[ 5510].coeff =  -4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   1   1 -sqrt(42)/14
    tc->cgs[ 5514].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   2   1 sqrt(42)/42
    tc->cgs[ 5520].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   3   1 sqrt(6)/6
    tc->cgs[ 5528].coeff =  -3.1209389196618000E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   4   1 -sqrt(2310)/154
    tc->cgs[ 5538].coeff =  -3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   5   1 -sqrt(210)/42
    tc->cgs[ 5550].coeff =   6.1545745489666400E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   1   6   1 5*sqrt(66)/66
    tc->cgs[ 5564].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   2   2   2 -sqrt(105)/21
    tc->cgs[ 5570].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   2   3   2 sqrt(6)/6
    tc->cgs[ 5578].coeff =   1.3957263155977100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   2   4   2 sqrt(462)/154
    tc->cgs[ 5588].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   2   5   2 -sqrt(3)/3
    tc->cgs[ 5600].coeff =   4.9236596391733100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   2   6   2 2*sqrt(66)/33
    tc->cgs[ 5620].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   3   3   3 -sqrt(6)/6
    tc->cgs[ 5628].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   3   4   3 3*sqrt(22)/22
    tc->cgs[ 5638].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   3   5   3 -sqrt(3)/3
    tc->cgs[ 5650].coeff =   3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3   0   3   3   6   3 sqrt(11)/11
    tc->cgs[ 5667].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   1   0   0   3   1 1
    tc->cgs[ 5676].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   1   1  -1   2   0 sqrt(14)/7
    tc->cgs[ 5682].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   1   1  -1   3   0 sqrt(2)/2
    tc->cgs[ 5690].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   1   1  -1   4   0 sqrt(42)/14
    tc->cgs[ 5702].coeff =  -6.1721339984836800E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   0   2   1 -2*sqrt(42)/21
    tc->cgs[ 5708].coeff =   2.8867513459481300E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   0   3   1 sqrt(3)/6
    tc->cgs[ 5716].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   0   4   1 sqrt(105)/14
    tc->cgs[ 5728].coeff =   2.1821789023599200E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   1   2   2 sqrt(21)/21
    tc->cgs[ 5734].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   1   3   2 -sqrt(15)/6
    tc->cgs[ 5742].coeff =   7.3192505471140000E-1; // j1 m1 j2 m2 J M coeff:   3   1   1   1   4   2 sqrt(105)/14
    tc->cgs[ 5746].coeff =   1.6903085094570300E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -2   1  -1 sqrt(35)/35
    tc->cgs[ 5750].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -2   2  -1 sqrt(42)/14
    tc->cgs[ 5756].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -2   3  -1 sqrt(10)/5
    tc->cgs[ 5764].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -2   4  -1 sqrt(14)/7
    tc->cgs[ 5774].coeff =   2.6726124191242400E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -2   5  -1 sqrt(14)/14
    tc->cgs[ 5783].coeff =  -4.7809144373375700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -1   1   0 -2*sqrt(70)/35
    tc->cgs[ 5787].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -1   2   0 -sqrt(7)/7
    tc->cgs[ 5793].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -1   3   0 sqrt(30)/30
    tc->cgs[ 5801].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -1   4   0 sqrt(70)/14
    tc->cgs[ 5811].coeff =   4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2  -1   5   0 sqrt(105)/21
    tc->cgs[ 5820].coeff =   4.1403933560541300E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   0   1   1 sqrt(210)/35
    tc->cgs[ 5824].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   0   2   1 -sqrt(7)/7
    tc->cgs[ 5830].coeff =  -3.8729833462074200E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   0   3   1 -sqrt(15)/10
    tc->cgs[ 5838].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   0   4   1 sqrt(21)/14
    tc->cgs[ 5848].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   0   5   1 sqrt(21)/7
    tc->cgs[ 5861].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   1   2   2 sqrt(42)/14
    tc->cgs[ 5867].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   1   3   2 -1/2
    tc->cgs[ 5875].coeff =  -1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   1   4   2 -sqrt(7)/14
    tc->cgs[ 5885].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   1   5   2 sqrt(2)/2
    tc->cgs[ 5904].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   2   3   3 sqrt(6)/6
    tc->cgs[ 5912].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   2   4   3 -sqrt(2)/2
    tc->cgs[ 5922].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   1   2   2   5   3 sqrt(3)/3
    tc->cgs[ 5929].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -3   2  -2 sqrt(210)/42
    tc->cgs[ 5935].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -3   3  -2 sqrt(3)/3
    tc->cgs[ 5943].coeff =   5.9215652546379200E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -3   4  -2 3*sqrt(231)/77
    tc->cgs[ 5953].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -3   5  -2 sqrt(6)/6
    tc->cgs[ 5965].coeff =   1.7407765595569800E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -3   6  -2 sqrt(33)/33
    tc->cgs[ 5975].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -2   1  -1 -sqrt(35)/14
    tc->cgs[ 5979].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -2   2  -1 -sqrt(35)/14
    tc->cgs[ 5993].coeff =   4.5584230583855200E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -2   4  -1 4*sqrt(77)/77
    tc->cgs[ 6003].coeff =   5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -2   5  -1 3*sqrt(7)/14
    tc->cgs[ 6015].coeff =   3.3709993123162100E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -2   6  -1 sqrt(55)/22
    tc->cgs[ 6023].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   0   0 sqrt(7)/7
    tc->cgs[ 6025].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   1   0 sqrt(7)/14
    tc->cgs[ 6029].coeff =  -3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   2   0 -sqrt(21)/14
    tc->cgs[ 6035].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   3   0 -sqrt(6)/6
    tc->cgs[ 6043].coeff =   8.0582296402538000E-2; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   4   0 sqrt(154)/154
    tc->cgs[ 6053].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   5   0 5*sqrt(21)/42
    tc->cgs[ 6065].coeff =   4.9346377121982700E-1; // j1 m1 j2 m2 J M coeff:   3   1   3  -1   6   0 5*sqrt(231)/154
    tc->cgs[ 6075].coeff =   4.6291004988627600E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   1   1 sqrt(42)/14
    tc->cgs[ 6079].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   2   1 sqrt(42)/42
    tc->cgs[ 6085].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   3   1 -sqrt(6)/6
    tc->cgs[ 6093].coeff =  -3.1209389196618000E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   4   1 -sqrt(2310)/154
    tc->cgs[ 6103].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   5   1 sqrt(210)/42
    tc->cgs[ 6115].coeff =   6.1545745489666400E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   0   6   1 5*sqrt(66)/66
    tc->cgs[ 6129].coeff =   5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   1   2   2 sqrt(14)/7
    tc->cgs[ 6143].coeff =  -5.0964719143762500E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   1   4   2 -2*sqrt(385)/77
    tc->cgs[ 6165].coeff =   6.7419986246324200E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   1   6   2 sqrt(55)/11
    tc->cgs[ 6185].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   2   3   3 sqrt(3)/3
    tc->cgs[ 6193].coeff =  -3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   2   4   3 -sqrt(11)/11
    tc->cgs[ 6203].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   2   5   3 -sqrt(6)/6
    tc->cgs[ 6215].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   2   6   3 3*sqrt(22)/22
    tc->cgs[ 6243].coeff =   5.2223296786709400E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   3   4   4 sqrt(33)/11
    tc->cgs[ 6253].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   3   5   4 -sqrt(2)/2
    tc->cgs[ 6265].coeff =   4.7673129462279600E-1; // j1 m1 j2 m2 J M coeff:   3   1   3   3   6   4 sqrt(110)/22
    tc->cgs[ 6282].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   2   0   0   3   2 1
    tc->cgs[ 6291].coeff =   6.9006555934235400E-1; // j1 m1 j2 m2 J M coeff:   3   2   1  -1   2   1 sqrt(210)/21
    tc->cgs[ 6297].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3   2   1  -1   3   1 sqrt(15)/6
    tc->cgs[ 6305].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3   2   1  -1   4   1 sqrt(21)/14
    tc->cgs[ 6317].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3   2   1   0   2   2 -sqrt(105)/21
    tc->cgs[ 6323].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   2   1   0   3   2 sqrt(3)/3
    tc->cgs[ 6331].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   2   1   0   4   2 sqrt(21)/7
    tc->cgs[ 6349].coeff =  -5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3   2   1   1   3   3 -1/2
    tc->cgs[ 6357].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   3   2   1   1   4   3 sqrt(3)/2
    tc->cgs[ 6361].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -2   1   0 sqrt(7)/7
    tc->cgs[ 6365].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -2   2   0 sqrt(70)/14
    tc->cgs[ 6371].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -2   3   0 sqrt(3)/3
    tc->cgs[ 6379].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -2   4   0 sqrt(7)/7
    tc->cgs[ 6389].coeff =   1.5430334996209200E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -2   5   0 sqrt(42)/42
    tc->cgs[ 6398].coeff =  -5.3452248382484900E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -1   1   1 -sqrt(14)/7
    tc->cgs[ 6408].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -1   3   1 1/2
    tc->cgs[ 6416].coeff =   5.9160797830996200E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -1   4   1 sqrt(35)/10
    tc->cgs[ 6426].coeff =   3.3806170189140700E-1; // j1 m1 j2 m2 J M coeff:   3   2   2  -1   5   1 2*sqrt(35)/35
    tc->cgs[ 6439].coeff =  -5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   0   2   2 -sqrt(70)/14
    tc->cgs[ 6453].coeff =   5.8554004376912000E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   0   4   2 2*sqrt(105)/35
    tc->cgs[ 6463].coeff =   5.4772255750516600E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   0   5   2 sqrt(30)/10
    tc->cgs[ 6482].coeff =  -6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   1   3   3 -sqrt(15)/6
    tc->cgs[ 6490].coeff =   2.2360679774997900E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   1   4   3 sqrt(5)/10
    tc->cgs[ 6500].coeff =   7.3029674334022100E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   1   5   3 2*sqrt(30)/15
    tc->cgs[ 6527].coeff =  -6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   2   4   4 -sqrt(10)/5
    tc->cgs[ 6537].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   3   2   2   2   5   4 sqrt(15)/5
    tc->cgs[ 6540].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   1  -1 sqrt(21)/14
    tc->cgs[ 6544].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   2  -1 5*sqrt(21)/42
    tc->cgs[ 6550].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   3  -1 sqrt(3)/3
    tc->cgs[ 6558].coeff =   4.4136741475237500E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   4  -1 sqrt(1155)/77
    tc->cgs[ 6568].coeff =   2.4397501823713300E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   5  -1 sqrt(105)/42
    tc->cgs[ 6580].coeff =   8.7038827977848900E-2; // j1 m1 j2 m2 J M coeff:   3   2   3  -3   6  -1 sqrt(33)/66
    tc->cgs[ 6588].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   0   0 -sqrt(7)/7
    tc->cgs[ 6590].coeff =  -3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   1   0 -sqrt(7)/7
    tc->cgs[ 6600].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   3   0 sqrt(6)/6
    tc->cgs[ 6608].coeff =   5.6407607481776600E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   4   0 sqrt(154)/22
    tc->cgs[ 6618].coeff =   4.3643578047198500E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   5   0 2*sqrt(21)/21
    tc->cgs[ 6630].coeff =   1.9738550848793100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -2   6   0 sqrt(231)/77
    tc->cgs[ 6640].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -1   1   1 -sqrt(35)/14
    tc->cgs[ 6644].coeff =  -4.2257712736425800E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -1   2   1 -sqrt(35)/14
    tc->cgs[ 6658].coeff =   4.5584230583855200E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -1   4   1 4*sqrt(77)/77
    tc->cgs[ 6668].coeff =   5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -1   5   1 3*sqrt(7)/14
    tc->cgs[ 6680].coeff =   3.3709993123162100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3  -1   6   1 sqrt(55)/22
    tc->cgs[ 6694].coeff =  -4.8795003647426700E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   0   2   2 -sqrt(105)/21
    tc->cgs[ 6700].coeff =  -4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   0   3   2 -sqrt(6)/6
    tc->cgs[ 6708].coeff =   1.3957263155977100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   0   4   2 sqrt(462)/154
    tc->cgs[ 6718].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   0   5   2 sqrt(3)/3
    tc->cgs[ 6730].coeff =   4.9236596391733100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   0   6   2 2*sqrt(66)/33
    tc->cgs[ 6750].coeff =  -5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   1   3   3 -sqrt(3)/3
    tc->cgs[ 6758].coeff =  -3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   1   4   3 -sqrt(11)/11
    tc->cgs[ 6768].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   1   5   3 sqrt(6)/6
    tc->cgs[ 6780].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   1   6   3 3*sqrt(22)/22
    tc->cgs[ 6808].coeff =  -6.7419986246324200E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   2   4   4 -sqrt(55)/11
    tc->cgs[ 6830].coeff =   7.3854894587599600E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   2   6   4 sqrt(66)/11
    tc->cgs[ 6868].coeff =  -7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   3   5   5 -sqrt(2)/2
    tc->cgs[ 6880].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   2   3   3   6   5 sqrt(2)/2
    tc->cgs[ 6897].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   3   0   0   3   3 1
    tc->cgs[ 6906].coeff =   8.4515425472851700E-1; // j1 m1 j2 m2 J M coeff:   3   3   1  -1   2   2 sqrt(35)/7
    tc->cgs[ 6912].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3   3   1  -1   3   2 1/2
    tc->cgs[ 6920].coeff =   1.8898223650461400E-1; // j1 m1 j2 m2 J M coeff:   3   3   1  -1   4   2 sqrt(7)/14
    tc->cgs[ 6938].coeff =   8.6602540378443900E-1; // j1 m1 j2 m2 J M coeff:   3   3   1   0   3   3 sqrt(3)/2
    tc->cgs[ 6946].coeff =   5.0000000000000000E-1; // j1 m1 j2 m2 J M coeff:   3   3   1   0   4   3 1/2
    tc->cgs[ 6972].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   3   1   1   4   4 1
    tc->cgs[ 6976].coeff =   6.5465367070797700E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -2   1   1 sqrt(21)/7
    tc->cgs[ 6980].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -2   2   1 sqrt(70)/14
    tc->cgs[ 6986].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -2   3   1 sqrt(6)/6
    tc->cgs[ 6994].coeff =   2.0701966780270600E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -2   4   1 sqrt(210)/70
    tc->cgs[ 7004].coeff =   6.9006555934235400E-2; // j1 m1 j2 m2 J M coeff:   3   3   2  -2   5   1 sqrt(210)/210
    tc->cgs[ 7017].coeff =   5.9761430466719700E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -1   2   2 sqrt(70)/14
    tc->cgs[ 7023].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -1   3   2 sqrt(15)/6
    tc->cgs[ 7031].coeff =   4.3915503282684000E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -1   4   2 3*sqrt(105)/70
    tc->cgs[ 7041].coeff =   1.8257418583505500E-1; // j1 m1 j2 m2 J M coeff:   3   3   2  -1   5   2 sqrt(30)/30
    tc->cgs[ 7060].coeff =   6.4549722436790300E-1; // j1 m1 j2 m2 J M coeff:   3   3   2   0   3   3 sqrt(15)/6
    tc->cgs[ 7068].coeff =   6.7082039324993700E-1; // j1 m1 j2 m2 J M coeff:   3   3   2   0   4   3 3*sqrt(5)/10
    tc->cgs[ 7078].coeff =   3.6514837167011100E-1; // j1 m1 j2 m2 J M coeff:   3   3   2   0   5   3 sqrt(30)/15
    tc->cgs[ 7105].coeff =   7.7459666924148300E-1; // j1 m1 j2 m2 J M coeff:   3   3   2   1   4   4 sqrt(15)/5
    tc->cgs[ 7115].coeff =   6.3245553203367600E-1; // j1 m1 j2 m2 J M coeff:   3   3   2   1   5   4 sqrt(10)/5
    tc->cgs[ 7152].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   3   2   2   5   5 1
    tc->cgs[ 7153].coeff =   3.7796447300922700E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   0   0 sqrt(7)/7
    tc->cgs[ 7155].coeff =   5.6694670951384100E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   1   0 3*sqrt(7)/14
    tc->cgs[ 7159].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   2   0 5*sqrt(21)/42
    tc->cgs[ 7165].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   3   0 sqrt(6)/6
    tc->cgs[ 7173].coeff =   2.4174688920761400E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   4   0 3*sqrt(154)/154
    tc->cgs[ 7183].coeff =   1.0910894511799600E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   5   0 sqrt(21)/42
    tc->cgs[ 7195].coeff =   3.2897584747988400E-2; // j1 m1 j2 m2 J M coeff:   3   3   3  -3   6   0 sqrt(231)/462
    tc->cgs[ 7205].coeff =   3.2732683535398900E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   1   1 sqrt(21)/14
    tc->cgs[ 7209].coeff =   5.4554472558998100E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   2   1 5*sqrt(21)/42
    tc->cgs[ 7215].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   3   1 sqrt(3)/3
    tc->cgs[ 7223].coeff =   4.4136741475237500E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   4   1 sqrt(1155)/77
    tc->cgs[ 7233].coeff =   2.4397501823713300E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   5   1 sqrt(105)/42
    tc->cgs[ 7245].coeff =   8.7038827977848900E-2; // j1 m1 j2 m2 J M coeff:   3   3   3  -2   6   1 sqrt(33)/66
    tc->cgs[ 7259].coeff =   3.4503277967117700E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -1   2   2 sqrt(210)/42
    tc->cgs[ 7265].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -1   3   2 sqrt(3)/3
    tc->cgs[ 7273].coeff =   5.9215652546379200E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -1   4   2 3*sqrt(231)/77
    tc->cgs[ 7283].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -1   5   2 sqrt(6)/6
    tc->cgs[ 7295].coeff =   1.7407765595569800E-1; // j1 m1 j2 m2 J M coeff:   3   3   3  -1   6   2 sqrt(33)/33
    tc->cgs[ 7315].coeff =   4.0824829046386300E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   0   3   3 sqrt(6)/6
    tc->cgs[ 7323].coeff =   6.3960214906683100E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   0   4   3 3*sqrt(22)/22
    tc->cgs[ 7333].coeff =   5.7735026918962600E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   0   5   3 sqrt(3)/3
    tc->cgs[ 7345].coeff =   3.0151134457776400E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   0   6   3 sqrt(11)/11
    tc->cgs[ 7373].coeff =   5.2223296786709400E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   1   4   4 sqrt(33)/11
    tc->cgs[ 7383].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   1   5   4 sqrt(2)/2
    tc->cgs[ 7395].coeff =   4.7673129462279600E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   1   6   4 sqrt(110)/22
    tc->cgs[ 7433].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   2   5   5 sqrt(2)/2
    tc->cgs[ 7445].coeff =   7.0710678118654800E-1; // j1 m1 j2 m2 J M coeff:   3   3   3   2   6   5 sqrt(2)/2
    tc->cgs[ 7495].coeff =   1.0000000000000000E+0; // j1 m1 j2 m2 J M coeff:   3   3   3   3   6   6 1
    // clang-format on

    return;
}

struct TestCase_math_clebschGordanM0 : public TestCase
{
    struct SingleClebschGordanM0
    {
        Int  j1;
        Int  j2;
        Int  J;
        Real coeff;
    };

    Int                                jmax;
    Real                               tolerance = 10 * std::numeric_limits<Real>::epsilon();
    std::vector<SingleClebschGordanM0> cgs;

    TestCase_math_clebschGordanM0(String name) : TestCase(name) {};
};

template<>
inline void TestCaseContainer<TestCase_math_clebschGordanM0>::setup()
{
    TestCase_math_clebschGordanM0* tc = nullptr;
    using SingleClebschGordanM0 = TestCase_math_clebschGordanM0::SingleClebschGordanM0;

    Int jmax = 6;
    testCases.push_back(TestCase_math_clebschGordanM0("j1max = j2max = " + std::to_string(jmax)));
    tc = &(testCases.back());
    tc->jmax = jmax;
    for (Int j1 = 0; j1 <= jmax; ++j1)
    {
        for (Int j2 = 0; j2 <= jmax; ++j2)
        {
            for (Int J = 0; J <= j1 + j2; ++J)
            {
                tc->cgs.push_back(SingleClebschGordanM0({j1, j2, J, 0.0}));
            }
        }
    }

    // clang-format off
    tc->cgs[    0].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   0   0 1
    tc->cgs[    2].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   1   1 1
    tc->cgs[    5].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   2   2 1
    tc->cgs[    9].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   3   3 1
    tc->cgs[   14].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   4   4 1
    tc->cgs[   20].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   5   5 1
    tc->cgs[   27].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   0   6   6 1
    tc->cgs[   29].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   1   0   1 1
    tc->cgs[   30].coeff =  -5.7735026918962600E-1; // j1 j2 J coeff:   1   1   0 -sqrt(3)/3
    tc->cgs[   32].coeff =   8.1649658092772600E-1; // j1 j2 J coeff:   1   1   2 sqrt(6)/3
    tc->cgs[   34].coeff =  -6.3245553203367600E-1; // j1 j2 J coeff:   1   2   1 -sqrt(10)/5
    tc->cgs[   36].coeff =   7.7459666924148300E-1; // j1 j2 J coeff:   1   2   3 sqrt(15)/5
    tc->cgs[   39].coeff =  -6.5465367070797700E-1; // j1 j2 J coeff:   1   3   2 -sqrt(21)/7
    tc->cgs[   41].coeff =   7.5592894601845500E-1; // j1 j2 J coeff:   1   3   4 2*sqrt(7)/7
    tc->cgs[   45].coeff =  -6.6666666666666700E-1; // j1 j2 J coeff:   1   4   3 -2/3
    tc->cgs[   47].coeff =   7.4535599249993000E-1; // j1 j2 J coeff:   1   4   5 sqrt(5)/3
    tc->cgs[   52].coeff =  -6.7419986246324200E-1; // j1 j2 J coeff:   1   5   4 -sqrt(55)/11
    tc->cgs[   54].coeff =   7.3854894587599600E-1; // j1 j2 J coeff:   1   5   6 sqrt(66)/11
    tc->cgs[   60].coeff =  -6.7936622048675700E-1; // j1 j2 J coeff:   1   6   5 -sqrt(78)/13
    tc->cgs[   62].coeff =   7.3379938570534300E-1; // j1 j2 J coeff:   1   6   7 sqrt(91)/13
    tc->cgs[   65].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   2   0   2 1
    tc->cgs[   67].coeff =  -6.3245553203367600E-1; // j1 j2 J coeff:   2   1   1 -sqrt(10)/5
    tc->cgs[   69].coeff =   7.7459666924148300E-1; // j1 j2 J coeff:   2   1   3 sqrt(15)/5
    tc->cgs[   70].coeff =   4.4721359549995800E-1; // j1 j2 J coeff:   2   2   0 sqrt(5)/5
    tc->cgs[   72].coeff =  -5.3452248382484900E-1; // j1 j2 J coeff:   2   2   2 -sqrt(14)/7
    tc->cgs[   74].coeff =   7.1713716560063600E-1; // j1 j2 J coeff:   2   2   4 3*sqrt(70)/35
    tc->cgs[   76].coeff =   5.0709255283711000E-1; // j1 j2 J coeff:   2   3   1 3*sqrt(35)/35
    tc->cgs[   78].coeff =  -5.1639777949432200E-1; // j1 j2 J coeff:   2   3   3 -2*sqrt(15)/15
    tc->cgs[   80].coeff =   6.9006555934235400E-1; // j1 j2 J coeff:   2   3   5 sqrt(210)/21
    tc->cgs[   83].coeff =   5.3452248382484900E-1; // j1 j2 J coeff:   2   4   2 sqrt(14)/7
    tc->cgs[   85].coeff =  -5.0964719143762500E-1; // j1 j2 J coeff:   2   4   4 -2*sqrt(385)/77
    tc->cgs[   87].coeff =   6.7419986246324200E-1; // j1 j2 J coeff:   2   4   6 sqrt(55)/11
    tc->cgs[   91].coeff =   5.5048188256318000E-1; // j1 j2 J coeff:   2   5   3 sqrt(330)/33
    tc->cgs[   93].coeff =  -5.0636968354183300E-1; // j1 j2 J coeff:   2   5   5 -sqrt(390)/39
    tc->cgs[   95].coeff =   6.6374651830306500E-1; // j1 j2 J coeff:   2   5   7 3*sqrt(1001)/143
    tc->cgs[  100].coeff =   5.6096819400507400E-1; // j1 j2 J coeff:   2   6   4 3*sqrt(715)/143
    tc->cgs[  102].coeff =  -5.0452497910951300E-1; // j1 j2 J coeff:   2   6   6 -sqrt(770)/55
    tc->cgs[  104].coeff =   6.5633012331389400E-1; // j1 j2 J coeff:   2   6   8 2*sqrt(455)/65
    tc->cgs[  108].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   3   0   3 1
    tc->cgs[  111].coeff =  -6.5465367070797700E-1; // j1 j2 J coeff:   3   1   2 -sqrt(21)/7
    tc->cgs[  113].coeff =   7.5592894601845500E-1; // j1 j2 J coeff:   3   1   4 2*sqrt(7)/7
    tc->cgs[  115].coeff =   5.0709255283711000E-1; // j1 j2 J coeff:   3   2   1 3*sqrt(35)/35
    tc->cgs[  117].coeff =  -5.1639777949432200E-1; // j1 j2 J coeff:   3   2   3 -2*sqrt(15)/15
    tc->cgs[  119].coeff =   6.9006555934235400E-1; // j1 j2 J coeff:   3   2   5 sqrt(210)/21
    tc->cgs[  120].coeff =  -3.7796447300922700E-1; // j1 j2 J coeff:   3   3   0 -sqrt(7)/7
    tc->cgs[  122].coeff =   4.3643578047198500E-1; // j1 j2 J coeff:   3   3   2 2*sqrt(21)/21
    tc->cgs[  124].coeff =  -4.8349377841522800E-1; // j1 j2 J coeff:   3   3   4 -3*sqrt(154)/77
    tc->cgs[  126].coeff =   6.5795169495976900E-1; // j1 j2 J coeff:   3   3   6 10*sqrt(231)/231
    tc->cgs[  128].coeff =  -4.3643578047198500E-1; // j1 j2 J coeff:   3   4   1 -2*sqrt(21)/21
    tc->cgs[  130].coeff =   4.2640143271122100E-1; // j1 j2 J coeff:   3   4   3 sqrt(22)/11
    tc->cgs[  132].coeff =  -4.6880723093849500E-1; // j1 j2 J coeff:   3   4   5 -2*sqrt(455)/91
    tc->cgs[  134].coeff =   6.3869038502658500E-1; // j1 j2 J coeff:   3   4   7 5*sqrt(3003)/429
    tc->cgs[  137].coeff =  -4.6524210519923500E-1; // j1 j2 J coeff:   3   5   2 -5*sqrt(462)/231
    tc->cgs[  139].coeff =   4.2405209564413200E-1; // j1 j2 J coeff:   3   5   4 6*sqrt(5005)/1001
    tc->cgs[  141].coeff =  -4.6056618647183800E-1; // j1 j2 J coeff:   3   5   6 -sqrt(231)/33
    tc->cgs[  143].coeff =   6.2578621877474400E-1; // j1 j2 J coeff:   3   5   8 2*sqrt(2002)/143
    tc->cgs[  147].coeff =  -4.8280454958526800E-1; // j1 j2 J coeff:   3   6   3 -10*sqrt(429)/429
    tc->cgs[  149].coeff =   4.2365927286816200E-1; // j1 j2 J coeff:   3   6   5 sqrt(273)/39
    tc->cgs[  151].coeff =  -4.5532635512896800E-1; // j1 j2 J coeff:   3   6   7 -6*sqrt(34034)/2431
    tc->cgs[  153].coeff =   6.1651479928510800E-1; // j1 j2 J coeff:   3   6   9 2*sqrt(4641)/221
    tc->cgs[  158].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   4   0   4 1
    tc->cgs[  162].coeff =  -6.6666666666666700E-1; // j1 j2 J coeff:   4   1   3 -2/3
    tc->cgs[  164].coeff =   7.4535599249993000E-1; // j1 j2 J coeff:   4   1   5 sqrt(5)/3
    tc->cgs[  167].coeff =   5.3452248382484900E-1; // j1 j2 J coeff:   4   2   2 sqrt(14)/7
    tc->cgs[  169].coeff =  -5.0964719143762500E-1; // j1 j2 J coeff:   4   2   4 -2*sqrt(385)/77
    tc->cgs[  171].coeff =   6.7419986246324200E-1; // j1 j2 J coeff:   4   2   6 sqrt(55)/11
    tc->cgs[  173].coeff =  -4.3643578047198500E-1; // j1 j2 J coeff:   4   3   1 -2*sqrt(21)/21
    tc->cgs[  175].coeff =   4.2640143271122100E-1; // j1 j2 J coeff:   4   3   3 sqrt(22)/11
    tc->cgs[  177].coeff =  -4.6880723093849500E-1; // j1 j2 J coeff:   4   3   5 -2*sqrt(455)/91
    tc->cgs[  179].coeff =   6.3869038502658500E-1; // j1 j2 J coeff:   4   3   7 5*sqrt(3003)/429
    tc->cgs[  180].coeff =   3.3333333333333300E-1; // j1 j2 J coeff:   4   4   0 1/3
    tc->cgs[  182].coeff =  -3.7986858819879300E-1; // j1 j2 J coeff:   4   4   2 -10*sqrt(77)/231
    tc->cgs[  184].coeff =   4.0229114064090700E-1; // j1 j2 J coeff:   4   4   4 9*sqrt(2002)/1001
    tc->cgs[  186].coeff =  -4.4946657497549500E-1; // j1 j2 J coeff:   4   4   6 -2*sqrt(55)/33
    tc->cgs[  188].coeff =   6.1703353290593600E-1; // j1 j2 J coeff:   4   4   8 7*sqrt(1430)/429
    tc->cgs[  190].coeff =   3.8924947208076100E-1; // j1 j2 J coeff:   4   5   1 sqrt(165)/33
    tc->cgs[  192].coeff =  -3.7397879600338300E-1; // j1 j2 J coeff:   4   5   3 -2*sqrt(715)/143
    tc->cgs[  194].coeff =   3.9223227027636800E-1; // j1 j2 J coeff:   4   5   5 sqrt(26)/13
    tc->cgs[  196].coeff =  -4.3813798950473400E-1; // j1 j2 J coeff:   4   5   7 -10*sqrt(102102)/7293
    tc->cgs[  198].coeff =   6.0234015052236400E-1; // j1 j2 J coeff:   4   5   9 21*sqrt(4862)/2431
    tc->cgs[  201].coeff =   4.1812100500354500E-1; // j1 j2 J coeff:   4   6   2 5*sqrt(143)/143
    tc->cgs[  203].coeff =  -3.7397879600338300E-1; // j1 j2 J coeff:   4   6   4 -2*sqrt(715)/143
    tc->cgs[  205].coeff =   3.8695299497594700E-1; // j1 j2 J coeff:   4   6   6 2*sqrt(1309)/187
    tc->cgs[  207].coeff =  -4.3069561387887500E-1; // j1 j2 J coeff:   4   6   8 -6*sqrt(38038)/2717
    tc->cgs[  209].coeff =   5.9167842041038500E-1; // j1 j2 J coeff:   4   6  10 7*sqrt(125970)/4199
    tc->cgs[  215].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   5   0   5 1
    tc->cgs[  220].coeff =  -6.7419986246324200E-1; // j1 j2 J coeff:   5   1   4 -sqrt(55)/11
    tc->cgs[  222].coeff =   7.3854894587599600E-1; // j1 j2 J coeff:   5   1   6 sqrt(66)/11
    tc->cgs[  226].coeff =   5.5048188256318000E-1; // j1 j2 J coeff:   5   2   3 sqrt(330)/33
    tc->cgs[  228].coeff =  -5.0636968354183300E-1; // j1 j2 J coeff:   5   2   5 -sqrt(390)/39
    tc->cgs[  230].coeff =   6.6374651830306500E-1; // j1 j2 J coeff:   5   2   7 3*sqrt(1001)/143
    tc->cgs[  233].coeff =  -4.6524210519923500E-1; // j1 j2 J coeff:   5   3   2 -5*sqrt(462)/231
    tc->cgs[  235].coeff =   4.2405209564413200E-1; // j1 j2 J coeff:   5   3   4 6*sqrt(5005)/1001
    tc->cgs[  237].coeff =  -4.6056618647183800E-1; // j1 j2 J coeff:   5   3   6 -sqrt(231)/33
    tc->cgs[  239].coeff =   6.2578621877474400E-1; // j1 j2 J coeff:   5   3   8 2*sqrt(2002)/143
    tc->cgs[  241].coeff =   3.8924947208076100E-1; // j1 j2 J coeff:   5   4   1 sqrt(165)/33
    tc->cgs[  243].coeff =  -3.7397879600338300E-1; // j1 j2 J coeff:   5   4   3 -2*sqrt(715)/143
    tc->cgs[  245].coeff =   3.9223227027636800E-1; // j1 j2 J coeff:   5   4   5 sqrt(26)/13
    tc->cgs[  247].coeff =  -4.3813798950473400E-1; // j1 j2 J coeff:   5   4   7 -10*sqrt(102102)/7293
    tc->cgs[  249].coeff =   6.0234015052236400E-1; // j1 j2 J coeff:   5   4   9 21*sqrt(4862)/2431
    tc->cgs[  250].coeff =  -3.0151134457776400E-1; // j1 j2 J coeff:   5   5   0 -sqrt(11)/11
    tc->cgs[  252].coeff =   3.4139437099945900E-1; // j1 j2 J coeff:   5   5   2 5*sqrt(858)/429
    tc->cgs[  254].coeff =  -3.5478743759345000E-1; // j1 j2 J coeff:   5   5   4 -3*sqrt(286)/143
    tc->cgs[  256].coeff =   3.7762745602468100E-1; // j1 j2 J coeff:   5   5   6 4*sqrt(2805)/561
    tc->cgs[  258].coeff =  -4.2467160232308200E-1; // j1 j2 J coeff:   5   5   8 -7*sqrt(27170)/2717
    tc->cgs[  260].coeff =   5.8627485133113200E-1; // j1 j2 J coeff:   5   5  10 126*sqrt(46189)/46189
    tc->cgs[  262].coeff =  -3.5478743759345000E-1; // j1 j2 J coeff:   5   6   1 -3*sqrt(286)/143
    tc->cgs[  264].coeff =   3.3796318470968700E-1; // j1 j2 J coeff:   5   6   3 7*sqrt(429)/429
    tc->cgs[  266].coeff =  -3.4736673714593700E-1; // j1 j2 J coeff:   5   6   5 -4*sqrt(3315)/663
    tc->cgs[  268].coeff =   3.6931844203655700E-1; // j1 j2 J coeff:   5   6   7 30*sqrt(323323)/46189
    tc->cgs[  270].coeff =  -4.1565419288457300E-1; // j1 j2 J coeff:   5   6   9 -2*sqrt(255255)/2431
    tc->cgs[  272].coeff =   5.7452466451886100E-1; // j1 j2 J coeff:   5   6  11 3*sqrt(646646)/4199
    tc->cgs[  279].coeff =   1.0000000000000000E+0; // j1 j2 J coeff:   6   0   6 1
    tc->cgs[  285].coeff =  -6.7936622048675700E-1; // j1 j2 J coeff:   6   1   5 -sqrt(78)/13
    tc->cgs[  287].coeff =   7.3379938570534300E-1; // j1 j2 J coeff:   6   1   7 sqrt(91)/13
    tc->cgs[  292].coeff =   5.6096819400507400E-1; // j1 j2 J coeff:   6   2   4 3*sqrt(715)/143
    tc->cgs[  294].coeff =  -5.0452497910951300E-1; // j1 j2 J coeff:   6   2   6 -sqrt(770)/55
    tc->cgs[  296].coeff =   6.5633012331389400E-1; // j1 j2 J coeff:   6   2   8 2*sqrt(455)/65
    tc->cgs[  300].coeff =  -4.8280454958526800E-1; // j1 j2 J coeff:   6   3   3 -10*sqrt(429)/429
    tc->cgs[  302].coeff =   4.2365927286816200E-1; // j1 j2 J coeff:   6   3   5 sqrt(273)/39
    tc->cgs[  304].coeff =  -4.5532635512896800E-1; // j1 j2 J coeff:   6   3   7 -6*sqrt(34034)/2431
    tc->cgs[  306].coeff =   6.1651479928510800E-1; // j1 j2 J coeff:   6   3   9 2*sqrt(4641)/221
    tc->cgs[  309].coeff =   4.1812100500354500E-1; // j1 j2 J coeff:   6   4   2 5*sqrt(143)/143
    tc->cgs[  311].coeff =  -3.7397879600338300E-1; // j1 j2 J coeff:   6   4   4 -2*sqrt(715)/143
    tc->cgs[  313].coeff =   3.8695299497594700E-1; // j1 j2 J coeff:   6   4   6 2*sqrt(1309)/187
    tc->cgs[  315].coeff =  -4.3069561387887500E-1; // j1 j2 J coeff:   6   4   8 -6*sqrt(38038)/2717
    tc->cgs[  317].coeff =   5.9167842041038500E-1; // j1 j2 J coeff:   6   4  10 7*sqrt(125970)/4199
    tc->cgs[  319].coeff =  -3.5478743759345000E-1; // j1 j2 J coeff:   6   5   1 -3*sqrt(286)/143
    tc->cgs[  321].coeff =   3.3796318470968700E-1; // j1 j2 J coeff:   6   5   3 7*sqrt(429)/429
    tc->cgs[  323].coeff =  -3.4736673714593700E-1; // j1 j2 J coeff:   6   5   5 -4*sqrt(3315)/663
    tc->cgs[  325].coeff =   3.6931844203655700E-1; // j1 j2 J coeff:   6   5   7 30*sqrt(323323)/46189
    tc->cgs[  327].coeff =  -4.1565419288457300E-1; // j1 j2 J coeff:   6   5   9 -2*sqrt(255255)/2431
    tc->cgs[  329].coeff =   5.7452466451886100E-1; // j1 j2 J coeff:   6   5  11 3*sqrt(646646)/4199
    tc->cgs[  330].coeff =   2.7735009811261500E-1; // j1 j2 J coeff:   6   6   0 sqrt(13)/13
    tc->cgs[  332].coeff =  -3.1289310938737200E-1; // j1 j2 J coeff:   6   6   2 -sqrt(2002)/143
    tc->cgs[  334].coeff =   3.2196435336464700E-1; // j1 j2 J coeff:   6   6   4 6*sqrt(17017)/2431
    tc->cgs[  336].coeff =  -3.3553079968086100E-1; // j1 j2 J coeff:   6   6   6 -20*sqrt(3553)/3553
    tc->cgs[  338].coeff =   3.5891301156572900E-1; // j1 j2 J coeff:   6   6   8 5*sqrt(38038)/2717
    tc->cgs[  340].coeff =  -4.0544662514408500E-1; // j1 j2 J coeff:   6   6  10 -126*sqrt(96577)/96577
    tc->cgs[  342].coeff =   5.6189620668849700E-1; // j1 j2 J coeff:   6   6  12 66*sqrt(676039)/96577
    // clang-format on

    return;
}

struct TestCase_math_gaussian : public TestCaseFunction1D
{

    TestCase_math_gaussian(String name) : TestCaseFunction1D(name) {};
};

template<>
inline void TestCaseContainer<TestCase_math_gaussian>::setup()
{
    TestCase_math_gaussian* tc = nullptr;

    Real width = 0.1;
    testCases.push_back(TestCase_math_gaussian("width = " + std::to_string(width)));
    tc = &(testCases.back());
    tc->parameters["width"] = width;
    // random x-values: seed = 1, xmin = -10, xmax = 10, n = 20
    // function parameters =  {'width': 0.1}
    // clang-format off
    tc->x.push_back(-1.6595599059485195E+00); tc->fx.push_back( 7.5925779519378211E-01); tc->dfx.push_back( 2.5200675903649472E-01);
    tc->x.push_back( 4.4064898688431615E+00); tc->fx.push_back( 1.4345766619192118E-01); tc->dfx.push_back(-1.2642895053651698E-01);
    tc->x.push_back(-9.9977125036531014E+00); tc->fx.push_back( 4.5608086092699085E-05); tc->dfx.push_back( 9.1195306519332956E-05);
    tc->x.push_back(-3.9533485473632046E+00); tc->fx.push_back( 2.0952829822253113E-01); tc->dfx.push_back( 1.6566767868190557E-01);
    tc->x.push_back(-7.0648821836577387E+00); tc->fx.push_back( 6.7971217606534583E-03); tc->dfx.push_back( 9.6041728853985871E-03);
    tc->x.push_back(-8.1532281046240449E+00); tc->fx.push_back( 1.2972445348442211E-03); tc->dfx.push_back( 2.1153461200123700E-03);
    tc->x.push_back(-6.2747957724465824E+00); tc->fx.push_back( 1.9500675000693614E-02); tc->dfx.push_back( 2.4472550610841411E-02);
    tc->x.push_back(-3.0887854591390447E+00); tc->fx.push_back( 3.8517420702387206E-01); tc->dfx.push_back( 2.3794409797814964E-01);
    tc->x.push_back(-2.0646505153866013E+00); tc->fx.push_back( 6.5293468670456734E-01); tc->dfx.push_back( 2.6961638748367484E-01);
    tc->x.push_back( 7.7633468006713890E-01); tc->fx.push_back( 9.4151071182413115E-01); tc->dfx.push_back(-1.4618548344875421E-01);
    tc->x.push_back(-1.6161097119341044E+00); tc->fx.push_back( 7.7014146509935510E-01); tc->dfx.push_back( 2.4892662026204559E-01);
    tc->x.push_back( 3.7043900079351904E+00); tc->fx.push_back( 2.5353572646276795E-01); tc->dfx.push_back(-1.8783904235265347E-01);
    tc->x.push_back(-5.9109550053696509E+00); tc->fx.push_back( 3.0380968355574790E-02); tc->dfx.push_back( 3.5916107393872361E-02);
    tc->x.push_back( 7.5623487278189074E+00); tc->fx.push_back( 3.2832817340272811E-03); tc->dfx.push_back(-4.9658642888784537E-03);
    tc->x.push_back(-9.4522481360414758E+00); tc->fx.push_back( 1.3176381972044134E-04); tc->dfx.push_back( 2.4909286387004935E-04);
    tc->x.push_back( 3.4093502035680441E+00); tc->fx.push_back( 3.1274507378162025E-01); tc->dfx.push_back(-2.1325149619245401E-01);
    tc->x.push_back(-1.6539039526574602E+00); tc->fx.push_back( 7.6068203895714259E-01); tc->dfx.push_back( 2.5161900618935085E-01);
    tc->x.push_back( 1.1737965689150336E+00); tc->fx.push_back( 8.7129049523638158E-01); tc->dfx.push_back(-2.0454355876734903E-01);
    tc->x.push_back(-7.1922612280953242E+00); tc->fx.push_back( 5.6683219506215351E-03); tc->dfx.push_back( 8.1536104387633854E-03);
    tc->x.push_back(-6.0379702183024246E+00); tc->fx.push_back( 2.6102911067737860E-02); tc->dfx.push_back( 3.1521719927599591E-02);
    // clang-format on

    width = 3.0;
    testCases.push_back(TestCase_math_gaussian("width = " + std::to_string(width)));
    tc = &(testCases.back());
    tc->parameters["width"] = width;
    // random x-values: seed = 2, xmin = -0.1, xmax = 0.1, n = 10
    // function parameters =  {'width': 3.0}
    // clang-format off
    tc->x.push_back(-1.2801019571599251E-02); tc->fx.push_back( 9.9950852250843170E-01); tc->dfx.push_back( 7.6768368951664112E-02);
    tc->x.push_back(-9.4814753634421739E-02); tc->fx.push_back( 9.7339091730782701E-01); tc->dfx.push_back( 5.5375092008715232E-01);
    tc->x.push_back( 9.9324955757418287E-03); tc->fx.push_back( 9.9970408038776171E-01); tc->dfx.push_back(-5.9577338133014981E-02);
    tc->x.push_back(-1.2935521476344622E-02); tc->fx.push_back( 9.9949814282486571E-01); tc->dfx.push_back( 7.7574178152465686E-02);
    tc->x.push_back(-1.5926439582502200E-02); tc->fx.push_back( 9.9923933501907480E-01); tc->dfx.push_back( 9.5485949385845814E-02);
    tc->x.push_back(-3.3933035799225172E-02); tc->fx.push_back( 9.9655160665687992E-01); tc->dfx.push_back( 2.0289612806677962E-01);
    tc->x.push_back(-5.9070273192431501E-02); tc->fx.push_back( 9.8958670617772648E-01); tc->dfx.push_back( 3.5073094248910047E-01);
    tc->x.push_back( 2.3854193270132745E-02); tc->fx.push_back( 9.9829438860147068E-01); tc->dfx.push_back(-1.4288104371713092E-01);
    tc->x.push_back(-4.0069065265095377E-02); tc->fx.push_back( 9.9519499119445809E-01); tc->dfx.push_back( 2.3925919832200054E-01);
    tc->x.push_back(-4.6634544979426679E-02); tc->fx.push_back( 9.9349689495368254E-01); tc->dfx.push_back( 2.7798765380782947E-01);
    // clang-format on

    return;
}

struct TestCase_math_sphericalBessel : public TestCaseFunction1D
{
    TestCase_math_sphericalBessel(String name) : TestCaseFunction1D(name) {};
};

template<>
inline void TestCaseContainer<TestCase_math_sphericalBessel>::setup()
{
    TestCase_math_sphericalBessel* tc = nullptr;

    Real nu = 0.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    //tc->dtolerance = 4.0E-08;
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 0}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 9.9998333341666645E-01); tc->dfx.push_back(-3.3333000001190523E-03);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 9.9994742846005780E-01); tc->dfx.push_back(-5.9199691754290833E-03);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 9.9983417701028687E-01); tc->dfx.push_back(-1.0513427678672796E-02);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 9.9947699361247344E-01); tc->dfx.push_back(-1.8668337014296450E-02);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 9.9835081655750690E-01); tc->dfx.push_back(-3.3133424925372648E-02);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 9.9480345598009989E-01); tc->dfx.push_back(-5.8721064792020827E-02);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 9.8366337615787558E-01); tc->dfx.push_back(-1.0359078875950513E-01);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 9.4901284484345705E-01); tc->dfx.push_back(-1.8009635527426940E-01);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 8.4447068545507120E-01); tc->dfx.push_back(-2.9876850997321608E-01);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 8.4147098480789650E-01); tc->dfx.push_back(-3.0116867893975707E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 4.7040002686622402E-02); tc->dfx.push_back(-3.4567749976235596E-01);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back(-1.9178485493262770E-01); tc->dfx.push_back( 9.5089408079170795E-02);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back( 9.3855228388398437E-02); tc->dfx.push_back( 9.4292432279272309E-02);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back( 4.5790942804639620E-02); tc->dfx.push_back(-1.0632457829881295E-01);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back(-9.0908200595518504E-02); tc->dfx.push_back( 8.6667180530517543E-03);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back( 3.2320541294356991E-02); tc->dfx.push_back( 6.7317403088910710E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back( 4.3352522677141118E-02); tc->dfx.push_back(-5.3536029035730827E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back(-5.6552793639973932E-02); tc->dfx.push_back(-1.2859443788918999E-02);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back( 7.8882741927869659E-03); tc->dfx.push_back( 5.1621912841783281E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back( 3.9840744692193147E-02); tc->dfx.push_back(-2.7979524043641023E-02);
    // clang-format on

    nu = 1.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    //tc->dtolerance = 3.0E-08;
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 1}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 3.3333000001190523E-03); tc->dfx.push_back( 3.3332333339285602E-01);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 5.9199691754290833E-03); tc->dfx.push_back( 3.3330179050412700E-01);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 1.0513427678672796E-02); tc->dfx.push_back( 3.3323384048232751E-01);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 1.8668337014296450E-02); tc->dfx.push_back( 3.3301953888103586E-01);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 3.3133424925372648E-02); tc->dfx.push_back( 3.3234391657647377E-01);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 5.8721064792020827E-02); tc->dfx.push_back( 3.3021633459748179E-01);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 1.0359078875950513E-01); tc->dfx.push_back( 3.2354056631385042E-01);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 1.8009635527426940E-01); tc->dfx.push_back( 3.0283193393396635E-01);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 2.9876850997321608E-01); tc->dfx.push_back( 2.4089793803443260E-01);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 3.0116867893975707E-01); tc->dfx.push_back( 2.3913362692838236E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 3.4567749976235596E-01); tc->dfx.push_back(-1.8341166382161492E-01);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back(-9.5089408079170795E-02); tc->dfx.push_back(-1.5374909170095938E-01);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back(-9.4292432279272309E-02); tc->dfx.push_back( 1.2079592332533338E-01);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back( 1.0632457829881295E-01); tc->dfx.push_back( 2.2163258738236741E-02);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back(-8.6667180530517543E-03); tc->dfx.push_back(-8.9332433676781828E-02);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back(-6.7317403088910710E-02); tc->dfx.push_back( 4.2677064846497102E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back( 5.3536029035730827E-02); tc->dfx.push_back( 3.6214385472377007E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back( 1.2859443788918999E-02); tc->dfx.push_back(-5.8065669379846759E-02);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back(-5.1621912841783281E-02); tc->dfx.push_back( 1.3322159755079943E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back( 2.7979524043641023E-02); tc->dfx.push_back( 3.7176028116608285E-02);
    // clang-format on

    nu = 2.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    //tc->dtolerance = 1.0E-08;
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 2}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 6.6666190477513343E-06); tc->dfx.push_back( 1.3333142857936521E-03);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 2.1028473837476034E-05); tc->dfx.push_back( 2.3679556578229542E-03);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 6.6327781651687783E-05); tc->dfx.push_back( 4.2051917367649854E-03);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 2.0918848468112506E-04); tc->dfx.push_back( 7.4663302458620427E-03);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 6.5953341404146411E-04); tc->dfx.push_back( 1.3247744279212173E-02);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 2.0772260938265274E-03); tc->dfx.push_back( 2.3456946667294205E-02);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 6.5208386081612672E-03); tc->dfx.push_back( 4.1260618923378935E-02);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 2.0258521520778203E-02); tc->dfx.push_back( 7.1065804350236378E-02);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 6.0888435675885845E-02); tc->dfx.push_back( 1.1425809883416804E-01);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 6.2035052011373916E-02); tc->dfx.push_back( 1.1506352290563532E-01);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 2.9863749707573356E-01); tc->dfx.push_back( 4.7040002686622395E-02);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back( 1.3473121008512523E-01); tc->dfx.push_back(-1.7592813413024594E-01);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back(-1.3426627079380085E-01); tc->dfx.push_back(-3.6749744796214796E-02);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back(-1.0349416705035301E-02); tc->dfx.push_back( 1.0977438386715804E-01);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back( 8.8544550217413476E-02); tc->dfx.push_back(-3.2815231748709976E-02);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back(-4.7855326622567154E-02); tc->dfx.push_back(-5.6273866176010601E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back(-3.2645316869994952E-02); tc->dfx.push_back( 6.0065092409729820E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back( 5.8822107249783165E-02); tc->dfx.push_back( 2.4790719213102046E-03);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back(-1.6039102536226431E-02); tc->dfx.push_back(-4.9089422967642266E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back(-3.5843669828815858E-02); tc->dfx.push_back( 3.3100048304900431E-02);
    // clang-format on

    nu = 3.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    tc->dtolerance = 20 * std::numeric_limits<Real>::epsilon();
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 3}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 9.5237566138768746E-09); tc->dfx.push_back( 2.8571164022005847E-06);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 5.3353914471854230E-08); tc->dfx.push_back( 9.0121429085375626E-06);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 2.9889117354698134E-07); tc->dfx.push_back( 2.8425593524055858E-05);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 1.6742664275692801E-06); tc->dfx.push_back( 8.9646252182793326E-05);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 9.3761515615140830E-06); tc->dfx.push_back( 2.8259793880452950E-04);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 5.2465415856863946E-05); tc->dfx.push_back( 8.8965090867105970E-04);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 2.9282763403854856E-04); tc->dfx.push_back( 2.7888040731164159E-03);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 1.6212295991189316E-03); tc->dfx.push_back( 8.6246645495286264E-03);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 8.7488419251972038E-03); tc->dfx.push_back( 2.5539579412462800E-02);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 9.0065811171125242E-03); tc->dfx.push_back( 2.6008727542923819E-02);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 1.5205166203053339E-01); tc->dfx.push_back( 9.5901947701689055E-02);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back( 2.2982061816429603E-01); tc->dfx.push_back(-4.9125284446311590E-02);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back(-1.6120468591568626E-03); tc->dfx.push_back(-1.3334510115999693E-01);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back(-1.1207425424605479E-01); tc->dfx.push_back( 3.9461362959877934E-02);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back( 4.8914240879148786E-02); tc->dfx.push_back( 7.0757553534086637E-02);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back( 4.8911508234077195E-02); tc->dfx.push_back(-6.2905021463821681E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back(-6.4417801325729149E-02); tc->dfx.push_back(-1.5467236516467180E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back( 4.4411759904289892E-03); tc->dfx.push_back( 5.7777124663799875E-02);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back( 4.7401096384881589E-02); tc->dfx.push_back(-2.6018280722517294E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back(-3.6513731145740039E-02); tc->dfx.push_back(-2.8888673420103470E-02);
    // clang-format on

    nu = 4.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    tc->tolerance = 30 * std::numeric_limits<Real>::epsilon();
    tc->dtolerance = 200 * std::numeric_limits<Real>::epsilon();
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 4}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 1.0581962482055022E-11); tc->dfx.push_back( 4.2327753728493632E-09);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 1.0528816625862799E-10); tc->dfx.push_back( 2.3712756433389543E-08);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 1.0475716681714829E-09); tc->dfx.push_back( 1.3283885267620602E-07);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 1.0422190955283725E-08); tc->dfx.push_back( 7.4408892281027208E-07);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 1.0366762317025884E-07); tc->dfx.push_back( 4.1666574863064705E-06);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 1.0304801955396851E-06); tc->dfx.push_back( 2.3308763588163295E-05);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 1.0221828167219311E-05); tc->dfx.push_back( 1.2998347549195984E-04);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 1.0072817890854618E-04); tc->dfx.push_back( 7.1770456244821701E-04);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 9.7206278510453498E-04); tc->dfx.push_back( 3.8394339196187446E-03);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 1.0110158084137538E-03); tc->dfx.push_back( 3.9515020750437550E-03);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 5.6149714328844191E-02); tc->dfx.push_back( 5.8468804815793071E-02);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back( 1.8701765534488918E-01); tc->dfx.push_back( 4.2802962819406848E-02);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back( 1.3265422393464399E-01); tc->dfx.push_back(-9.6365063955331134E-02);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back(-7.6819447708562860E-02); tc->dfx.push_back(-6.9396783296853198E-02);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back(-5.7417306021591522E-02); tc->dfx.push_back( 7.5013016343508571E-02);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back( 7.4192292594762566E-02); tc->dfx.push_back( 2.0376011082245439E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back( 2.5836762513213492E-03); tc->dfx.push_back(-6.5279026742836269E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back(-5.6993387724312401E-02); tc->dfx.push_back( 2.1203937085814993E-02);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back( 3.3502664362235438E-02); tc->dfx.push_back( 3.8584605763240687E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back( 2.3672426113569178E-02); tc->dfx.push_back(-4.2150023077542224E-02);
    // clang-format on

    nu = 5.0;
    testCases.push_back(TestCase_math_sphericalBessel("nu = " + std::to_string(nu)));
    tc = &(testCases.back());
    tc->parameters["nu"] = nu;
    tc->tolerance = 300 * std::numeric_limits<Real>::epsilon();
    tc->dtolerance = 2000 * std::numeric_limits<Real>::epsilon();
    // logarithmic grid x-values: xmin = 0.01, xmax = 0.99, n = 9
    // regular grid x-values: xmin = 1.0, xmax = 21, n = 11
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 5}
    // clang-format off
    tc->x.push_back( 1.0000000000000000E-02); tc->fx.push_back( 9.6199726200342818E-15); tc->dfx.push_back( 4.8099789100344535E-12);
    tc->x.push_back( 1.7760467745895393E-02); tc->fx.push_back( 1.6999738233754120E-13); tc->dfx.push_back( 4.7858130708974457E-11);
    tc->x.push_back( 3.1543421455299037E-02); tc->fx.push_back( 3.0040204141031691E-12); tc->dfx.push_back( 4.7616496423603804E-10);
    tc->x.push_back( 5.6022591935202322E-02); tc->fx.push_back( 5.3080996933479789E-11); tc->dfx.push_back( 4.7372347510215110E-09);
    tc->x.push_back( 9.9498743710661974E-02); tc->fx.push_back( 9.3777385963941542E-10); tc->dfx.push_back( 4.7117731704307011E-08);
    tc->x.push_back( 1.7671442284303238E-01); tc->fx.push_back( 1.6558226797225274E-08); tc->dfx.push_back( 4.6827729673272526E-07);
    tc->x.push_back( 3.1385308071381968E-01); tc->fx.push_back( 2.9185134531099393E-07); tc->dfx.push_back( 4.6424402960431329E-06);
    tc->x.push_back( 5.5741775169676944E-01); tc->fx.push_back( 5.1154668883558095E-06); tc->dfx.push_back( 4.5665703347431509E-05);
    tc->x.push_back( 9.8999999999999999E-01); tc->fx.push_back( 8.8092484844012306E-05); tc->dfx.push_back( 4.3816893756506645E-04);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.push_back( 9.2561158611258252E-05); tc->dfx.push_back( 4.5564885674620427E-04);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.push_back( 1.6397480955999116E-02); tc->dfx.push_back( 2.3354752416845960E-02);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.push_back( 1.0681116145650460E-01); tc->dfx.push_back( 5.8844261597083669E-02);
    tc->x.push_back( 7.0000000000000000E+00); tc->fx.push_back( 1.7216747763227055E-01); tc->dfx.push_back(-1.4917899750159336E-02);
    tc->x.push_back( 9.0000000000000000E+00); tc->fx.push_back( 3.5254806537491926E-02); tc->dfx.push_back(-1.0032265206689081E-01);
    tc->x.push_back( 1.1000000000000000E+01); tc->fx.push_back(-9.5892036714996395E-02); tc->dfx.push_back(-5.1125587225025818E-03);
    tc->x.push_back( 1.3000000000000000E+01); tc->fx.push_back( 2.4523866392199695E-03); tc->dfx.push_back( 7.3060421838199496E-02);
    tc->x.push_back( 1.5000000000000000E+01); tc->fx.push_back( 6.5968007076521964E-02); tc->dfx.push_back(-2.3803526579287439E-02);
    tc->x.push_back( 1.7000000000000000E+01); tc->fx.push_back(-3.4614145962123789E-02); tc->dfx.push_back(-4.4776630325915771E-02);
    tc->x.push_back( 1.9000000000000000E+01); tc->fx.push_back(-3.1531413265927961E-02); tc->dfx.push_back( 4.3459952762002166E-02);
    tc->x.push_back( 2.1000000000000000E+01); tc->fx.push_back( 4.6659056622983974E-02); tc->dfx.push_back( 1.0341267078430900E-02);
    // clang-format on

    return;
}

struct TestCase_math_sphericalBesselRoots : public TestCase
{
    Real tolerance;
    /// Given order of spherical Bessel function.
    Int order;
    /// First few roots of spherical Bessel function.
    Vec1Real roots;

    TestCase_math_sphericalBesselRoots(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
inline void TestCaseContainer<TestCase_math_sphericalBesselRoots>::setup()
{
    TestCase_math_sphericalBesselRoots* tc = nullptr;

    testCases.push_back(TestCase_math_sphericalBesselRoots("order = 0"));
    tc = &(testCases.back());
    tc->order = 0;
    tc->roots.push_back(3.1415926535613856);
    tc->roots.push_back(6.2831853072158941);
    tc->roots.push_back(9.4247779606841302);
    tc->roots.push_back(12.566370614338638);
    tc->roots.push_back(15.707963267993136);
    tc->roots.push_back(18.849555921461427);
    tc->roots.push_back(21.991148575115997);
    tc->roots.push_back(25.132741228770552);

    testCases.push_back(TestCase_math_sphericalBesselRoots("order = 1"));
    tc = &(testCases.back());
    tc->order = 1;
    tc->roots.push_back(4.4934094578959050);
    tc->roots.push_back(7.7252518369816121);
    tc->roots.push_back(10.904121659416685);
    tc->roots.push_back(14.066193912830171);
    tc->roots.push_back(17.220755271892973);
    tc->roots.push_back(20.371302959229819);
    tc->roots.push_back(23.519452498760142);
    tc->roots.push_back(26.666054258775073);

    testCases.push_back(TestCase_math_sphericalBesselRoots("order = 2"));
    tc = &(testCases.back());
    tc->order = 2;
    tc->roots.push_back(5.7634591969661377);
    tc->roots.push_back(9.0950113304890525);
    tc->roots.push_back(12.322940970491594);
    tc->roots.push_back(15.514603010844397);
    tc->roots.push_back(18.689036355447012);
    tc->roots.push_back(21.853874222654891);
    tc->roots.push_back(25.012803202215668);
    tc->roots.push_back(28.167829707916948);

    testCases.push_back(TestCase_math_sphericalBesselRoots("order = 3"));
    tc = &(testCases.back());
    tc->order = 3;
    tc->roots.push_back(6.9879320005886170);
    tc->roots.push_back(10.417118547391125);
    tc->roots.push_back(13.698023153189533);
    tc->roots.push_back(16.923621285241072);
    tc->roots.push_back(20.121806174423565);
    tc->roots.push_back(23.304246988985742);
    tc->roots.push_back(26.476763664465491);
    tc->roots.push_back(29.642604540381729);

    testCases.push_back(TestCase_math_sphericalBesselRoots("order = 4"));
    tc = &(testCases.back());
    tc->order = 4;
    tc->roots.push_back(8.1825614525936459);
    tc->roots.push_back(11.704907154571238);
    tc->roots.push_back(15.039664707612209);
    tc->roots.push_back(18.301255959551785);
    tc->roots.push_back(21.525417733471873);
    tc->roots.push_back(24.727565547917120);
    tc->roots.push_back(27.915576199349143);
    tc->roots.push_back(31.093933214154273);

    return;
}

struct TestCase_math_modifiedSphericalBessel : public TestCase
{

    /// Tolerance for function values.
    Real tolerance;
    /// Vector with x-values.
    Vec1Real x;
    /// Expected function values for given x-values and all Bessel function orders.
    Vec2Real fx;
    /// Expected function derivatives for given x-values and all Bessel function orders.
    Vec2Real dfx;
    /// Maximum order of Bessel function.
    UInt lMax;

    TestCase_math_modifiedSphericalBessel(String name) :
        TestCase(name),
        //tolerance(1000 * std::numeric_limits<Real>::epsilon()){};
        tolerance(1.0E-10) {};
};

template<>
inline void TestCaseContainer<TestCase_math_modifiedSphericalBessel>::setup()
{
    TestCase_math_modifiedSphericalBessel* tc = nullptr;

    UInt lMax = 5;
    testCases.push_back(TestCase_math_modifiedSphericalBessel("lMax = " + std::to_string(lMax)));
    tc = &(testCases.back());
    tc->lMax = lMax;
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 0}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->x.push_back( 0.0000000000000000E+00); tc->fx.back().push_back( 1.0000000000000000E+00); tc->dfx.back().push_back( 0.0000000000000000E+00);
    tc->x.push_back( 5.0000000000000000E-01); tc->fx.back().push_back( 1.0421906109874948E+00); tc->dfx.back().push_back( 1.7087070843777216E-01);
    tc->x.push_back( 1.0000000000000000E+00); tc->fx.back().push_back( 1.1752011936438014E+00); tc->dfx.back().push_back( 3.6787944117144228E-01);
    tc->x.push_back( 1.5000000000000000E+00); tc->fx.back().push_back( 1.4195196367298784E+00); tc->dfx.back().push_back( 6.2192665234224576E-01);
    tc->x.push_back( 2.0000000000000000E+00); tc->fx.back().push_back( 1.8134302039235091E+00); tc->dfx.back().push_back( 9.7438274358006138E-01);
    tc->x.push_back( 2.5000000000000000E+00); tc->fx.back().push_back( 2.4200817924159153E+00); tc->dfx.back().push_back( 1.4848830748991089E+00);
    tc->x.push_back( 3.0000000000000000E+00); tc->fx.back().push_back( 3.3392916424699681E+00); tc->dfx.back().push_back( 2.2427901177692662E+00);
    tc->x.push_back( 3.5000000000000000E+00); tc->fx.back().push_back( 4.7264649393242850E+00); tc->dfx.back().push_back( 3.3846742090665796E+00);
    tc->x.push_back( 4.0000000000000000E+00); tc->fx.back().push_back( 6.8224792992819383E+00); tc->dfx.back().push_back( 5.1214383841836373E+00);
    tc->x.push_back( 4.5000000000000000E+00); tc->fx.back().push_back( 1.0000669144887061E+01); tc->dfx.back().push_back( 7.7807668896984374E+00);
    tc->x.push_back( 5.0000000000000000E+00); tc->fx.back().push_back( 1.4840642115557754E+01); tc->dfx.back().push_back( 1.1873861281846018E+01);
    tc->x.push_back( 5.5000000000000000E+00); tc->fx.back().push_back( 2.2244349590252909E+01); tc->dfx.back().push_back( 1.8200665441377549E+01);
    tc->x.push_back( 6.0000000000000000E+00); tc->fx.back().push_back( 3.3618859561713208E+01); tc->dfx.back().push_back( 2.8016129426790453E+01);
    // clang-format on
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 1}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->fx.back().push_back( 0.0000000000000000E+00); tc->dfx.back().push_back( 3.3333333333333333E-01);
    tc->fx.back().push_back( 1.7087070843777216E-01); tc->dfx.back().push_back( 3.5870777723640612E-01);
    tc->fx.back().push_back( 3.6787944117144228E-01); tc->dfx.back().push_back( 4.3944231130091682E-01);
    tc->fx.back().push_back( 6.2192665234224576E-01); tc->dfx.back().push_back( 5.9028410027355072E-01);
    tc->fx.back().push_back( 9.7438274358006138E-01); tc->dfx.back().push_back( 8.3904746034344768E-01);
    tc->fx.back().push_back( 1.4848830748991089E+00); tc->dfx.back().push_back( 1.2321753324966283E+00);
    tc->fx.back().push_back( 2.2427901177692662E+00); tc->dfx.back().push_back( 1.8440982306237907E+00);
    tc->fx.back().push_back( 3.3846742090665796E+00); tc->dfx.back().push_back( 2.7923653912862396E+00);
    tc->fx.back().push_back( 5.1214383841836373E+00); tc->dfx.back().push_back( 4.2617601071901197E+00);
    tc->fx.back().push_back( 7.7807668896984374E+00); tc->dfx.back().push_back( 6.5425505272433107E+00);
    tc->fx.back().push_back( 1.1873861281846018E+01); tc->dfx.back().push_back( 1.0091097602819346E+01);
    tc->fx.back().push_back( 1.8200665441377549E+01); tc->dfx.back().push_back( 1.5625925793388348E+01);
    tc->fx.back().push_back( 2.8016129426790453E+01); tc->dfx.back().push_back( 2.4280149752783057E+01);
    // clang-format on
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 2}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->fx.back().push_back( 0.0000000000000000E+00); tc->dfx.back().push_back( 0.0000000000000000E+00);
    tc->fx.back().push_back( 1.6966360360861979E-02); tc->dfx.back().push_back( 6.9072546272600288E-02);
    tc->fx.back().push_back( 7.1562870129474501E-02); tc->dfx.back().push_back( 1.5319083078301876E-01);
    tc->fx.back().push_back( 1.7566633204538637E-01); tc->dfx.back().push_back( 2.7059398825147302E-01);
    tc->fx.back().push_back( 3.5185608855341766E-01); tc->dfx.back().push_back( 4.4659861074993490E-01);
    tc->fx.back().push_back( 6.3822210253698475E-01); tc->dfx.back().push_back( 7.1901655185472713E-01);
    tc->fx.back().push_back( 1.0965015247007011E+00); tc->dfx.back().push_back( 1.1462885930685651E+00);
    tc->fx.back().push_back( 1.8253156172672160E+00); tc->dfx.back().push_back( 1.8201179656946800E+00);
    tc->fx.back().push_back( 2.9814005111442103E+00); tc->dfx.back().push_back( 2.8853880008254795E+00);
    tc->fx.back().push_back( 4.8134912184214400E+00); tc->dfx.back().push_back( 4.5717727440841447E+00);
    tc->fx.back().push_back( 7.7163253464501427E+00); tc->dfx.back().push_back( 7.2440660739759322E+00);
    tc->fx.back().push_back( 1.2316713894956060E+01); tc->dfx.back().push_back( 1.1482457862310607E+01);
    tc->fx.back().push_back( 1.9610794848317980E+01); tc->dfx.back().push_back( 1.8210732002631463E+01);
    // clang-format on
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 3}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->fx.back().push_back( 0.0000000000000000E+00); tc->dfx.back().push_back( 0.0000000000000000E+00);
    tc->fx.back().push_back( 1.2071048291523293E-03); tc->dfx.back().push_back( 7.3095217276433443E-03);
    tc->fx.back().push_back( 1.0065090524069861E-02); tc->dfx.back().push_back( 3.1302508033195058E-02);
    tc->fx.back().push_back( 3.6372212190958188E-02); tc->dfx.back().push_back( 7.8673766202831202E-02);
    tc->fx.back().push_back( 9.4742522196516493E-02); tc->dfx.back().push_back( 1.6237104416038467E-01);
    tc->fx.back().push_back( 2.0843886982513898E-01); tc->dfx.back().push_back( 3.0471991081676236E-01);
    tc->fx.back().push_back( 4.1528757660143101E-01); tc->dfx.back().push_back( 5.4278475589879305E-01);
    tc->fx.back().push_back( 7.7708047011341419E-01); tc->dfx.back().push_back( 9.3722365142331410E-01);
    tc->fx.back().push_back( 1.3946877452533744E+00); tc->dfx.back().push_back( 1.5867127658908360E+00);
    tc->fx.back().push_back( 2.4324433136746157E+00); tc->dfx.back().push_back( 2.6513193840440037E+00);
    tc->fx.back().push_back( 4.1575359353958765E+00); tc->dfx.back().push_back( 4.3902965981334416E+00);
    tc->fx.back().push_back( 7.0036528095993189E+00); tc->dfx.back().push_back( 7.2231482152474635E+00);
    tc->fx.back().push_back( 1.1673800386525460E+01); tc->dfx.back().push_back( 1.1828261257301005E+01);
    // clang-format on
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 4}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->fx.back().push_back( 0.0000000000000000E+00); tc->dfx.back().push_back( 0.0000000000000000E+00);
    tc->fx.back().push_back( 6.6892752729369567E-05); tc->dfx.back().push_back( 5.3817730185863372E-04);
    tc->fx.back().push_back( 1.1072364609854644E-03); tc->dfx.back().push_back( 4.5289082191425391E-03);
    tc->fx.back().push_back( 5.9293418209147315E-03); tc->dfx.back().push_back( 1.6607739454575752E-02);
    tc->fx.back().push_back( 2.0257260865610173E-02); tc->dfx.back().push_back( 4.4099370032491064E-02);
    tc->fx.back().push_back( 5.4593267026595678E-02); tc->dfx.back().push_back( 9.9252335771947614E-02);
    tc->fx.back().push_back( 1.2749717929736218E-01); tc->dfx.back().push_back( 2.0279227777249406E-01);
    tc->fx.back().push_back( 2.7115467704038754E-01); tc->dfx.back().push_back( 3.8971664577000348E-01);
    tc->fx.back().push_back( 5.4069695695080522E-01); tc->dfx.back().push_back( 7.1881654906486780E-01);
    tc->fx.back().push_back( 1.0296905082609265E+00); tc->dfx.back().push_back( 1.2883427489402530E+00);
    tc->fx.back().push_back( 1.8957750368959121E+00); tc->dfx.back().push_back( 2.2617608984999644E+00);
    tc->fx.back().push_back( 3.4029739554660181E+00); tc->dfx.back().push_back( 3.9100401228120298E+00);
    tc->fx.back().push_back( 5.9913610640382720E+00); tc->dfx.back().push_back( 6.6809994998269007E+00);
    // clang-format on
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // using analytic gradient (user-defined function)
    // function parameters =  {'n': 5}
    tc->fx.push_back(Vec1Real());
    tc->dfx.push_back(Vec1Real());
    // clang-format off
    tc->fx.back().push_back( 0.0000000000000000E+00); tc->dfx.back().push_back( 0.0000000000000000E+00);
    tc->fx.back().push_back( 3.0352800236771830E-06); tc->dfx.back().push_back( 3.0469392445243372E-05);
    tc->fx.back().push_back( 9.9962375200682578E-05); tc->dfx.back().push_back( 5.0746220978136890E-04);
    tc->fx.back().push_back( 7.9616126546980855E-04); tc->dfx.back().push_back( 2.7446967590354972E-03);
    tc->fx.back().push_back( 3.5848483012706567E-03); tc->dfx.back().push_back( 9.5027159617982020E-03);
    tc->fx.back().push_back( 1.1903108529394531E-02); tc->dfx.back().push_back( 2.6025806556048803E-02);
    tc->fx.back().push_back( 3.2796038709344513E-02); tc->dfx.back().push_back( 6.1905101878673152E-02);
    tc->fx.back().push_back( 7.9825586295274836E-02); tc->dfx.back().push_back( 1.3431081481991639E-01);
    tc->fx.back().push_back( 1.7811959211406250E-01); tc->dfx.back().push_back( 2.7351756877971145E-01);
    tc->fx.back().push_back( 3.7306229715276246E-01); tc->dfx.back().push_back( 5.3227411205724318E-01);
    tc->fx.back().push_back( 7.4514086898323662E-01); tc->dfx.back().push_back( 1.0016059941160282E+00);
    tc->fx.back().push_back( 1.4351499733821962E+00); tc->dfx.back().push_back( 1.8373558026854404E+00);
    tc->fx.back().push_back( 2.6867587904680557E+00); tc->dfx.back().push_back( 3.3046022735702159E+00);
    // clang-format on

    return;
}

struct TestCase_math_expModSphBessel : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    /// Scalar @f$r@f$ value (argument of @f$h_{nl}@f$ function).
    Real r;
    /// @f$r'@f$-values for which the function should be computed.
    Vec1Real rp;
    /// Expected function values for given grid and for each l.
    Vec2Real f;
    /// Maximum order of Bessel function.
    UInt lMax;
    /// Sigma used in @f$h_{nl}@f$ expression.
    Real sigma;

    TestCase_math_expModSphBessel(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
inline void TestCaseContainer<TestCase_math_expModSphBessel>::setup()
{
    TestCase_math_expModSphBessel* tc = nullptr;

    Real r = 3.0;
    UInt lMax = 5;
    Real sigma = 0.5;
    testCases.push_back(TestCase_math_expModSphBessel("r = " + std::to_string(r)
                                                      + " lMax = " + std::to_string(lMax)
                                                      + " sigma = " + std::to_string(sigma)));
    tc = &(testCases.back());
    tc->r = r;
    tc->lMax = lMax;
    tc->sigma = sigma;
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // function parameters =  {'n': 0, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->rp.push_back( 0.0000000000000000E+00); tc->f.back().push_back( 1.5229979744712629E-08);
    tc->rp.push_back( 5.0000000000000000E-01); tc->f.back().push_back( 3.1055252289418455E-07);
    tc->rp.push_back( 1.0000000000000000E+00); tc->f.back().push_back( 1.3977609495410312E-05);
    tc->rp.push_back( 1.5000000000000000E+00); tc->f.back().push_back( 3.0858323717339732E-04);
    tc->rp.push_back( 2.0000000000000000E+00); tc->f.back().push_back( 2.8194850674294306E-03);
    tc->rp.push_back( 2.5000000000000000E+00); tc->f.back().push_back( 1.0108844328543891E-02);
    tc->rp.push_back( 3.0000000000000000E+00); tc->f.back().push_back( 1.3888888888888890E-02);
    tc->rp.push_back( 3.5000000000000000E+00); tc->f.back().push_back( 7.2206030918170660E-03);
    tc->rp.push_back( 4.0000000000000000E+00); tc->f.back().push_back( 1.4097425337147153E-03);
    tc->rp.push_back( 4.5000000000000000E+00); tc->f.back().push_back( 1.0286107905779910E-04);
    tc->rp.push_back( 5.0000000000000000E+00); tc->f.back().push_back( 2.7955218991875997E-06);
    tc->rp.push_back( 5.5000000000000000E+00); tc->f.back().push_back( 2.8232221000595999E-08);
    tc->rp.push_back( 6.0000000000000000E+00); tc->f.back().push_back( 1.0576374822717104E-10);
    // clang-format on
    // function parameters =  {'n': 1, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 2.5879758530322990E-07);
    tc->f.back().push_back( 1.2812808705181474E-05);
    tc->f.back().push_back( 2.9143972399709757E-04);
    tc->f.back().push_back( 2.7020065229532040E-03);
    tc->f.back().push_back( 9.7718828509257612E-03);
    tc->f.back().push_back( 1.3503086419753084E-02);
    tc->f.back().push_back( 7.0486839705833273E-03);
    tc->f.back().push_back( 1.3803728975956590E-03);
    tc->f.back().push_back( 1.0095624426043253E-04);
    tc->f.back().push_back( 2.7489298675344724E-06);
    tc->f.back().push_back( 2.7804460076344535E-08);
    tc->f.back().push_back( 1.0429480727957141E-10);
    // clang-format on
    // function parameters =  {'n': 2, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 1.8115373024256960E-07);
    tc->f.back().push_back( 1.0774407319114947E-05);
    tc->f.back().push_back( 2.6000994984054764E-04);
    tc->f.back().push_back( 2.4817342520602808E-03);
    tc->f.back().push_back( 9.1316560434513121E-03);
    tc->f.back().push_back( 1.2763631687242797E-02);
    tc->f.back().push_back( 6.7171256653468263E-03);
    tc->f.back().push_back( 1.3234692276149868E-03);
    tc->f.back().push_back( 9.7252398821108391E-05);
    tc->f.back().push_back( 2.6580754058108749E-06);
    tc->f.back().push_back( 2.6968381906216694E-08);
    tc->f.back().push_back( 1.0141813125718883E-10);
    // clang-format on
    // function parameters =  {'n': 3, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 1.0783614343442181E-07);
    tc->f.back().push_back( 8.3234723222169092E-06);
    tc->f.back().push_back( 2.1921473793027880E-04);
    tc->f.back().push_back( 2.1849785537739795E-03);
    tc->f.back().push_back( 8.2499401770172089E-03);
    tc->f.back().push_back( 1.1730359796524911E-02);
    tc->f.back().push_back( 6.2490261532801294E-03);
    tc->f.back().push_back( 1.2425115197190976E-03);
    tc->f.back().push_back( 9.1951392517737254E-05);
    tc->f.back().push_back( 2.5274235837168989E-06);
    tc->f.back().push_back( 2.5761400841025085E-08);
    tc->f.back().push_back( 9.7251881497822231E-11);
    // clang-format on
    // function parameters =  {'n': 4, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 5.5344896235744116E-08);
    tc->f.back().push_back( 5.9190484644884179E-06);
    tc->f.back().push_back( 1.7475977397877265E-04);
    tc->f.back().push_back( 1.8444488405428695E-03);
    tc->f.back().push_back( 7.2066700021472989E-03);
    tc->f.back().push_back( 1.0482728393474064E-02);
    tc->f.back().push_back( 5.6756213064668038E-03);
    tc->f.back().push_back( 1.1422696309892854E-03);
    tc->f.back().push_back( 8.5332773865105449E-05);
    tc->f.back().push_back( 2.3632093210439025E-06);
    tc->f.back().push_back( 2.4236112120047374E-08);
    tc->f.back().push_back( 9.1963087222678368E-11);
    // clang-format on
    // function parameters =  {'n': 5, 'rij': 3.0, 'sigma': 0.5}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 2.4818799080805664E-08);
    tc->f.back().push_back( 3.8841859738505946E-06);
    tc->f.back().push_back( 1.3183485094089246E-04);
    tc->f.back().push_back( 1.4933102385704035E-03);
    tc->f.back().push_back( 6.0879391763730183E-03);
    tc->f.back().push_back( 9.1096776981564077E-03);
    tc->f.back().push_back( 5.0328215876086733E-03);
    tc->f.back().push_back( 1.0283359639086068E-03);
    tc->f.back().push_back( 7.7729263540219636E-05);
    tc->f.back().push_back( 2.1729421855603136E-06);
    tc->f.back().push_back( 2.2456476461018645E-08);
    tc->f.back().push_back( 8.5756495594987404E-11);
    // clang-format on

    r = 1.75;
    lMax = 5;
    sigma = 0.8;
    testCases.push_back(TestCase_math_expModSphBessel("r = " + std::to_string(r)
                                                      + " lMax = " + std::to_string(lMax)
                                                      + " sigma = " + std::to_string(sigma)));
    tc = &(testCases.back());
    tc->r = r;
    tc->lMax = lMax;
    tc->sigma = sigma;
    // regular grid x-values: xmin = 0.0, xmax = 6.0, n = 13
    // function parameters =  {'n': 0, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->rp.push_back( 0.0000000000000000E+00); tc->f.back().push_back( 9.1393755356047282E-02);
    tc->rp.push_back( 5.0000000000000000E-01); tc->f.back().push_back( 1.0088794855764620E-01);
    tc->rp.push_back( 1.0000000000000000E+00); tc->f.back().push_back( 1.1733424641753087E-01);
    tc->rp.push_back( 1.5000000000000000E+00); tc->f.back().push_back( 1.1606357956394764E-01);
    tc->rp.push_back( 2.0000000000000000E+00); tc->f.back().push_back( 8.7069976528707696E-02);
    tc->rp.push_back( 2.5000000000000000E+00); tc->f.back().push_back( 4.7132378031693845E-02);
    tc->rp.push_back( 3.0000000000000000E+00); tc->f.back().push_back( 1.7982331980689742E-02);
    tc->rp.push_back( 3.5000000000000000E+00); tc->f.back().push_back( 4.7748573994402634E-03);
    tc->rp.push_back( 4.0000000000000000E+00); tc->f.back().push_back( 8.7575642656337316E-04);
    tc->rp.push_back( 4.5000000000000000E+00); tc->f.back().push_back( 1.1040770963374001E-04);
    tc->rp.push_back( 5.0000000000000000E+00); tc->f.back().push_back( 9.5359555535901176E-06);
    tc->rp.push_back( 5.5000000000000000E+00); tc->f.back().push_back( 5.6292114073362367E-07);
    tc->rp.push_back( 6.0000000000000000E+00); tc->f.back().push_back( 2.2671943062659105E-08);
    // clang-format on
    // function parameters =  {'n': 1, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 4.1107723299210863E-02);
    tc->f.back().push_back( 7.5417105685843630E-02);
    tc->f.back().push_back( 8.7829746535856776E-02);
    tc->f.back().push_back( 7.1151705458303738E-02);
    tc->f.back().push_back( 4.0237693270669887E-02);
    tc->f.back().push_back( 1.5790202778167783E-02);
    tc->f.back().push_back( 4.2759335707194882E-03);
    tc->f.back().push_back( 7.9568726811692536E-04);
    tc->f.back().push_back( 1.0143489260137146E-04);
    tc->f.back().push_back( 8.8384685188386589E-06);
    tc->f.back().push_back( 5.2549054020701811E-07);
    tc->f.back().push_back( 2.1290034152173468E-08);
    // clang-format on
    // function parameters =  {'n': 2, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 1.0685858575377773E-02);
    tc->f.back().push_back( 3.4590907607919545E-02);
    tc->f.back().push_back( 5.1822393526292358E-02);
    tc->f.back().push_back( 4.8038183820152493E-02);
    tc->f.back().push_back( 2.9473778927765568E-02);
    tc->f.back().push_back( 1.2207629250388388E-02);
    tc->f.back().push_back( 3.4344831209453378E-03);
    tc->f.back().push_back( 6.5751077587987352E-04);
    tc->f.back().push_back( 8.5676916770929419E-05);
    tc->f.back().push_back( 7.5965430328849506E-06);
    tc->f.back().push_back( 4.5809601479102903E-07);
    tc->f.back().push_back( 1.8778908246261663E-08);
    // clang-format on
    // function parameters =  {'n': 3, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 2.0280119378292849E-03);
    tc->f.back().push_back( 1.2165160345647848E-02);
    tc->f.back().push_back( 2.4655781094281324E-02);
    tc->f.back().push_back( 2.7231080251307178E-02);
    tc->f.back().push_back( 1.8679729254932803E-02);
    tc->f.back().push_back( 8.3493620922167628E-03);
    tc->f.back().push_back( 2.4815913687562087E-03);
    tc->f.back().push_back( 4.9511091342898370E-04);
    tc->f.back().push_back( 6.6620145469057301E-05);
    tc->f.back().push_back( 6.0603042096693074E-06);
    tc->f.back().push_back( 3.7318848853623442E-07);
    tc->f.back().push_back( 1.5566938305693714E-08);
    // clang-format on
    // function parameters =  {'n': 4, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 3.0243745369183331E-04);
    tc->f.back().push_back( 3.4480971230610646E-03);
    tc->f.back().push_back( 9.7431937920522214E-03);
    tc->f.back().push_back( 1.3182401098479292E-02);
    tc->f.back().push_back( 1.0345736170714378E-02);
    tc->f.back().push_back( 5.0828402650300777E-03);
    tc->f.back().push_back( 1.6193762912265109E-03);
    tc->f.back().push_back( 3.4063979128532400E-04);
    tc->f.back().push_back( 4.7777456237421283E-05);
    tc->f.back().push_back( 4.4936672775342631E-06);
    tc->f.back().push_back( 2.8439373649052723E-07);
    tc->f.back().push_back( 1.2137014569165677E-08);
    // clang-format on
    // function parameters =  {'n': 5, 'rij': 1.75, 'sigma': 0.8}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 0.0000000000000000E+00);
    tc->f.back().push_back( 3.7109385526473308E-05);
    tc->f.back().push_back( 8.1599495774399317E-04);
    tc->f.back().push_back( 3.2764301448638742E-03);
    tc->f.back().push_back( 5.5366144435241111E-03);
    tc->f.back().push_back( 5.0588286050322730E-03);
    tc->f.back().push_back( 2.7727602014409065E-03);
    tc->f.back().push_back( 9.5871668508850164E-04);
    tc->f.back().push_back( 2.1481302802848834E-04);
    tc->f.back().push_back( 3.1674348906829151E-05);
    tc->f.back().push_back( 3.1021872361153223E-06);
    tc->f.back().push_back( 2.0299545765982516E-07);
    tc->f.back().push_back( 8.9089188848942571E-09);
    // clang-format on

    r = 0.9;
    lMax = 5;
    sigma = 0.3;
    testCases.push_back(TestCase_math_expModSphBessel("r = " + std::to_string(r)
                                                      + " lMax = " + std::to_string(lMax)
                                                      + " sigma = " + std::to_string(sigma)));
    tc = &(testCases.back());
    tc->r = r;
    tc->lMax = lMax;
    tc->sigma = sigma;
    tc->tolerance = 2.0E-10;
    // logarithmic grid x-values: xmin = 1e-06, xmax = 1.0, n = 11
    // function parameters =  {'n': 0, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->rp.push_back( 9.9999999999999995E-07); tc->f.back().push_back( 1.1108996538365732E-02);
    tc->rp.push_back( 3.9810717055349691E-06); tc->f.back().push_back( 1.1108996540198583E-02);
    tc->rp.push_back( 1.5848931924611141E-05); tc->f.back().push_back( 1.1108996569247340E-02);
    tc->rp.push_back( 6.3095734448019293E-05); tc->f.back().push_back( 1.1108997029639101E-02);
    tc->rp.push_back( 2.5118864315095795E-04); tc->f.back().push_back( 1.1109004326356968E-02);
    tc->rp.push_back( 1.0000000000000000E-03); tc->f.back().push_back( 1.1109119971605748E-02);
    tc->rp.push_back( 3.9810717055349691E-03); tc->f.back().push_back( 1.1110952841353348E-02);
    tc->rp.push_back( 1.5848931924611141E-02); tc->f.back().push_back( 1.1140005900505966E-02);
    tc->rp.push_back( 6.3095734448019303E-02); tc->f.back().push_back( 1.1601455042914979E-02);
    tc->rp.push_back( 2.5118864315095768E-01); tc->f.back().push_back( 1.9074020723857915E-02);
    tc->rp.push_back( 1.0000000000000000E+00); tc->f.back().push_back( 4.7297973347849805E-02);
    // clang-format on
    // function parameters =  {'n': 1, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 3.7029988460972255E-08);
    tc->f.back().push_back( 1.4741903932799198E-07);
    tc->f.back().push_back( 5.8688576693967196E-07);
    tc->f.back().push_back( 2.3364343598763937E-06);
    tc->f.back().push_back( 9.3015151657469277E-06);
    tc->f.back().push_back( 3.7030153038370892E-05);
    tc->f.back().push_back( 1.4742942330685918E-04);
    tc->f.back().push_back( 5.8754079868249438E-04);
    tc->f.back().push_back( 2.3776101473225748E-03);
    tc->f.back().push_back( 1.1733179718342509E-02);
    tc->f.back().push_back( 4.2568176208041611E-02);
    // clang-format on
    // function parameters =  {'n': 2, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 7.4059976921732868E-14);
    tc->f.back().push_back( 1.1737715325984797E-12);
    tc->f.back().push_back( 1.8603025122149083E-11);
    tc->f.back().push_back( 2.9483808049834511E-10);
    tc->f.back().push_back( 4.6728691050695480E-09);
    tc->f.back().push_back( 7.4060094476807776E-08);
    tc->f.back().push_back( 1.1738010598802595E-06);
    tc->f.back().push_back( 1.8610437191615199E-05);
    tc->f.back().push_back( 2.9668063246093444E-04);
    tc->f.back().push_back( 5.0608317860628502E-03);
    tc->f.back().push_back( 3.4527520485437328E-02);
    // clang-format on
    // function parameters =  {'n': 3, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 1.0579996703087905E-19);
    tc->f.back().push_back( 6.6755266243895271E-18);
    tc->f.back().push_back( 4.2119725519316568E-16);
    tc->f.back().push_back( 2.6575750163902555E-14);
    tc->f.back().push_back( 1.6768164750984853E-12);
    tc->f.back().push_back( 1.0579996703058225E-10);
    tc->f.back().push_back( 6.6755266196847381E-09);
    tc->f.back().push_back( 4.2119718064730996E-07);
    tc->f.back().push_back( 2.6574573710610120E-05);
    tc->f.back().push_back( 1.6594126033590802E-03);
    tc->f.back().push_back( 2.5304415965322943E-02);
    // clang-format on
    // function parameters =  {'n': 4, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 1.1755551892308015E-25);
    tc->f.back().push_back( 2.9528611292752291E-23);
    tc->f.back().push_back( 7.4172518029953329E-21);
    tc->f.back().push_back( 1.8631294092839395E-18);
    tc->f.back().push_back( 4.6799692038774601E-16);
    tc->f.back().push_back( 1.1755540018012383E-13);
    tc->f.back().push_back( 2.9528138561028741E-11);
    tc->f.back().push_back( 7.4153693512436785E-09);
    tc->f.back().push_back( 1.8556056954200439E-06);
    tc->f.back().push_back( 4.3646339193417148E-04);
    tc->f.back().push_back( 1.6814429309711269E-02);
    // clang-format on
    // function parameters =  {'n': 5, 'rij': 0.9, 'sigma': 0.3}
    tc->f.push_back(Vec1Real());
    // clang-format off
    tc->f.back().push_back( 1.0686865356636176E-31);
    tc->f.back().push_back( 1.0686865356364875E-28);
    tc->f.back().push_back( 1.0686865352065710E-25);
    tc->f.back().push_back( 1.0686865283927625E-22);
    tc->f.back().push_back( 1.0686864204012381E-19);
    tc->f.back().push_back( 1.0686847088513625E-16);
    tc->f.back().push_back( 1.0686575827328219E-13);
    tc->f.back().push_back( 1.0682276926294331E-10);
    tc->f.back().push_back( 1.0614219954042002E-07);
    tc->f.back().push_back( 9.5579749241222112E-05);
    tc->f.back().push_back( 1.0171429586582804E-02);
    // clang-format on

    return;
}

struct TestCase_math_calculateMean : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    Vec1Real data;
    Real result;

    TestCase_math_calculateMean(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
void TestCaseContainer<TestCase_math_calculateMean>::setup()
{
    TestCase_math_calculateMean* tc = nullptr;
    testCases.push_back(TestCase_math_calculateMean( "math_mean_calculateMean" ));
    tc = &(testCases.back());
    tc->data   = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    tc->result = 4.0;
    return;
}

struct TestCase_math_calculateVariance : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    Vec1Real data;
    Real result;
    TestCase_math_calculateVariance(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
void TestCaseContainer<TestCase_math_calculateVariance>::setup()
{
    TestCase_math_calculateVariance* tc = nullptr;
    testCases.push_back(TestCase_math_calculateVariance( "math_calculateVariance" ));
    tc = &(testCases.back());
    tc->data   = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    tc->result = 4.0;
}


struct TestCase_math_calculateRMSE : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    Vec1Real dataA;
    Vec1Real dataB;
    Real result;
    TestCase_math_calculateRMSE(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
void TestCaseContainer<TestCase_math_calculateRMSE>::setup()
{
    TestCase_math_calculateRMSE* tc = nullptr;
    testCases.push_back(TestCase_math_calculateRMSE( "math_calculateRMSE" ));
    tc = &(testCases.back());
    tc->dataA   = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    tc->dataB   = {7.0,6.0,5.0,4.0,3.0,2.0,1.0};
    tc->result = 4.0;
}

struct TestCase_math_calculateMAE : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    Vec1Real dataA;
    Vec1Real dataB;
    Real result;
    TestCase_math_calculateMAE(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
void TestCaseContainer<TestCase_math_calculateMAE>::setup()
{
    TestCase_math_calculateMAE* tc = nullptr;
    testCases.push_back(TestCase_math_calculateMAE( "math_calculateMAE" ));
    tc = &(testCases.back());
    tc->dataA   =  {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    tc->dataB   =  {7.0,6.0,5.0,4.0,3.0,2.0,1.0};
    tc->result  =  3.4285714285714284;
}


struct TestCase_math_pearsonCorr : public TestCase
{
    /// Tolerance for function values.
    Real tolerance;
    Vec1Real dataA;
    Vec1Real dataB;
    Real result;
    TestCase_math_pearsonCorr(String name) : TestCase(name), tolerance(1.0E-14) {}
};

template<>
void TestCaseContainer<TestCase_math_pearsonCorr>::setup()
{
    TestCase_math_pearsonCorr* tc = nullptr;
    testCases.push_back(TestCase_math_pearsonCorr( "math_pearsonCorr" ));
    tc = &(testCases.back());
    tc->dataA   =  {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    tc->dataB   =  {7.0,6.0,5.0,4.0,3.0,2.0,1.0};
    tc->result  =  -1.0;
}

} //namespace vaspml

#endif
