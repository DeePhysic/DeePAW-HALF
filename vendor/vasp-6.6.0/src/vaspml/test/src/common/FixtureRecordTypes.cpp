#include "FixtureRecordTypes.hpp"

using namespace vaspml;

FixtureRecordTypes::FixtureRecordTypes() :
    Fixture("Record types"),
    real1(3.1415),
    real2(0.0),
    int1(42),
    int2(0),
    string1("abc def"),
    string2(""),
    bool1(true),
    bool2(false),
    vec1real1({1.1, 2.2, 3.3, 4.4, 5.5}),
    vec1real2(Vec1Real(5)),
    vec1int1({4, 3, 2, 1}),
    vec1int2(Vec1Int(4)),
    vec1string1({"FeO", "H2O", "Ar ", " Xe", "B  ", " N "}),
    vec1string2(Vec1String(6)),
    vec2real1({{1.1, -1.1}, {2.2, -2.2, 4.4}, {3.3}}),
    vec2real2({Vec1Real(2), Vec1Real(3), Vec1Real(1)}),
    vec2int1({{3, -3}, {2, -2, 8, 9, 10}, {1, -1, 2}}),
    vec2int2({Vec1Int(2), Vec1Int(5), Vec1Int(3)}),
    vec2string1({{"ML_RCUT1", "5.0123", "from_incar"},
                 {"ENERGY", "1.2345", "ab_initio", "TOTEN"},
                 {"descriptions", "This is a longer description of a test system"}}),
    vec2string2({Vec1String(3), Vec1String(4), Vec1String(2)}),
    buf(),
    pos(0)
{}
