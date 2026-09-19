#include "FixtureStructure.hpp"

#include "types.hpp"

using namespace vaspml;

FixtureStructure_CsPbBr3_40::FixtureStructure_CsPbBr3_40() :
    Fixture("CsPbBr3_40 structure"),
    structure(SampleStructure("CsPbBr3_40"))
{}

FixtureStructure_CaO_16::FixtureStructure_CaO_16() :
    Fixture("CaO_16 structure"),
    structure(SampleStructure("CaO_16"))
{}

FixtureStructure_HCN_24::FixtureStructure_HCN_24() :
    Fixture("HCN_24 structure"),
    structure(SampleStructure("HCN_24"))
{}
