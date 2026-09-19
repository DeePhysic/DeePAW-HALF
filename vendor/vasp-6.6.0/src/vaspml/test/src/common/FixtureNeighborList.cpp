#include "FixtureNeighborList.hpp"

#include "types.hpp"

using namespace vaspml;

FixtureNeighborList_CsPbBr3_40::FixtureNeighborList_CsPbBr3_40() :
    Fixture("CsPbBr3_40 structure with neighbor list (cutoff = 5.0, type-sorted)"),
    neighborList(5.0, true, false, "CsPbBr3_40")
{}

FixtureNeighborList_CaO_16::FixtureNeighborList_CaO_16() :
    Fixture("CaO_16 structure with neighbor list (cutoff = 5.0, type-sorted)"),
    neighborList(5.0, true, false, "CaO_16")
{}

FixtureNeighborList_HCN_24::FixtureNeighborList_HCN_24() :
    Fixture("HCN_24 structure with neighbor list (cutoff = 5.0, type-sorted)"),
    neighborList(5.0, true, false, "HCN_24")
{}
