#ifndef FIXTURENEIGHBORLIST_HPP
#define FIXTURENEIGHBORLIST_HPP

#include "Fixture.hpp"
#include "SampleNeighborList.hpp"
#include "SampleStructure.hpp"

namespace vaspml
{

class Structure;

struct FixtureNeighborList_CsPbBr3_40 : public Fixture
{
    FixtureNeighborList_CsPbBr3_40();

    SampleNeighborList neighborList;
    Structure&         structure = neighborList.structure;
};

struct FixtureNeighborList_CaO_16 : public Fixture
{
    FixtureNeighborList_CaO_16();

    SampleNeighborList neighborList;
    Structure&         structure = neighborList.structure;
};

struct FixtureNeighborList_HCN_24 : public Fixture
{
    FixtureNeighborList_HCN_24();

    SampleNeighborList neighborList;
    Structure&         structure = neighborList.structure;
};

} //namespace vaspml

#endif
