#ifndef FIXTURESTRUCTURE_HPP
#define FIXTURESTRUCTURE_HPP

#include "Fixture.hpp"
#include "SampleStructure.hpp"

namespace vaspml
{

struct FixtureStructure_CsPbBr3_40 : public Fixture
{
    FixtureStructure_CsPbBr3_40();

    SampleStructure structure;
};

struct FixtureStructure_CaO_16 : public Fixture
{
    FixtureStructure_CaO_16();

    SampleStructure structure;
};

struct FixtureStructure_HCN_24 : public Fixture
{
    FixtureStructure_HCN_24();

    SampleStructure structure;
};

} //namespace vaspml

#endif
