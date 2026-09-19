#ifndef FIXTURERECORDTYPES_HPP
#define FIXTURERECORDTYPES_HPP

#include "Fixture.hpp"

#include "types.hpp"

namespace vaspml
{

struct FixtureRecordTypes : public Fixture
{
    FixtureRecordTypes();

    Real       real1;
    Real       real2;
    Int        int1;
    Int        int2;
    String     string1;
    String     string2;
    bool       bool1;
    bool       bool2;
    Vec1Real   vec1real1;
    Vec1Real   vec1real2;
    Vec1Int    vec1int1;
    Vec1Int    vec1int2;
    Vec1String vec1string1;
    Vec1String vec1string2;
    Vec2Real   vec2real1;
    Vec2Real   vec2real2;
    Vec2Int    vec2int1;
    Vec2Int    vec2int2;
    Vec2String vec2string1;
    Vec2String vec2string2;

    // For testing buffer functions.
    Buffer buf;
    Size   pos;
};

} // namespace vaspml

#endif
