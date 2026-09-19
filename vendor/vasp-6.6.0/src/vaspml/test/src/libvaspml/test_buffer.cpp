#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE buffer

#include "boost_helpers.hpp"
#include "FixtureRecordTypes.hpp"

#include "buffer.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;
using namespace vaspml::io::detail;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerReal_CorrectResult, FixtureRecordTypes)
{
    process(buf, real1, pos, true);
    process(buf, real2, pos, false);
    REQUIRE_EQUAL_MESSAGE(real1, real2, "Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerInt_CorrectResult, FixtureRecordTypes)
{
    process(buf, int1, pos, true);
    process(buf, int2, pos, false);
    REQUIRE_EQUAL_MESSAGE(int1, int2, "Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerString_CorrectResult, FixtureRecordTypes)
{
    process(buf, string1, pos, true);
    process(buf, string2, pos, false);
    REQUIRE_EQUAL_MESSAGE(string1, string2, "String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerBool_CorrectResult, FixtureRecordTypes)
{
    process(buf, bool1, pos, true);
    process(buf, bool2, pos, false);
    REQUIRE_EQUAL_MESSAGE(bool1, bool2, "bool");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec1Real_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1real1, pos, true);
    process(buf, vec1real2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1real1, vec1real2, "Vec1Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec1Int_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1int1, pos, true);
    process(buf, vec1int2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1int1, vec1int2, "Vec1Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec1String_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1string1, pos, true);
    process(buf, vec1string2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1string1, vec1string2, "Vec1String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec2Real_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2real1, pos, true);
    process(buf, vec2real2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2real1, vec2real2, "Vec2Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec2Int_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2int1, pos, true);
    process(buf, vec2int2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2int1, vec2int2, "Vec2Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerVec2String_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2string1, pos, true);
    process(buf, vec2string2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2string1, vec2string2, "Vec2String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferSerializerAll_CorrectResult, FixtureRecordTypes)
{
    process(buf, real1, pos, true);
    process(buf, int1, pos, true);
    process(buf, string1, pos, true);
    process(buf, bool1, pos, true);
    process(buf, vec1real1, pos, true);
    process(buf, vec1int1, pos, true);
    process(buf, vec1string1, pos, true);
    process(buf, vec2real1, pos, true);
    process(buf, vec2int1, pos, true);
    process(buf, vec2string1, pos, true);

    process(buf, real2, pos, false);
    REQUIRE_EQUAL_MESSAGE(real1, real2, "Real");
    process(buf, int2, pos, false);
    REQUIRE_EQUAL_MESSAGE(int1, int2, "Int");
    process(buf, string2, pos, false);
    REQUIRE_EQUAL_MESSAGE(string1, string2, "String");
    process(buf, bool2, pos, false);
    REQUIRE_EQUAL_MESSAGE(bool1, bool2, "bool");
    process(buf, vec1real2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1real1, vec1real2, "Vec1Real");
    process(buf, vec1int2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1int1, vec1int2, "Vec1Int");
    process(buf, vec1string2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS(vec1string1, vec1string2, "Vec1String");
    process(buf, vec2real2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2real1, vec2real2, "Vec2Real");
    process(buf, vec2int2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2int1, vec2int2, "Vec2Int");
    process(buf, vec2string2, pos, false);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2string1, vec2string2, "Vec2String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawReal_CorrectResult, FixtureRecordTypes)
{
    process(buf, real1, pos, true, 0);
    process(buf, real2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(real1, real2, "Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawInt_CorrectResult, FixtureRecordTypes)
{
    process(buf, int1, pos, true, 0);
    process(buf, int2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(int1, int2, "Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawString_CorrectResult, FixtureRecordTypes)
{
    process(buf, string1, pos, true, 10);
    process(buf, string2, pos, false, 10);
    string1.resize(10, ' ');
    REQUIRE_EQUAL_MESSAGE(string1, string2, "String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawBool_CorrectResult, FixtureRecordTypes)
{
    process(buf, bool1, pos, true, 0);
    process(buf, bool2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(bool1, bool2, "bool");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec1Real_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1real1, pos, true, 0);
    process(buf, vec1real2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS(vec1real1, vec1real2, "Vec1Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec1Int_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1int1, pos, true, 0);
    process(buf, vec1int2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS(vec1int1, vec1int2, "Vec1Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec1String_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec1string1, pos, true, 10);
    process(buf, vec1string2, pos, false, 10);
    for (String& s : vec1string1) s.resize(10, ' ');
    REQUIRE_EQUAL_COLLECTIONS(vec1string1, vec1string2, "Vec1String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec2Real_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2real1, pos, true, 0);
    process(buf, vec2real2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2real1, vec2real2, "Vec2Real");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec2Int_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2int1, pos, true, 0);
    process(buf, vec2int2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2int1, vec2int2, "Vec2Int");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawVec2String_CorrectResult, FixtureRecordTypes)
{
    process(buf, vec2string1, pos, true, 50);
    process(buf, vec2string2, pos, false, 50);
    for (Vec1String& v : vec2string1) for (String& s: v) s.resize(50, ' ');
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2string1, vec2string2, "Vec2String");
}

BOOST_FIXTURE_TEST_CASE(ProcessBufferRawAll_CorrectResult, FixtureRecordTypes)
{
    process(buf, real1, pos, true, 0);
    process(buf, int1, pos, true, 0);
    process(buf, string1, pos, true, 10);
    process(buf, bool1, pos, true, 0);
    process(buf, vec1real1, pos, true, 0);
    process(buf, vec1int1, pos, true, 0);
    process(buf, vec1string1, pos, true, 10);
    process(buf, vec2real1, pos, true, 0);
    process(buf, vec2int1, pos, true, 0);
    process(buf, vec2string1, pos, true, 50);

    process(buf, real2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(real1, real2, "Real");
    process(buf, int2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(int1, int2, "Int");
    process(buf, string2, pos, false, 10);
    string1.resize(10, ' ');
    REQUIRE_EQUAL_MESSAGE(string1, string2, "String");
    process(buf, bool2, pos, false, 0);
    REQUIRE_EQUAL_MESSAGE(bool1, bool2, "bool");
    process(buf, vec1real2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS(vec1real1, vec1real2, "Vec1Real");
    process(buf, vec1int2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS(vec1int1, vec1int2, "Vec1Int");
    process(buf, vec1string2, pos, false, 10);
    for (String& s : vec1string1) s.resize(10, ' ');
    REQUIRE_EQUAL_COLLECTIONS(vec1string1, vec1string2, "Vec1String");
    process(buf, vec2real2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2real1, vec2real2, "Vec2Real");
    process(buf, vec2int2, pos, false, 0);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2int1, vec2int2, "Vec2Int");
    for (Vec1String& v : vec2string1) for (String& s: v) s.resize(50, ' ');
    process(buf, vec2string2, pos, false, 50);
    REQUIRE_EQUAL_COLLECTIONS_2D(vec2string1, vec2string2, "Vec2String");
}

BOOST_AUTO_TEST_SUITE_END()
