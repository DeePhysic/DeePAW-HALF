#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Serializer

#include "boost_helpers.hpp"

#include "Record.hpp"
#include "rec.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(EmptyVec1ShRecThroughBuffer_CorrectRecovery)
{
    Buffer buffer;
    Record original;

    original.add("v", "Vec1ShRec");

    rec::toBuffer(original, buffer);

    Record newRecord;
    rec::fromBuffer(buffer, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("v"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("v"), "Vec1ShRec", "Checking type of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.vcget<Vec1ShRec>("v").size(),
                          0,
                          "Checking size of empty vector");
}

BOOST_AUTO_TEST_CASE(Vec1ShRecThroughBuffer_CorrectRecovery)
{
    Buffer buffer;
    Record original;

    original.add("v", "Vec1ShRec");
    original.get<Vec1ShRec>("v").push_back(std::make_shared<Record>());
    original.get<Vec1ShRec>("v").push_back(std::make_shared<Record>());

    rec::toBuffer(original, buffer);

    Record newRecord;
    rec::fromBuffer(buffer, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("v"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("v"), "Vec1ShRec", "Checking type of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.vcget<Vec1ShRec>("v").size(),
                          2,
                          "Checking size of vector");
}

BOOST_AUTO_TEST_SUITE_END()
