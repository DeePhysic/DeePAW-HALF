#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DeepCopy

#include "boost_helpers.hpp"
#include "record_helpers.hpp"

#include "Record.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(CopyAssignment_AllItemsCopied)
{
    Record original = generateSampleRecord("Nested 1");
    // Add empty Vec1ShRec, requires special care when copying.
    original.add("sub4", "Vec1ShRec");
    Record newRecord = original;

    BOOST_REQUIRE_MESSAGE(newRecord.contains("sub4"),
                          "Checking existence of empty vector of sub-records.");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub"), "ShRec", "Checking existence of sub-record.");
    Record& sub = newRecord.dget<ShRec>("sub");
    REQUIRE_EQUAL_MESSAGE(sub.typeOf("subsub"),
                          "Vec1ShRec",
                          "Checking existence of sub-sub-record.");
    Vec1ShRec& subsub = sub.get<Vec1ShRec>("subsub");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub2"), "ShRec", "Checking existence of sub-record.");
    Record& sub2 = newRecord.dget<ShRec>("sub2");
    REQUIRE_EQUAL_MESSAGE(sub2.typeOf("sub2sub"),
                          "Vec1ShRec",
                          "Checking existence of sub-sub-record.");
    Vec1ShRec& sub2sub = sub2.get<Vec1ShRec>("sub2sub");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub3"), "ShRec", "Checking existence of sub-record.");
    Record& sub3 = newRecord.dget<ShRec>("sub3");

    std::vector<Record*> recordList{&newRecord,
                                    &sub,
                                    &*subsub[0],
                                    &*subsub[1],
                                    &sub2,
                                    &*sub2sub[0],
                                    &*sub2sub[1],
                                    &*sub2sub[2],
                                    &sub3};
    Vec1String recordNames{"root",
                           "sub",
                           "subsub0",
                           "subsub1",
                           "sub2",
                           "sub2sub0",
                           "sub2sub1",
                           "sub2sub2",
                           "sub3"};
    Vec1String::const_iterator iname = recordNames.begin();

    for (const Record* const& curRecord : recordList)
    {
        BOOST_TEST_CONTEXT("Record: " + *iname)
        {
            BOOST_REQUIRE_EQUAL(original.template cget<Real>("energy"),
                                curRecord->cget<Real>("energy"));
            BOOST_REQUIRE_EQUAL(original.template cget<Real>("cutoff"),
                                curRecord->cget<Real>("cutoff"));
            BOOST_REQUIRE_EQUAL(original.template cget<Int>("numAtoms"),
                                curRecord->cget<Int>("numAtoms"));
            BOOST_REQUIRE_EQUAL(original.template cget<Int>("numTypes"),
                                curRecord->cget<Int>("numTypes"));
            BOOST_REQUIRE_EQUAL(original.template cget<String>("system"),
                                curRecord->cget<String>("system"));
            BOOST_REQUIRE_EQUAL(original.template cget<String>("subsystem"),
                                curRecord->cget<String>("subsystem"));
            BOOST_REQUIRE_EQUAL(original.template cget<bool>("typeSort"),
                                curRecord->cget<bool>("typeSort"));
            BOOST_REQUIRE_EQUAL(original.template cget<bool>("distSort"),
                                curRecord->cget<bool>("distSort"));
            REQUIRE_EQUAL_COLLECTIONS(original.template cget<Vec1Real>("dist"),
                                      curRecord->cget<Vec1Real>("dist"), "dist");
            REQUIRE_EQUAL_COLLECTIONS(original.template cget<Vec1Int>("atoms"),
                                      curRecord->cget<Vec1Int>("atoms"), "atoms");
            REQUIRE_EQUAL_COLLECTIONS(original.template cget<Vec1String>("types"),
                                      curRecord->cget<Vec1String>("types"), "types");
            REQUIRE_EQUAL_COLLECTIONS_2D(original.template cget<Vec2Real>("desc"),
                                         curRecord->cget<Vec2Real>("desc"), "desc");
            REQUIRE_EQUAL_COLLECTIONS_2D(original.template cget<Vec2Int>("nlm"),
                                         curRecord->cget<Vec2Int>("nlm"), "nlm");
            REQUIRE_EQUAL_COLLECTIONS_2D(original.template cget<Vec2String>("info"),
                                         curRecord->cget<Vec2String>("info"), "info");
            iname++;
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
