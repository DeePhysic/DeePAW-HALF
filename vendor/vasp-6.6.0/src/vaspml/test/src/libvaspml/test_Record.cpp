#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Record

#include "test_Record.hpp"

#include "boost_helpers.hpp"

#include "rec.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseRecord> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(ConvertToByteBuffer_CorrectPackUnpack,
                     bdata::make(container.testCases),
                     testCase)
{
    Buffer buffer;
    rec::toBuffer(testCase.record, buffer);
    
    Record newRecord;
    rec::fromBuffer(buffer, newRecord);

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub"), "ShRec", "Checking existence of sub-record.");
    Record& sub = newRecord.dget<ShRec>("sub");
    REQUIRE_EQUAL_MESSAGE(sub.typeOf("subsub"),
                          "Vec1ShRec",
                          "Checking existence of sub-sub-record.");
    Vec1ShRec& subsub = sub.get<Vec1ShRec>("subsub");
    REQUIRE_EQUAL_MESSAGE(subsub.size(), 2, "Size of subsub");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub2"), "ShRec", "Checking existence of sub-record.");
    Record& sub2 = newRecord.dget<ShRec>("sub2");
    REQUIRE_EQUAL_MESSAGE(sub2.typeOf("sub2sub"),
                          "Vec1ShRec",
                          "Checking existence of sub-sub-record.");
    Vec1ShRec& sub2sub = sub2.get<Vec1ShRec>("sub2sub");
    REQUIRE_EQUAL_MESSAGE(sub2sub.size(), 3, "Size of sub2sub");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("sub3"), "ShRec", "Checking existence of sub-record.");
    Record& sub3 = newRecord.dget<ShRec>("sub3");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("vsub"),
                          "Vec1ShRec",
                          "Checking existence of sub-record vector.");
    Vec1ShRec& vsub = newRecord.get<Vec1ShRec>("vsub");
    REQUIRE_EQUAL_MESSAGE(vsub.size(), 2, "Size of vsub");

    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("vsub2"),
                          "Vec1ShRec",
                          "Checking existence of empty sub-record vector.");
    Vec1ShRec& vsub2 = newRecord.get<Vec1ShRec>("vsub2");
    REQUIRE_EQUAL_MESSAGE(vsub2.size(), 0, "Size of vsub2");

    std::vector<Record*> recordList{&newRecord,
                                    &sub,
                                    &*subsub[0],
                                    &*subsub[1],
                                    &sub2,
                                    &*sub2sub[0],
                                    &*sub2sub[1],
                                    &*sub2sub[2],
                                    &sub3,
                                    &*vsub[0],
                                    &*vsub[1]};
    Vec1String recordNames{"root",
                           "sub",
                           "subsub0",
                           "subsub1",
                           "sub2",
                           "sub2sub0",
                           "sub2sub1",
                           "sub2sub2",
                           "sub3",
                           "vsub0",
                           "vsub1"};
    Vec1String::const_iterator iname = recordNames.begin();

    for (const Record* const& curRecord : recordList)
    {
        BOOST_TEST_CONTEXT("Record: " + *iname)
        {
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<Real>("energy"),
                                curRecord->cget<Real>("energy"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<Real>("cutoff"),
                                curRecord->cget<Real>("cutoff"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<Int>("numAtoms"),
                                curRecord->cget<Int>("numAtoms"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<Int>("numTypes"),
                                curRecord->cget<Int>("numTypes"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<String>("system"),
                                curRecord->cget<String>("system"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<String>("subsystem"),
                                curRecord->cget<String>("subsystem"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<bool>("typeSort"),
                                curRecord->cget<bool>("typeSort"));
            BOOST_REQUIRE_EQUAL(testCase.record.template cget<bool>("distSort"),
                                curRecord->cget<bool>("distSort"));
            REQUIRE_EQUAL_COLLECTIONS(testCase.record.template cget<Vec1Real>("dist"),
                                      curRecord->cget<Vec1Real>("dist"), "dist");
            REQUIRE_EQUAL_COLLECTIONS(testCase.record.template cget<Vec1Int>("atoms"),
                                      curRecord->cget<Vec1Int>("atoms"), "atoms");
            REQUIRE_EQUAL_COLLECTIONS(testCase.record.template cget<Vec1String>("types"),
                                      curRecord->cget<Vec1String>("types"), "types");
            REQUIRE_EQUAL_COLLECTIONS_2D(testCase.record.template cget<Vec2Real>("desc"),
                                         curRecord->cget<Vec2Real>("desc"), "desc");
            REQUIRE_EQUAL_COLLECTIONS_2D(testCase.record.template cget<Vec2Int>("nlm"),
                                         curRecord->cget<Vec2Int>("nlm"), "nlm");
            REQUIRE_EQUAL_COLLECTIONS_2D(testCase.record.template cget<Vec2String>("info"),
                                         curRecord->cget<Vec2String>("info"), "info");
            iname++;
        }
    }

    BOOST_TEST_CONTEXT("Record merge")
    {
        Record recordA;
        Record recordB;

        recordA.put( "vector", Vec1Int{ 1,2,3,4,5 } );
        recordA.put( "int", (Int) 1 );
        recordA.put( "vector2D", Vec2Real{ {1.0,2.0,3.0},
                                          {4.0,5.0,6.0},
                                          {7.0,8.0,9.0}});

        recordB.put( "vector", Vec1Int{ 6,7,8,9,10 } );
        recordB.put( "int", (Int) 2 );
        recordB.put( "vector2D", Vec2Real{ {1.0,2.0,3.0},
                                           {4.0,5.0,6.0},
                                           {7.0,8.0,9.0}} );

        Vec1Int totalInt       = {1,2};
        Vec1Int totalVector   = {1,2,3,4,5,6,7,8,9,10};
        Vec2Real totalVector2D = { {1.0,2.0,3.0},
                                   {4.0,5.0,6.0},
                                   {7.0,8.0,9.0},
                                   {1.0,2.0,3.0},
                                   {4.0,5.0,6.0},
                                   {7.0,8.0,9.0}};
        rec::merge( recordA, recordB );
        
        REQUIRE_EQUAL_COLLECTIONS( recordA.cget<Vec1Int>( "int" ), totalInt, "integer merge" );
        REQUIRE_EQUAL_COLLECTIONS( recordA.cget<Vec1Int>( "vector" ), totalVector, "Vec1Int merge" );
        REQUIRE_EQUAL_COLLECTIONS_2D( recordA.cget<Vec2Real>( "vector2D" ), totalVector2D, "Vec2Real merge" );
    }



};

BOOST_AUTO_TEST_SUITE_END()
