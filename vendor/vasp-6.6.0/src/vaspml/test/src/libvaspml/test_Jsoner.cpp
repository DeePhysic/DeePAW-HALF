#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Jsoner

#include "test_Record.hpp" // IWYU pragma: associated

#include "boost_helpers.hpp"

#include "rec.hpp"

#include "boost_unit_test_case.hpp"

#include <limits>
#include <regex>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

const Real tol = 10.0 * std::numeric_limits<Real>::epsilon();
TestCaseContainer<TestCaseRecord> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(EmptyRecord_CorrectJsonAndDejson)
{
    Record original;
    const char* jsonExpected =
R"raw({
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.empty(), "Checking for empty record");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Real_CorrectJsonAndDejson)
{
    Record original;
    Real item = 1.234567890123;
    original.put("item", item);
//    const char* jsonExpected =
//R"raw({
//  "item" :   1.2345678901229999E+00
//}
//)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    BOOST_REQUIRE_SMALL(newRecord.cget<Real>("item") - item, tol);
    // Can only check for part of the string because string representing number may depend on
    // compiler.
    BOOST_REQUIRE_MESSAGE(json.find("  \"item\" :   1.23456789012") != json.npos,
                          "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Int_CorrectJsonAndDejson)
{
    Record original;
    Int item = 123456789;
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : 123456789
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.cget<Int>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(String_CorrectJsonAndDejson)
{
    Record original;
    String item = "my string";
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : "my string"
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.cget<String>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(StringSpecial_CorrectJsonAndDejson)
{
    Record original;
    String item = "my \"special\" string: backslash \\ slash / escaped slash newline \n";
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : "my \"special\" string: backslash \\ slash / escaped slash \/ newline \n"
}
)raw";
    String json = rec::toJson(original, 2);
    // Enforce escaped forward slash in JSON, should be converted back to Record string the same way
    // as non-escaped slash.
    json = std::regex_replace(json, std::regex("escaped slash"), "escaped slash \\/");
    item = std::regex_replace(item, std::regex("escaped slash"), "escaped slash /");
     
    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.cget<String>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(BoolTrue_CorrectJsonAndDejson)
{
    Record original;
    bool item = true;
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : true
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.cget<bool>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(BoolFalse_CorrectJsonAndDejson)
{
    Record original;
    bool item = false;
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : false
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.cget<bool>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec1Real_CorrectJsonAndDejson)
{
    Record original;
    Vec1Real item{1.111, 2.222, 3.333, 4.444};
    original.put("item", item);
//    const char* jsonExpected =
//R"raw({
//  "item" : [
//      1.1110000000000000E+00,
//      2.2220000000000000E+00,
//      3.3330000000000000E+00,
//      4.4440000000000000E+00
//  ]
//}
//)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_CLOSE_COLLECTIONS(newRecord.cget<Vec1Real>("item"), item, tol, "Checking item value");
    // Can only check for part of the string because string representing number may depend on
    // compiler.
    BOOST_REQUIRE_MESSAGE(json.find("  \"item\" : [") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("    1.11") != json.npos, "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("    2.22") != json.npos, "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("    3.33") != json.npos, "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("    4.44") != json.npos, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec1Int_CorrectJsonAndDejson)
{
    Record original;
    Vec1Int item{1, 2, 3, 4};
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : [
    1,
    2,
    3,
    4
  ]
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_COLLECTIONS(newRecord.cget<Vec1Int>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec1String_CorrectJsonAndDejson)
{
    Record original;
    Vec1String item{"abc", "def", "hijklm", "nop"};
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : [
    "abc",
    "def",
    "hijklm",
    "nop"
  ]
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_COLLECTIONS(newRecord.cget<Vec1String>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec2Real_CorrectJsonAndDejson)
{
    Record original;
    Vec2Real item{{1.111, 2.222, 3.333}, {4.444, 5.555}, {6.666, 7.777, 8.888, 9.999}};
    original.put("item", item);
//    const char* jsonExpected =
//R"raw({
//  "item" : [
//    [
//        1.1110000000000000E+00,
//        2.2220000000000000E+00,
//        3.3330000000000002E+00
//    ],
//    [
//        4.4440000000000000E+00,
//        5.5549999999999997E+00
//    ],
//    [
//        6.6660000000000004E+00,
//        7.7770000000000001E+00,
//        8.8879999999999999E+00,
//        9.9990000000000006E+00
//    ]
//  ]
//}
//)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_CLOSE_COLLECTIONS_2D(newRecord.cget<Vec2Real>("item"),
                                 item,
                                 tol,
                                 "Checking item value");
    // Can only check for part of the string because string representing number may depend on
    // compiler.
    BOOST_REQUIRE_MESSAGE(json.find("  \"item\" : [") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        1.11") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        2.22") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        3.33") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        4.44") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        5.55") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        6.66") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        7.77") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        8.88") != json.npos,
                          "Checking Json string represenation");
    BOOST_REQUIRE_MESSAGE(json.find("        9.99") != json.npos,
                          "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec2Int_CorrectJsonAndDejson)
{
    Record original;
    Vec2Int item{{1, 2, 3}, {4, 5}, {6, 7, 8, 9}};
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : [
    [
      1,
      2,
      3
    ],
    [
      4,
      5
    ],
    [
      6,
      7,
      8,
      9
    ]
  ]
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_COLLECTIONS_2D(newRecord.cget<Vec2Int>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec2String_CorrectJsonAndDejson)
{
    Record original;
    Vec2String item{{"abc"}, {"def", "hijklm", "nop"}, {"qrst", "uvw"}};
    original.put("item", item);
    const char* jsonExpected =
R"raw({
  "item" : [
    [
      "abc"
    ],
    [
      "def",
      "hijklm",
      "nop"
    ],
    [
      "qrst",
      "uvw"
    ]
  ]
}
)raw";
    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_COLLECTIONS_2D(newRecord.cget<Vec2String>("item"), item, "Checking item value");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Record_CorrectJsonAndDejson)
{
    Record original;
    original.add("item", "ShRec");
    const char* jsonExpected =
R"raw({
  "item" : {
  }
}
)raw";

    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("item"), "ShRec", "Checking type of item");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(EmptyVec1ShRec_CorrectJsonAndDejson)
{
    Record original;
    original.add("item", "Vec1ShRec");
    const char* jsonExpected =
R"raw({
  "item" : [
  ]
}
)raw";

    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("item"), "Vec1ShRec", "Checking type of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.vcget<Vec1ShRec>("item").size(),
                          0,
                          "Checking size of vector");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_AUTO_TEST_CASE(Vec1ShRec_CorrectJsonAndDejson)
{
    Record original;
    original.add("item", "Vec1ShRec");
    original.get<Vec1ShRec>("item").push_back(std::make_shared<Record>());
    original.get<Vec1ShRec>("item").push_back(std::make_shared<Record>());
    const char* jsonExpected =
R"raw({
  "item" : [
    {
    },
    {
    }
  ]
}
)raw";

    String json = rec::toJson(original, 2);

    Record newRecord;
    rec::fromJson(json, newRecord);

    BOOST_REQUIRE_MESSAGE(newRecord.contains("item"), "Checking existence of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("item"), "Vec1ShRec", "Checking type of item");
    REQUIRE_EQUAL_MESSAGE(newRecord.vcget<Vec1ShRec>("item").size(),
                          2,
                          "Checking size of vector");
    REQUIRE_EQUAL_MESSAGE(json, jsonExpected, "Checking Json string represenation");
}

BOOST_DATA_TEST_CASE(FullRecord_CorrectJsonAndDejson,
                     bdata::make(container.testCases),
                     testCase)
{
    String json = rec::toJson(testCase.record);
    
    Record newRecord;
    rec::fromJson(json, newRecord);

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

    //REQUIRE_EQUAL_MESSAGE(newRecord.typeOf("vsub2"),
    //                      "Vec1ShRec",
    //                      "Checking existence of empty sub-record vector.");
    //Vec1ShRec& vsub2 = newRecord.get<Vec1ShRec>("vsub2");
    //REQUIRE_EQUAL_MESSAGE(vsub2.size(), 0, "Size of vsub2");

    std::vector<Record*> recordList{&newRecord,
                                    &sub,
                                    &*subsub[0],
                                    &*subsub[1],
                                    &sub2,
                                    &*sub2sub[0],
                                    &*sub2sub[1],
                                    &*sub2sub[2],
                                    &sub3,
                                    &*vsub[0]};//,
                                    //&*vsub[1]};
    Vec1String recordNames{"root",
                           "sub",
                           "subsub0",
                           "subsub1",
                           "sub2",
                           "sub2sub0",
                           "sub2sub1",
                           "sub2sub2",
                           "sub3",
                           "vsub0"};//,
                           //"vsub1"};
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
}

BOOST_AUTO_TEST_SUITE_END()
