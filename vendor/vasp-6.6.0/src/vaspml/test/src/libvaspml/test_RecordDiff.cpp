#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE RecordDiff

#include "test_RecordDiff.hpp"

#include "boost_helpers.hpp"

#include "Item.hpp"
#include "RecordDiff.hpp"
#include "rec.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseRecordDiff> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(CreateDiffRecord_CorrectDifferencesReported,
                     bdata::make(container.testCases),
                     testCase)
{
    const Record result = rec::diff(testCase.lhs, testCase.rhs);

    BOOST_REQUIRE(result.contains("only-in-1"));
    BOOST_REQUIRE(result.contains("only-in-2"));
    BOOST_REQUIRE(result.contains("local-diff"));
    BOOST_REQUIRE(result.contains("_type"));
    BOOST_REQUIRE_EQUAL(result.cget<String>("_type"), "RecordDiffResult");
    BOOST_REQUIRE(rec::checkType(result, "RecordDiffResult", "test suite"));

    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("energy"), "Real");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-2").cget<String>("system"), "String");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("sub4"), "ShRec");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("vsub3"), "Vec1ShRec");

    const Record* localDiff = &result.dcget<ShRec>("local-diff");
    BOOST_REQUIRE(localDiff->contains("KEYS"));
    const Vec1String& localDiffKeys = localDiff->cget<Vec1String>("KEYS");
    for (const auto& k : localDiffKeys)
    {
        BOOST_TEST_INFO("Checking local-diff key contents: " + k);
        BOOST_REQUIRE(localDiff->contains(k + "-1"));
        BOOST_REQUIRE(localDiff->contains(k + "-2"));
    }
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:cutoff-1"), "Real");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:cutoff-2"), "Vec1Real");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:dist-1"), 5);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:dist-2"), 4);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE2:desc;1-1"), 3);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE2:desc;1-2"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("numAtoms-1"), 100);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("numAtoms-2"), 50);

    Vec1Int index1Atoms{2, 3};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX1:atoms"), index1Atoms, "index1Atoms");
    Vec1Int atoms1{6, 8};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("atoms-1"), atoms1, "atoms1");
    Vec1Int atoms2{5, 10};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("atoms-2"), atoms2, "atoms2");

    Vec1Int index1Types{3};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX1:types"), index1Types, "index1Types");
    Vec1String types1{"B"};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1String>("types-1"), types1, "types1");
    Vec1String types2{"Br"};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1String>("types-2"), types2, "types2");

    Vec1Int index1Nlm{2, 3};
    Vec1Int index2Nlm{1, 3};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX1:nlm"), index1Nlm, "index1Nlm");
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX2:nlm"), index2Nlm, "index2Nlm");
    Vec1Int nlm1{0, 2};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("nlm-1"), nlm1, "nlm1");
    Vec1Int nlm2{4, 8};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("nlm-2"), nlm2, "nlm2");

    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:sub3-1"), "Vec1ShRec");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:sub3-2"), "ShRec");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:vsub-1"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:vsub-2"), 1);

    localDiff = &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub2").dcget<ShRec>("local-diff");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:sub2sub-1"), 3);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:sub2sub-2"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("system-1"), "Test system");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("system-2"), "other_system");
    const Record* onlyIn =
        &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub2").dcget<ShRec>("only-in-1");
    BOOST_REQUIRE_MESSAGE(onlyIn->contains("system new"), "system new");
    BOOST_REQUIRE_EQUAL(onlyIn->cget<String>("system new"), "String");

    onlyIn = &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub5").dcget<ShRec>("only-in-2");
    BOOST_REQUIRE_MESSAGE(onlyIn->contains("date"), "sub5/date");
    BOOST_REQUIRE_EQUAL(onlyIn->cget<String>("date"), "String");

    localDiff = &result.dcget<ShRec>("common-sub")
                     .dcget<ShRec>("sub")
                     .dcget<ShRec>("common-sub")
                     .vcget<Vec1ShRec>("subsub")[0]
                     ->dcget<ShRec>("local-diff");
    index1Atoms = {1};
    atoms1 = {4};
    atoms2 = {14};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX1:atoms"), index1Atoms, "index1Atoms");
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("atoms-1"), atoms1, "atoms1");
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("atoms-2"), atoms2, "atoms2");

    localDiff = &result.dcget<ShRec>("common-sub")
                     .dcget<ShRec>("sub")
                     .dcget<ShRec>("common-sub")
                     .vcget<Vec1ShRec>("subsub")[1]
                     ->dcget<ShRec>("local-diff");
    Vec1Int index1Dist{2, 4};
    Vec1Real dist1{3.0, 5.0};
    Vec1Real dist2{-0.2, -0.4};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Int>("INDEX1:dist"), index1Dist, "index1Dist");
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Real>("dist-1"), dist1, "dist1");
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1Real>("dist-2"), dist2, "dist2");

    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("distSort-1"), true, "distSort-1");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("distSort-2"), false, "distSort-2");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("typeSort-1"), false, "typeSort-1");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("typeSort-2"), true, "typeSort-2");
}

BOOST_DATA_TEST_CASE(CreateDiffRecordWithoutVectorValues_CorrectDifferencesReported,
                     bdata::make(container.testCases),
                     testCase)
{
    const Record result = rec::diff(testCase.lhs, testCase.rhs, false);

    BOOST_REQUIRE(result.contains("only-in-1"));
    BOOST_REQUIRE(result.contains("only-in-2"));
    BOOST_REQUIRE(result.contains("local-diff"));
    BOOST_REQUIRE(result.contains("_type"));
    BOOST_REQUIRE_EQUAL(result.cget<String>("_type"), "RecordDiffResult");
    BOOST_REQUIRE(rec::checkType(result, "RecordDiffResult", "test suite"));

    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("energy"), "Real");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-2").cget<String>("system"), "String");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("sub4"), "ShRec");
    BOOST_REQUIRE_EQUAL(result.dcget<ShRec>("only-in-1").cget<String>("vsub3"), "Vec1ShRec");

    const Record* localDiff = &result.dcget<ShRec>("local-diff");
    BOOST_REQUIRE(localDiff->contains("KEYS"));
    const Vec1String& localDiffKeys = localDiff->cget<Vec1String>("KEYS");
    for (const auto& k : localDiffKeys)
    {
        if (k == "atoms" || k == "nlm" || k == "types") continue;
        BOOST_TEST_INFO("Checking local-diff key contents: " + k);
        BOOST_REQUIRE(localDiff->contains(k + "-1"));
        BOOST_REQUIRE(localDiff->contains(k + "-2"));
    }
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:cutoff-1"), "Real");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:cutoff-2"), "Vec1Real");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:dist-1"), 5);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:dist-2"), 4);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE2:desc;1-1"), 3);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE2:desc;1-2"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("numAtoms-1"), 100);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("numAtoms-2"), 50);

    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:sub3-1"), "Vec1ShRec");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("TYPE:sub3-2"), "ShRec");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:vsub-1"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:vsub-2"), 1);

    localDiff = &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub2").dcget<ShRec>("local-diff");
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:sub2sub-1"), 3);
    BOOST_REQUIRE_EQUAL(localDiff->cget<Int>("SIZE1:sub2sub-2"), 2);
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("system-1"), "Test system");
    BOOST_REQUIRE_EQUAL(localDiff->cget<String>("system-2"), "other_system");
    const Record* onlyIn =
        &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub2").dcget<ShRec>("only-in-1");
    BOOST_REQUIRE_MESSAGE(onlyIn->contains("system new"), "system new");
    BOOST_REQUIRE_EQUAL(onlyIn->cget<String>("system new"), "String");

    onlyIn = &result.dcget<ShRec>("common-sub").dcget<ShRec>("sub5").dcget<ShRec>("only-in-2");
    BOOST_REQUIRE_MESSAGE(onlyIn->contains("date"), "sub5/date");
    BOOST_REQUIRE_EQUAL(onlyIn->cget<String>("date"), "String");

    localDiff = &result.dcget<ShRec>("common-sub")
                     .dcget<ShRec>("sub")
                     .dcget<ShRec>("common-sub")
                     .vcget<Vec1ShRec>("subsub")[0]
                     ->dcget<ShRec>("local-diff");
    Vec1String expectedKeys{"atoms"};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1String>("KEYS"), expectedKeys, "keys");

    localDiff = &result.dcget<ShRec>("common-sub")
                     .dcget<ShRec>("sub")
                     .dcget<ShRec>("common-sub")
                     .vcget<Vec1ShRec>("subsub")[1]
                     ->dcget<ShRec>("local-diff");
    expectedKeys = {"dist", "distSort", "typeSort"};
    REQUIRE_EQUAL_COLLECTIONS(localDiff->cget<Vec1String>("KEYS"), expectedKeys, "keys");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("distSort-1"), true, "distSort-1");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("distSort-2"), false, "distSort-2");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("typeSort-1"), false, "typeSort-1");
    REQUIRE_EQUAL_MESSAGE(localDiff->cget<bool>("typeSort-2"), true, "typeSort-2");
}

BOOST_AUTO_TEST_CASE(CheckEmptyDiff_CorrectResult)
{
    namespace rd = rec::detail;

    Record r1;
    Record r2;
    Record diff = rec::diff(r1, r2);
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(rd::emptyRecordDiffRecursive(diff));

    r1.put("number", 3.1415);
    diff = rec::diff(r1, r2);
    BOOST_REQUIRE(!rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(!rd::emptyRecordDiffRecursive(diff));

    r1.erase("number");
    diff = rec::diff(r1, r2);
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(rd::emptyRecordDiffRecursive(diff));

    r1.add("sub", ItemIndex::SHREC);
    r2.add("sub", ItemIndex::SHREC);
    r2.dget<ShRec>("sub").put("mystring", String("This is a string"));
    diff = rec::diff(r1, r2);
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(!rd::emptyRecordDiffRecursive(diff));

    r2.dget<ShRec>("sub").erase("mystring");
    r1.add("vsub", ItemIndex::VEC1SHREC);
    r2.add("vsub", ItemIndex::VEC1SHREC);
    Vec1ShRec& v1 = r1.get<Vec1ShRec>("vsub");
    Vec1ShRec& v2 = r2.get<Vec1ShRec>("vsub");
    v1.push_back(std::make_shared<Record>());
    v1.push_back(std::make_shared<Record>());
    v2.push_back(std::make_shared<Record>());
    v2.push_back(std::make_shared<Record>());
    diff = rec::diff(r1, r2);
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->cget<ShRec>("sub")));
    BOOST_REQUIRE(
        rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->vcget<Vec1ShRec>("vsub")[0]));
    BOOST_REQUIRE(
        rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->vcget<Vec1ShRec>("vsub")[1]));
    BOOST_REQUIRE(rd::emptyRecordDiffRecursive(diff));

    v1[1]->put("item", Int(2));
    v2[1]->put("item", String("StringItem"));
    diff = rec::diff(r1, r2);
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(diff));
    BOOST_REQUIRE(rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->cget<ShRec>("sub")));
    BOOST_REQUIRE(
        rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->vcget<Vec1ShRec>("vsub")[0]));
    BOOST_REQUIRE(
        !rd::emptyRecordDiffLocal(*diff.cget<ShRec>("common-sub")->vcget<Vec1ShRec>("vsub")[1]));
    BOOST_REQUIRE(!rd::emptyRecordDiffRecursive(diff));
}

BOOST_AUTO_TEST_SUITE_END()
