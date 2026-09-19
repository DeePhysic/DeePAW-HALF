#ifndef TEST_RECORDDIFF_HPP
#define TEST_RECORDDIFF_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"
#include "record_helpers.hpp"

#include "Record.hpp"
#include "types.hpp"

namespace vaspml
{

struct TestCaseRecordDiff : public TestCase
{
    Record lhs;
    Record rhs;

    TestCaseRecordDiff(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCaseRecordDiff>::setup()
{
    String              testName = "";
    TestCaseRecordDiff* tc = nullptr;

    testName = "Nested 1";
    testCases.push_back(TestCaseRecordDiff(testName));
    tc = &(testCases.back());
    tc->lhs = generateSampleRecord(testName);

    tc->rhs = tc->lhs;

    /*============================================================================================+
     | Artificial differences
     +============================================================================================*/
    // Highest level changes to simple types.
    tc->lhs.erase("system");
    tc->rhs.erase("energy");
    tc->rhs.erase("cutoff");
    tc->rhs.put("cutoff", Vec1Real({5.0, 6.0}), true);
    tc->rhs.get<Vec1Real>("dist").pop_back();
    tc->rhs.get<Vec2Real>("desc")[1].pop_back();
    tc->rhs.get<Int>("numAtoms") = 50;
    tc->rhs.get<Vec1Int>("atoms")[2] = 5;
    tc->rhs.get<Vec1Int>("atoms")[3] = 10;
    tc->rhs.get<Vec1String>("types")[3] = "Br";
    tc->rhs.get<Vec2Int>("nlm")[2][1] = 4;
    tc->rhs.get<Vec2Int>("nlm")[3][3] = 8;
    // Highest level changes to sub-records.
    tc->lhs.add("sub4", "ShRec");
    tc->lhs.add("vsub3", "Vec1ShRec");
    tc->lhs.erase("sub3");
    tc->lhs.add("sub3", "Vec1ShRec");
    tc->rhs.get<Vec1ShRec>("vsub").pop_back();
    tc->lhs.add("sub5", "ShRec");
    tc->rhs.add("sub5", "ShRec");
    tc->rhs.get<ShRec>("sub5")->put("date", String("1970-01-01T00:00:00"));
    // Sub-record changes.
    tc->rhs.get<ShRec>("sub2")->get<Vec1ShRec>("sub2sub").pop_back();
    tc->rhs.get<ShRec>("sub2")->get<String>("system") = "other_system";
    tc->lhs.get<ShRec>("sub2")->put("system new", String("new_system"));
    tc->rhs.get<ShRec>("sub")->get<Vec1ShRec>("subsub")[0]->get<Vec1Int>("atoms")[1] = 14;
    tc->rhs.get<ShRec>("sub")->get<Vec1ShRec>("subsub")[1]->get<Vec1Real>("dist")[2] = -0.2;
    tc->rhs.get<ShRec>("sub")->get<Vec1ShRec>("subsub")[1]->get<Vec1Real>("dist")[4] = -0.4;
    tc->rhs.get<ShRec>("sub")->get<Vec1ShRec>("subsub")[1]->get<bool>("distSort") = false;
    tc->rhs.get<ShRec>("sub")->get<Vec1ShRec>("subsub")[1]->get<bool>("typeSort") = true;

    return;
}

} //namespace vaspml

#endif
