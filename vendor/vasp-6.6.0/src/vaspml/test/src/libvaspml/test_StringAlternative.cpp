#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE StringAlternative

#include "test_StringAlternative.hpp"

#include "boost_helpers.hpp"

#include "StringAlternative.hpp"
#include "Record.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseStringAlternative> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(Setup_CorrectDataContents,
                     bdata::make(container.testCases),
                     testCase)
{
    ShRec data = std::make_shared<Record>();
    StringAlternative sa(data);
    sa.add(testCase.alternatives);

    for (const auto& alt : testCase.alternatives)
    {
        for (const auto& s : alt)
        {
            REQUIRE_EQUAL_MESSAGE(alt[0], sa(s), "Compare " + alt[0] + " to " + s + ".");
        }
    }

}

BOOST_AUTO_TEST_SUITE_END()
