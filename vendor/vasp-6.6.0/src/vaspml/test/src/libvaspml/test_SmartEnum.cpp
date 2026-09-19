#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE SmartEnum

#include "boost_helpers.hpp"

#include "SmartEnum.hpp"
// Using the enum CutoffType from here because otherwise we would need to
// define one for testing purposes only. This would need to be added to
// SmartEnum header to enable smart enum functionality. However, we do not
// want to mention an enum for testing only in the main code.
#include "cutoff.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <sstream>
#include <stdexcept>

using namespace vaspml;
using CT = vaspml::math::CutoffType;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(ConvertEnumToString_CorrectStringRepresentation)
{
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::NONE), "NONE");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::HARD), "HARD");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::BP), "BP");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::MIWA), "MIWA");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::JINNOUCHI), "JINNOUCHI");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::WMC), "WMC");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::POLY2), "POLY2");
    BOOST_REQUIRE_EQUAL(SmartEnum<CT>::toString(CT::BUMP), "BUMP");
}

BOOST_AUTO_TEST_CASE(ConvertEnumToStringDirect_CorrectStringRepresentation)
{
    BOOST_REQUIRE_EQUAL(toString(CT::NONE), "NONE");
    BOOST_REQUIRE_EQUAL(toString(CT::HARD), "HARD");
    BOOST_REQUIRE_EQUAL(toString(CT::BP), "BP");
    BOOST_REQUIRE_EQUAL(toString(CT::MIWA), "MIWA");
    BOOST_REQUIRE_EQUAL(toString(CT::JINNOUCHI), "JINNOUCHI");
    BOOST_REQUIRE_EQUAL(toString(CT::WMC), "WMC");
    BOOST_REQUIRE_EQUAL(toString(CT::POLY2), "POLY2");
    BOOST_REQUIRE_EQUAL(toString(CT::BUMP), "BUMP");
}

BOOST_AUTO_TEST_CASE(ConvertStringToEnum_CorrectEnumReturned)
{
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("NONE") == CT::NONE);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("HARD") == CT::HARD);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("BP") == CT::BP);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("MIWA") == CT::MIWA);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("JINNOUCHI") == CT::JINNOUCHI);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("WMC") == CT::WMC);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("POLY2") == CT::POLY2);
    BOOST_REQUIRE(SmartEnum<CT>::toEnum("BUMP") == CT::BUMP);
}

BOOST_AUTO_TEST_CASE(ConcatEnumToString_CorrectlyConcatenatedString)
{
    std::ostringstream ss;
    ss << CT::NONE << " " << CT::HARD << " " << CT::BP << " " << CT::MIWA << " " << CT::JINNOUCHI
       << " " << CT::WMC << " " << CT::POLY2 << " " << CT::BUMP;
    BOOST_REQUIRE_EQUAL(ss.str(), "NONE HARD BP MIWA JINNOUCHI WMC POLY2 BUMP");
}

BOOST_AUTO_TEST_CASE(ConvertWrongString_ThrowException)
{
    BOOST_REQUIRE_THROW(SmartEnum<CT>::toEnum("BESTCUTOFF"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(GetListOfStrings_CorrectlyOrderedList)
{
    Vec1String result = SmartEnum<CT>::listNames();
    Vec1String expect = {"NONE", "HARD", "BP", "MIWA", "JINNOUCHI", "WMC", "POLY2", "BUMP"};

    REQUIRE_EQUAL_COLLECTIONS(result, expect, "Color names list");
}

BOOST_AUTO_TEST_CASE(GetListOfEnums_CorrectlyOrderedList)
{
    std::vector<CT> result = SmartEnum<CT>::listEnums();
    std::vector<CT> expect =
        {CT::NONE, CT::HARD, CT::BP, CT::MIWA, CT::JINNOUCHI, CT::WMC, CT::POLY2, CT::BUMP};
    // Workaround: Here, result and expect vectors cannot be compared directly because I found no
    // way of making BOOST test aware of the overloaded << operator of SmartEnums, at least for
    // enums in a nested namespace (like vaspml::math). Instead, I compare them here one by one in a
    // test with a given failure message.
    for (Size i = 0; i < result.size(); ++i)
    {
        // Note: For some reason, the extra parenthesis around the == expression are required to
        // compile!
        BOOST_TEST((result[i] == expect[i]),
                   "Wrong enum list order, expected \"" + SmartEnum<CT>::toString(expect[i])
                       + "\", got \"" + SmartEnum<CT>::toString(result[i]) + "\".");
    }
    // This does not work (as of BOOST 1.71):
    //REQUIRE_EQUAL_COLLECTIONS(result, expect, "Cutoff enums list order");
}

BOOST_AUTO_TEST_SUITE_END()
