#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE TypeMap

#include "TypeMap.hpp"
#include "types.hpp"

#include "boost_helpers.hpp"

#include "boost_unit_test.hpp"

#include <stdexcept>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(RandomOrderMapping_CorrectMapping)
{
    TypeMap tm({"H", "O", "Na", "Cl", "S"}, {"Cl", "O", "S"});

    BOOST_REQUIRE_EQUAL(tm.countForceFieldTypes(), 5);
    BOOST_REQUIRE_EQUAL(tm.countStructureTypes(), 3);

    Vec1Int       toType;
    const Vec1Int toTypeExpected = {3, 1, 4};
    for (Size i = 0; i < tm.countStructureTypes(); ++i) toType.push_back(tm.toType(i));
    REQUIRE_EQUAL_COLLECTIONS(toType, toTypeExpected, "toType");

    Vec1Int       toSubType;
    const Vec1Int toSubTypeExpected = {-1, 1, -1, 0, 2};
    for (Size i = 0; i < tm.countForceFieldTypes(); ++i) { toSubType.push_back(tm.toSubType(i)); }
    REQUIRE_EQUAL_COLLECTIONS(toSubType, toSubTypeExpected, "toSubType");

    Vec1String       nameType;
    const Vec1String nameTypeExpected = {"H", "O", "Na", "Cl", "S"};
    for (Size i = 0; i < tm.countForceFieldTypes(); ++i) { nameType.push_back(tm.nameType(i)); }
    REQUIRE_EQUAL_COLLECTIONS(nameType, nameTypeExpected, "nameType");

    Vec1String       nameSubType;
    const Vec1String nameSubTypeExpected = {"Cl", "O", "S"};
    for (Size i = 0; i < tm.countStructureTypes(); ++i)
    {
        nameSubType.push_back(tm.nameSubType(i));
    }
    REQUIRE_EQUAL_COLLECTIONS(nameSubType, nameSubTypeExpected, "nameSubType");

    Vec1Int toTypeFromString;
    Vec1Int numbers1 = {0, 1, 2, 3, 4};
    for (auto const& s : nameTypeExpected) toTypeFromString.push_back(tm.toType(s));
    REQUIRE_EQUAL_COLLECTIONS(toTypeFromString, numbers1, "toTypeFromString");

    Vec1Int toSubTypeFromString;
    Vec1Int numbers2 = {0, 1, 2};
    for (auto const& s : nameSubTypeExpected) toSubTypeFromString.push_back(tm.toSubType(s));
    REQUIRE_EQUAL_COLLECTIONS(toSubTypeFromString, numbers2, "toSubTypeFromString");
}

BOOST_AUTO_TEST_CASE(IllegalTypeCombination_ThrowError)
{
    BOOST_REQUIRE_THROW(TypeMap({"H", "O", "Na", "Cl", "S"}, {"Ca", "O", "S"}), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
