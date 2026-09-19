#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE SemanticVersion

#include "SemanticVersion.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <stdexcept>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(SetupVersions_CheckInternals)
{
    SemanticVersion version1(1, 2, 3);
    BOOST_REQUIRE_EQUAL(version1.major(), 1);
    BOOST_REQUIRE_EQUAL(version1.minor(), 2);
    BOOST_REQUIRE_EQUAL(version1.patch(), 3);
    BOOST_REQUIRE_EQUAL(version1.suffix(), "");
    BOOST_REQUIRE_EQUAL(version1.prerelease(), "");
    BOOST_REQUIRE_EQUAL(version1.build(), "");
    BOOST_REQUIRE_EQUAL(version1.kind(), "release");

    SemanticVersion version2(4, 5, 6, "-alpha");
    BOOST_REQUIRE_EQUAL(version2.major(), 4);
    BOOST_REQUIRE_EQUAL(version2.minor(), 5);
    BOOST_REQUIRE_EQUAL(version2.patch(), 6);
    BOOST_REQUIRE_EQUAL(version2.suffix(), "-alpha");
    BOOST_REQUIRE_EQUAL(version2.prerelease(), "alpha");
    BOOST_REQUIRE_EQUAL(version2.build(), "");
    BOOST_REQUIRE_EQUAL(version2.kind(), "pre-release");

    SemanticVersion version3(7, 8, 9, "+abcdefgh");
    BOOST_REQUIRE_EQUAL(version3.major(), 7);
    BOOST_REQUIRE_EQUAL(version3.minor(), 8);
    BOOST_REQUIRE_EQUAL(version3.patch(), 9);
    BOOST_REQUIRE_EQUAL(version3.suffix(), "+abcdefgh");
    BOOST_REQUIRE_EQUAL(version3.prerelease(), "");
    BOOST_REQUIRE_EQUAL(version3.build(), "abcdefgh");
    BOOST_REQUIRE_EQUAL(version3.kind(), "dev");

    SemanticVersion version4(10, 11, 12, "-alpha+abcdefgh");
    BOOST_REQUIRE_EQUAL(version4.major(), 10);
    BOOST_REQUIRE_EQUAL(version4.minor(), 11);
    BOOST_REQUIRE_EQUAL(version4.patch(), 12);
    BOOST_REQUIRE_EQUAL(version4.suffix(), "-alpha+abcdefgh");
    BOOST_REQUIRE_EQUAL(version4.prerelease(), "alpha");
    BOOST_REQUIRE_EQUAL(version4.build(), "abcdefgh");
    BOOST_REQUIRE_EQUAL(version4.kind(), "pre-release-dev");

    SemanticVersion version5(13, 14, 15, "-1+2");
    BOOST_REQUIRE_EQUAL(version5.major(), 13);
    BOOST_REQUIRE_EQUAL(version5.minor(), 14);
    BOOST_REQUIRE_EQUAL(version5.patch(), 15);
    BOOST_REQUIRE_EQUAL(version5.suffix(), "-1+2");
    BOOST_REQUIRE_EQUAL(version5.prerelease(), "1");
    BOOST_REQUIRE_EQUAL(version5.build(), "2");
    BOOST_REQUIRE_EQUAL(version5.kind(), "pre-release-dev");
}

BOOST_AUTO_TEST_CASE(SetupVersions_RefuseInvalid)
{
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "abc"), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(-1, 2, 3), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, -2, 3), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, -3), std::runtime_error);

    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "-+"), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "-1+"), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "-+1"), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "++"), std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion(1, 2, 3, "+abc+def"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(CheckVersionEquality_CorrectResults)
{
    BOOST_REQUIRE(SemanticVersion(1, 2, 3) == SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") == SemanticVersion(1, 2, 3, "-alpha"));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha+abc") == SemanticVersion(1, 2, 3, "-alpha+abc"));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha+abc") == SemanticVersion(1, 2, 3, "-alpha+def"));

    BOOST_REQUIRE(SemanticVersion(1, 2, 3) != SemanticVersion(1, 2, 0));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3) != SemanticVersion(1, 0, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3) != SemanticVersion(0, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") != SemanticVersion(1, 2, 3, "-beta"));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha+abc") != SemanticVersion(1, 2, 3, "-beta+abc"));
}

BOOST_AUTO_TEST_CASE(CheckVersionSmaller_CorrectResults)
{
    BOOST_REQUIRE(SemanticVersion(1, 2, 0) < SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") < SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") < SemanticVersion(1, 2, 3, "-beta"));
    BOOST_REQUIRE(SemanticVersion(1, 2, 0, "-beta") < SemanticVersion(1, 2, 3, "-alpha"));
}

BOOST_AUTO_TEST_CASE(CheckVersionSmallerEqual_CorrectResults)
{
    BOOST_REQUIRE(SemanticVersion(1, 2, 0) <= SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") <= SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") <= SemanticVersion(1, 2, 3, "-beta"));
    BOOST_REQUIRE(SemanticVersion(1, 2, 0, "-beta") <= SemanticVersion(1, 2, 3, "-alpha"));

    BOOST_REQUIRE(SemanticVersion(1, 2, 3) <= SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion(1, 2, 3, "-alpha") <= SemanticVersion(1, 2, 3, "-alpha"));
}

BOOST_AUTO_TEST_CASE(CheckVersionStringParsing_CorrectResults)
{
    Int    major;
    Int    minor;
    Int    patch;
    String suffix;

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3-alpha+abc more stuff", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "-alpha+abc");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3-alpha+abc", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "-alpha+abc");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3+alpha-abc", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "+alpha-abc");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3+abc", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "+abc");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3-abc", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "-abc");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3+", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "+");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3-", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "-");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.3", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 3);
    BOOST_REQUIRE_EQUAL(suffix, "");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2.", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 0);
    BOOST_REQUIRE_EQUAL(suffix, "");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.2", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 2);
    BOOST_REQUIRE_EQUAL(patch, 0);
    BOOST_REQUIRE_EQUAL(suffix, "");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1.", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 0);
    BOOST_REQUIRE_EQUAL(patch, 0);
    BOOST_REQUIRE_EQUAL(suffix, "");

    major = 999;
    minor = 999;
    patch = 999;
    suffix = "???";
    SemanticVersion::parseVersionString("1", major, minor, patch, suffix);
    BOOST_REQUIRE_EQUAL(major, 1);
    BOOST_REQUIRE_EQUAL(minor, 0);
    BOOST_REQUIRE_EQUAL(patch, 0);
    BOOST_REQUIRE_EQUAL(suffix, "");

    BOOST_REQUIRE_THROW(SemanticVersion::parseVersionString("1..", major, minor, patch, suffix),
                        std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion::parseVersionString(".", major, minor, patch, suffix),
                        std::runtime_error);
    BOOST_REQUIRE_THROW(SemanticVersion::parseVersionString("", major, minor, patch, suffix),
                        std::runtime_error);
}

BOOST_AUTO_TEST_CASE(CheckStringVsIntConstructor_SameResults)
{
    BOOST_REQUIRE(SemanticVersion("1.2.3-alpha more stuff")
                  == SemanticVersion(1, 2, 3, "-alpha+ef342df other stuff"));
    BOOST_REQUIRE(SemanticVersion("1.2.3+alpha-abc") == SemanticVersion(1, 2, 3, "+alpha-abc"));
    BOOST_REQUIRE(SemanticVersion("1.2.3-alpha+abc") == SemanticVersion(1, 2, 3, "-alpha+abc"));
    BOOST_REQUIRE(SemanticVersion("1.2.3-alpha+abc").build()
                  == SemanticVersion(1, 2, 3, "-alpha+abc").build());
    BOOST_REQUIRE(SemanticVersion("1.2.3+abc") == SemanticVersion(1, 2, 3, "+abc"));
    BOOST_REQUIRE(SemanticVersion("1.2.3+abc").build() == SemanticVersion(1, 2, 3, "+abc").build());
    BOOST_REQUIRE(SemanticVersion("1.2.3") == SemanticVersion(1, 2, 3));
    BOOST_REQUIRE(SemanticVersion("1.2") == SemanticVersion(1, 2));
    BOOST_REQUIRE(SemanticVersion("1") == SemanticVersion(1));
}

BOOST_AUTO_TEST_CASE(UpdateInstance_CorrectValues)
{
    SemanticVersion v("9.9.9-alpha+abc");
    BOOST_REQUIRE(v.update("1.2.3") == SemanticVersion("1.2.3"));
    v.update("9.9.9-alpha+abc");
    BOOST_REQUIRE(v.build() == "abc");
    BOOST_REQUIRE(v.update("1.2.3").build() == "");
}

BOOST_AUTO_TEST_CASE(ConvertToString_CorrectString)
{
    BOOST_REQUIRE(SemanticVersion("9").toString() == "9.0.0");
    BOOST_REQUIRE(SemanticVersion("9.10.11-alpha+abc bla bla").toString() == "9.10.11-alpha+abc");
}

BOOST_AUTO_TEST_SUITE_END()
