#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE utils

#include "boost_helpers.hpp"

#include "debug.hpp"
#include "utils.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <stdexcept>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(ScanInString_FindNextNonWhitespace)
{
    String input = "   \t \n   \f\v \r N \n ex \tt \t \n";
    String::const_iterator pos = input.begin();
    String::const_iterator result;

    result = scan(input, pos, flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'N');
    pos = result + 1;
    result = scan(input, pos, flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'e');
    pos = result + 1;
    result = scan(input, pos, flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'x');
    pos = result + 1;
    result = scan(input, pos, flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 't');
    pos = result + 1;
    BOOST_REQUIRE_THROW(scan(input, pos, flf(VASPML_FLF)), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ScanInStringAndExpect_FindExpectedNonWhitespaceOrThrow)
{
    String input = "   \t \n   \f\v \r N \n ex \tt \t \n";
    String::const_iterator pos = input.begin();
    String::const_iterator result;

    BOOST_REQUIRE_THROW(scan(input, pos, flf(VASPML_FLF), "n"), std::runtime_error);
    result = scan(input, pos, flf(VASPML_FLF), "N");
    BOOST_REQUIRE_EQUAL(*result, 'N');
    pos = result + 1;
    result = scan(input, pos, flf(VASPML_FLF), "e");
    BOOST_REQUIRE_EQUAL(*result, 'e');
    pos = result + 1;
    BOOST_REQUIRE_THROW(scan(input, pos, flf(VASPML_FLF), "X"), std::runtime_error);
    result = scan(input, pos, flf(VASPML_FLF), "x");
    BOOST_REQUIRE_EQUAL(*result, 'x');
    pos = result + 1;
    result = scan(input, pos, flf(VASPML_FLF), "t");
    BOOST_REQUIRE_EQUAL(*result, 't');
    pos = result + 1;
    BOOST_REQUIRE_THROW(scan(input, pos, flf(VASPML_FLF)), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ScanInStringDontExpectEndOfString_FindNonWhitespaceOrThrowAtEnd)
{
    String input = "   \t \n   \f\v \r N \n ";
    String::const_iterator pos = input.begin();
    String::const_iterator result;

    result = scan(input, pos, flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'N');
    pos = result + 1;
    String::const_iterator result1 = scan(input, pos, flf(VASPML_FLF), "", false);
    BOOST_REQUIRE(result1 == input.end());
    String::const_iterator result2 = scan(input, pos, flf(VASPML_FLF), "e", false);
    BOOST_REQUIRE(result2 == input.end());
}

BOOST_AUTO_TEST_CASE(FindInString_FindCorrectPosition)
{
    String input = "After the next \nword = there is an equal sign.";
    String::const_iterator pos = input.begin();
    String::const_iterator result;

    result = find(input, pos, "n=", flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'n');
    BOOST_REQUIRE_EQUAL(String(pos, result), "After the ");
    pos = result + 1;
    result = find(input, pos, "n=", flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, '=');
    BOOST_REQUIRE_EQUAL(String(pos, result), "ext \nword ");
    pos = result + 1;
    BOOST_REQUIRE_THROW(find(input, pos, "?", flf(VASPML_FLF)), std::runtime_error);
    result = find(input, pos, "?", flf(VASPML_FLF), false);
    BOOST_REQUIRE(result == input.end());
}

BOOST_AUTO_TEST_CASE(ReverseFindInString_FindCorrectPosition)
{
    String input = "After the next \nword = there is an equal sign.";
    String::const_iterator start = input.begin();
    String::const_iterator pos = input.end();
    String::const_iterator result;

    result = rfind(input, pos, "q=", flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, 'q');
    BOOST_REQUIRE_EQUAL(String(start, result), "After the next \nword = there is an e");
    pos = result - 1;
    result = rfind(input, pos, "q=", flf(VASPML_FLF));
    BOOST_REQUIRE_EQUAL(*result, '=');
    BOOST_REQUIRE_EQUAL(String(start, result), "After the next \nword ");
    pos = result - 1;
    BOOST_REQUIRE_THROW(rfind(input, pos, "?", flf(VASPML_FLF)), std::runtime_error);
    result = rfind(input, pos, "?", flf(VASPML_FLF), false);
    BOOST_REQUIRE(result == input.end());
}

BOOST_AUTO_TEST_CASE(ParseQuotedString_FindClosingQuotes)
{
    String input = "\"";
    String::const_iterator pos = input.begin();
    String result;
    REQUIRE_EXCEPTION_WHAT(string_tools::parseQuotedString(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");

    input = "\"wrong end quotes'";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(string_tools::parseQuotedString(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");

    input = "\"\"empty string, words after end quotes";
    pos = input.begin();
    result = string_tools::parseQuotedString(input, pos);
    BOOST_REQUIRE_EQUAL(result, "");
    BOOST_REQUIRE_EQUAL(*pos, 'e');

    input = "\"\\\"";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(string_tools::parseQuotedString(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");

    input = "\"\\\\\"";
    pos = input.begin();
    result = string_tools::parseQuotedString(input, pos);
    BOOST_REQUIRE_EQUAL(result, "\\\\");
    BOOST_REQUIRE(pos == input.end());

    input = "\"\\\\ word1 \\\\word2\n\"";
    pos = input.begin();
    result = string_tools::parseQuotedString(input, pos);
    BOOST_REQUIRE_EQUAL(result, "\\\\ word1 \\\\word2\n");
    BOOST_REQUIRE(pos == input.end());
}

BOOST_AUTO_TEST_CASE(TrimDifferentStrings_TrimmedStrings)
{
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("  abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("abc def ghi  "), "abc def ghi  ");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("  abc def ghi  "), "abc def ghi  ");

    BOOST_REQUIRE_EQUAL(string_tools::rtrim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("  abc def ghi"), "  abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("abc def ghi  "), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("  abc def ghi  "), "  abc def ghi");

    BOOST_REQUIRE_EQUAL(string_tools::trim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("  abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("abc def ghi  "), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("  abc def ghi  "), "abc def ghi");
}

BOOST_AUTO_TEST_CASE(CreateStringsFromPrintfStyle_CorrectSizesAndContents)
{
    BOOST_REQUIRE_EQUAL(str("%04d", 44), "0044");
    BOOST_REQUIRE_EQUAL(str("%04d", 44).size(), 4);
    BOOST_REQUIRE_EQUAL(str("%4s", "ab"), "  ab");
    BOOST_REQUIRE_EQUAL(str("%4s", "ab").size(), 4);
    BOOST_REQUIRE_EQUAL(str("%-4s", "ab"), "ab  ");
    BOOST_REQUIRE_EQUAL(str(("%04d" + String(STR_MAX - 5, '?')).c_str(), 44),
                        "0044" + String(STR_MAX - 5, '?'));
}

BOOST_AUTO_TEST_CASE(SplitString_CheckCorrectParts)
{
    Vec1String result = string_tools::split("bacadae", "a");
    Vec1String expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "bacadae");

    result = string_tools::split("aabaacaadaaeaa", "aa", false);
    expect = {"", "b", "c", "d", "e", ""};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "aabaacaadaaeaa");

    result = string_tools::split("aabaacaadaaeaa", "aa");
    expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "aabaacaadaaeaa");

    result = string_tools::split("a ; b ; c ; d ; e", ";");
    expect = {"a ", " b ", " c ", " d ", " e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "a ; b ; c ; d ; e");

    result = string_tools::split("abacadaea", "a");
    expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "abacadaea");

    result = string_tools::split("a;b:a,c;:ad,aea:", ";:,", false, true);
    expect = {"a", "b", "a", "c", "", "ad", "aea", ""};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "abacadaea");

    result = string_tools::split("a;b:a,c;:ad,aea:", ";:,", true, true);
    expect = {"a", "b", "a", "c", "ad", "aea"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "abacadaea");
}

BOOST_AUTO_TEST_SUITE_END()
