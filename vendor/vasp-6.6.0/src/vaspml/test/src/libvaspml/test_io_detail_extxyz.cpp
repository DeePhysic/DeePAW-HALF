#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE io_detail_extxyz

#include "Record.hpp"
#include "SmartEnum.hpp"
#include "io_detail_extxyz.hpp"
#include "types.hpp"

#include "boost_helpers.hpp"

#include "boost_unit_test.hpp"

#include <stdexcept>

using namespace vaspml;
using namespace vaspml::io::detail::extxyz;

struct FixturePrimitiveTypes
{
    Vec1String arrayMarkers{"\"\" ", "{} ", "[], "};

    String bs{"bare_string"};

    String qs1{"\"quoted_string\""};
    String p_qs1{"quoted_string"};
    String qs2{"\"qu\\\\oted \\\"=\\\" \tstr\\i\\ng\""};
    String p_qs2{"qu\\oted \"=\" \tstri\ng"};

    String r1{"-0.12345E-03"};
    Real p_r1 = -0.12345E-03;
    String r2{".123456"};
    Real p_r2 = 0.123456;

    String i1{"10"};
    Int p_i1 = 10;
    String i2{"-1234567890"};
    Int p_i2 = -1234567890;

    String b1{"True"};
    bool p_b1 = true;
    String b2{"FALSE"};
    bool p_b2 = false;

    String newKey(Int& i, String extra = "")
    {
        ++i;
        return String{"key-" + std::to_string(i) + (extra.empty() ? "" : "_" + extra)};
    };
};


BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(DetectValueType_CorrectType)
{
    using VT = ValueType;
    BOOST_REQUIRE_EQUAL(toString(VT::Integer), toString(valueType("10")));
    BOOST_REQUIRE_EQUAL(toString(VT::Float), toString(valueType("-0.1")));
    BOOST_REQUIRE_EQUAL(toString(VT::Bool), toString(valueType("TRUE")));
    BOOST_REQUIRE_EQUAL(toString(VT::BareString), toString(valueType("bare")));
    BOOST_REQUIRE_EQUAL(toString(VT::QuotedString), toString(valueType("\"quoted string\"")));
}

BOOST_AUTO_TEST_CASE(ParseBareString_CorrectStringOrException)
{
    String input = "myBareString abc";
    String::const_iterator pos = input.begin();
    String result = parseBareString(input, pos, "=");
    BOOST_REQUIRE_EQUAL(result, "myBareString");
    BOOST_REQUIRE_EQUAL(*pos, ' ');

    input = "myBareString";
    pos = input.begin();
    result = parseBareString(input, pos, "=");
    BOOST_REQUIRE_EQUAL(result, "myBareString");
    BOOST_REQUIRE(pos == input.end());

    input = "";
    pos = input.begin();
    result = parseBareString(input, pos, "=");
    BOOST_REQUIRE_EQUAL(result, "");
    BOOST_REQUIRE(pos == input.begin());
    BOOST_REQUIRE(pos == input.end());

    input = "myBareString";
    pos = input.begin();
    result = parseBareString(input, pos, "r");
    BOOST_REQUIRE_EQUAL(result, "myBa");
    BOOST_REQUIRE_EQUAL(*pos, 'r');

    input = "myBa\"reString = 123";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseBareString(input, pos, "="),
                           std::runtime_error,
                           "1 character(s) which require double quotes");

    input = "my,Ba\"{re}[Str]i=ng 123";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseBareString(input, pos, "="),
                           std::runtime_error,
                           "6 character(s) which require double quotes");
}

BOOST_AUTO_TEST_CASE(RemoveMeaninglessEscapes_CorrectString)
{
    String input = "";
    String result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "");

    input = "\\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "");

    input = " \\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " ");

    input = "\\ ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " ");

    input = "\\\\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\");

    input = " \\\\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\");

    input = "\\\\ ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\ ");

    input = "\\\\\\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\");

    input = " \\\\\\";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\");

    input = "\\\\\\ ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\ ");

    input = " \\\\\\ ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\ ");


    input = "\\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\"");

    input = " \\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\"");

    input = "\\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\" ");

    input = " \\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\" ");

    input = "\\\\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\"");

    input = " \\\\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\"");

    input = "\\\\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\" ");

    input = " \\\\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\" ");

    input = "\\\\\\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\\\"");

    input = " \\\\\\\"";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\\\"");

    input = "\\\\\\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\\\" ");

    input = " \\\\\\\" ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\\\" ");


    input = "\\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\n");

    input = " \\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\n");

    input = "\\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\n ");

    input = " \\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\n ");

    input = "\\\\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\n");

    input = " \\\\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\n");

    input = "\\\\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\n ");

    input = " \\\\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\n ");

    input = "\\\\\\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\\n");

    input = " \\\\\\n";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\\n");

    input = "\\\\\\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\\\n ");

    input = " \\\\\\n ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\\\n ");


    input = "\\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "a");

    input = " \\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " a");

    input = "\\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "a ");

    input = " \\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " a ");

    input = "\\\\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\a");

    input = " \\\\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\a");

    input = "\\\\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\a ");

    input = " \\\\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\a ");

    input = "\\\\\\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\a");

    input = " \\\\\\a";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\a");

    input = "\\\\\\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "\\\\a ");

    input = " \\\\\\a ";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, " \\\\a ");

    input = "abc\\\"def\\ghi\\\\jklm\\nop";
    result = removeMeaninglessEscapes(input);
    BOOST_REQUIRE_EQUAL(result, "abc\\\"defghi\\\\jklm\\nop");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueLegacyQuotesInvalid_CorrectExceptionThrown)
{
    String input;
    String::const_iterator pos;

    input = "\" n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueLegacyQuotes_CorrectIntermediateStructure)
{
    using clm = CommentLineMode;

    String input = "\"\"";
    String::const_iterator pos = input.begin();
    CommentLineValue result = parseCommentLineValue(input, pos);
    Vec2String expect{{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyQuotes));
    BOOST_REQUIRE_EQUAL(result.dim, -1);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy quotes: empty");

    input = "\"   \"";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyQuotes));
    BOOST_REQUIRE_EQUAL(result.dim, -1);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy quotes: empty (spaces)");

    input = "\"1.0  0.0    2.0\"";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"1.0", "0.0", "2.0"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyQuotes));
    BOOST_REQUIRE_EQUAL(result.dim, -1);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy quotes: vector with extra spaces");

    input = "\"  abc \n123 \t b\ra ? !  \"n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"abc", "123", "b", "a", "?", "!"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyQuotes));
    BOOST_REQUIRE_EQUAL(result.dim, -1);
    BOOST_REQUIRE(*pos == 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy quotes: special chars");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueLegacyBracesInvalid_CorrectExceptionThrown)
{
    String input;
    String::const_iterator pos;

    input = "{ n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Unexpected end of content encountered");

    input = "{ \"abc  }n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");

    input = "{ \"abc\"\"def\"  }n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "back-to-back quoted strings without delimiter encountered");

    input = "{ \"123\"a}n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "additional content after closing double quote detected");

    input = "{a\"123\" }n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");

    input = "{\"123\" a\"456\"}n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueLegacyBraces_CorrectIntermediateStructure)
{
    using clm = CommentLineMode;

    String input = "{}";
    String::const_iterator pos = input.begin();
    CommentLineValue result = parseCommentLineValue(input, pos);
    Vec2String expect{{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: empty (pos at end)");

    input = "{}n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: empty");

    input = "{  }n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: empty (spaces)");

    input = "{abc}n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"abc"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: 1 item");

    input = "{ abc    def}n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"abc", "def"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: 2 items (bare)");

    input = "{ \"abc\" \"def\" }n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"\"abc\"", "\"def\""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: 2 items (quoted)");

    input = "{ ?!\\ \"d}ef}\n{}\" }n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"?!\\", "\"d}ef}\n{}\""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: 2 item (special chars)");

    input = "{1.0  0.0  0.0}n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"1.0", "0.0", "0.0"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::LegacyBraces));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Legacy braces: unit vector");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueNewArray1DInvalid_CorrectExceptionThrown)
{
    String input;
    String::const_iterator pos;

    input = "[ n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Unexpected end of content encountered");

    input = "[ \"abc  ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "unable to find closing double quotes");

    input = "[ \"abc\"\"def\"  ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "back-to-back quoted strings without delimiter encountered");

    input = "[ \"abc\" abc ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(
        parseCommentLineValue(input, pos),
        std::runtime_error,
        "additional content after closing double quote detected (before closing enclosure)");

    input = "[a\"123\" ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");

    input = "[\"123\", a\"456\"]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");

    input = "[ \"123\" 123, \"abc\" ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(
        parseCommentLineValue(input, pos),
        std::runtime_error,
        "additional content after closing double quote detected (before delimiter)");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueNewArray1D_CorrectIntermediateStructure)
{
    using clm = CommentLineMode;

    String input = "[]";
    String::const_iterator pos = input.begin();
    CommentLineValue result = parseCommentLineValue(input, pos);
    Vec2String expect{{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty (pos at end)");

    input = "[]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty");

    input = "[  ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty (spaces)");

    input = "[,]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"", ""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: 1 x 2");

    input = "[  ,   ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"", ""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: 1 x 2 (spaces)");

    input = "[a,bc,def]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"a", "bc", "def"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: 1 x 3 (bare strings)");

    input = "[\"][abc?!\",  123  456 ? , \"quoted\" ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"\"][abc?!\"", "123  456 ?", "\"quoted\""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: 1 x 3 (mixed strings)");

    input = "[1.0, 0.0, 0.0]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"1.0", "0.0", "0.0"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 1);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: 1 x 3 (unit vector)");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueNewArray2DInvalid_CorrectExceptionThrown)
{
    String input;
    String::const_iterator pos;

    input = "[ [ ";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Unexpected end of content encountered");

    input = "[ [ ]";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Unexpected end of content encountered");

    input = "[   ,[  ] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "switch from 1D to 2D");

    input = "[1[  ] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before opening enclosure");

    input = "[[  ],1[ ]]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before opening enclosure");

    input = "[[  ] 1[ ]]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "missing comma before opening enclosure");

    input = "[[  ] [ ]]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "missing comma before opening enclosure");


    input = "[[\"abc\"\"def\"]]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "back-to-back quoted strings without delimiter encountered");

    input = "[[  ], [ ], ,]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "expected opening enclosure after delimiter, got another delimiter");

    input = "[ [abc   \"123\"] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");

    input = "[ [a\"123\"] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "Non-whitespace characters found before start of quoted string");

    input = "[ [\"abc\" abc] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "additional content after closing double quote detected");

    input = "[ [\"abc\" abc, def] ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "additional content after closing double quote detected");

    input = "[[  ], [ ], ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "expected opening enclosure after delimiter, got closing enclosure");

    input = "[[  ], [ ], \"abc\" ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "switch from 2D to 1D, got quoted string");

    input = "[[  ], [ ], abc ]n";
    pos = input.begin();
    REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                           std::runtime_error,
                           "expected opening enclosure after delimiter, got closing enclosure");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueNewArray2D_CorrectIntermediateStructure)
{
    using clm = CommentLineMode;

    String input = "[[]]";
    String::const_iterator pos = input.begin();
    CommentLineValue result = parseCommentLineValue(input, pos);
    Vec2String expect{{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty (pos at end)");

    input = "[[]]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty");

    input = "[   [  ] ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: empty (spaces)");

    input = "[   [  ], [,], [,,] ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{""}, {"", ""}, {"", "", ""}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: multiple empty");

    input = "[ [ 123, \"456\" ], [ ? , , \"\\[]]\"], [a,b,c,d,e,f] ]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"123", "\"456\""}, {"?", "", "\"\\[]]\""}, {"a", "b", "c", "d", "e", "f"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: different sizes, special chars");

    input = "[[1.0, 0.0, 0.0], [0.0,1.0,0.0], [0.0, 0.0, 1.0]]n";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"1.0", "0.0", "0.0"}, {"0.0", "1.0", "0.0"}, {"0.0", "0.0", "1.0"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::NewArray));
    BOOST_REQUIRE_EQUAL(result.dim, 2);
    BOOST_REQUIRE_EQUAL(*pos, 'n');
    BOOST_REQUIRE_EQUAL(result.full, input.substr(0, input.size() - 1));
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "New array: unit matrix");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueBareStringInvalid_CorrectExceptionThrown)
{
    String input;
    String::const_iterator pos;

    for (auto const& s : illegalCharsInBareString())
    {
        input = "a" + String{s};
        pos = input.begin();
        BOOST_TEST_INFO("Special char: \"" + String{s} + "\"");
        REQUIRE_EXCEPTION_WHAT(parseCommentLineValue(input, pos),
                               std::runtime_error,
                               "character(s) which require double quotes");
    }
}

BOOST_AUTO_TEST_CASE(ParseCommentLineValueBareString_CorrectIntermediateStructure)
{
    using clm = CommentLineMode;

    String input = "abc";
    String::const_iterator pos = input.begin();
    CommentLineValue result = parseCommentLineValue(input, pos);
    Vec2String expect{{"abc"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::BareString));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE(pos == input.end());
    BOOST_REQUIRE_EQUAL(result.full, input);
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Bare String: simple string");

    input = "abc def";
    pos = input.begin();
    result = parseCommentLineValue(input, pos);
    expect = {{"abc"}};
    BOOST_REQUIRE_EQUAL(toString(result.mode), toString(clm::BareString));
    BOOST_REQUIRE_EQUAL(result.dim, 0);
    BOOST_REQUIRE(*pos == ' ');
    BOOST_REQUIRE_EQUAL(result.full, "abc");
    REQUIRE_EQUAL_COLLECTIONS_2D(result.parts, expect, "Bare String: simple string ends at space");
}

BOOST_AUTO_TEST_CASE(ParseCommentLineBooleansWithoutValue_CorrectItemRetrieved)
{
    Record record;
    String line = "hasEnergy \"no\\\"Stress\\\"\" 123 \"a \\\\= b\" \"flag_with_\\escape\"";

    parseCommentLine(record, line);

    REQUIRE_RECORD_ENTRY(record, "hasEnergy", true);
    REQUIRE_RECORD_ENTRY(record, "no\"Stress\"", true);
    REQUIRE_RECORD_ENTRY(record, "123", true);
    REQUIRE_RECORD_ENTRY(record, "a \\= b", true);
    REQUIRE_RECORD_ENTRY(record, "flag_with_escape", true);
}

BOOST_AUTO_TEST_CASE(ParseCommentLineBooleansWithValue_CorrectItemRetrieved)
{
    Record record;
    String line = "hasEnergy=TRUE \"no\\\"Stress\\\"\" =  false 123  = T \"a \\\\= b\"=False";

    parseCommentLine(record, line);

    REQUIRE_RECORD_ENTRY(record, "hasEnergy", true);
    REQUIRE_RECORD_ENTRY(record, "no\"Stress\"", false);
    REQUIRE_RECORD_ENTRY(record, "123", true);
    REQUIRE_RECORD_ENTRY(record, "a \\= b", false);
}

BOOST_FIXTURE_TEST_CASE(ParseCommentLineVarious_CorrectRecord, FixturePrimitiveTypes)
{
    Record record;
    Int i = 0;
    String line;

    line += " " + newKey(i);
    line += " " + newKey(i) + " = " + bs;
    line += " " + newKey(i) + " = " + qs1;
    line += " " + newKey(i) + " = " + qs2;
    line += " " + newKey(i) + " = " + r1;
    line += " " + newKey(i) + " = " + r2;
    line += " " + newKey(i) + " = " + i1;
    line += " " + newKey(i) + " = " + i2;
    line += " " + newKey(i) + " = " + b1;
    line += " " + newKey(i) + " = " + b2;
    for (const String& a : arrayMarkers)
    {
        String s = a.substr(2);
        line += " " + newKey(i) + " = " + a[0] + r1 + s + r2 + s + r1 + s + r2 + a[1];
        line += " " + newKey(i) + " = " + a[0] + i1 + s + i2 + s + i1 + s + i2 + a[1];
        line += " " + newKey(i) + " = " + a[0] + bs + s + r2 + s + i1 + a[1];
        line += " " + newKey(i) + " = " + a[0] + r1 + s + r2 + s + i1 + a[1];
        if (a[0] == '[')
        {
            // Vec2Real
            line += " " + newKey(i) + " = " + a[0];
            line += a[0] + r1 + s + r2 + s + r1 + s + r2 + a[1] + s;
            line += a[0] + r2 + s + r1 + s + r2 + a[1] + s;
            line += a[0] + r1 + s + r2 + s + r2 + s + r1 + s + r1 + a[1] + a[1];
            // Vec2Int
            line += " " + newKey(i) + " = " + a[0];
            line += a[0] + i1 + s + i2 + s + i1 + s + i2 + a[1] + s;
            line += a[0] + i2 + s + i1 + s + i2 + a[1] + s;
            line += a[0] + i1 + s + i2 + s + i2 + s + i1 + s + i1 + a[1] + a[1];
            // mixed, up to BareString => Vec2String
            line += " " + newKey(i) + " = " + a[0];
            line += a[0] + bs + s + r2 + s + i1 + s + bs + a[1] + s;
            line += a[0] + i2 + s + bs + s + i2 + a[1] + s;
            line += a[0] + i1 + s + r2 + s + r2 + s + i1 + s + bs + a[1] + a[1];
            // mixed, including QuotedString => Vec2String
            line += " " + newKey(i) + " = " + a[0];
            line += a[0] + bs + s + r2 + s + i1 + s + bs + a[1] + s;
            line += a[0] + i2 + s + qs1 + s + i2 + a[1] + s;
            line += a[0] + i1 + s + r2 + s + qs2 + s + i1 + s + bs + a[1] + a[1];
        }
    }

    parseCommentLine(record, line);

    i = 0;
    REQUIRE_RECORD_ENTRY(record, newKey(i), true);
    REQUIRE_RECORD_ENTRY(record, newKey(i), bs);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_qs1);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_qs2);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_r1);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_r2);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_i1);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_i2);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_b1);
    REQUIRE_RECORD_ENTRY(record, newKey(i), p_b2);
    for (const String& a : arrayMarkers)
    {
        String s = a.substr(2);
        REQUIRE_RECORD_ENTRY(record, newKey(i), Vec1Real{p_r1, p_r2, p_r1, p_r2});
        REQUIRE_RECORD_ENTRY(record, newKey(i), Vec1Int{p_i1, p_i2, p_i1, p_i2});
        if (a[0] == '"') REQUIRE_RECORD_ENTRY(record, newKey(i), bs + s + r2 + s + i1);
        else REQUIRE_RECORD_ENTRY(record, newKey(i), Vec1String{bs, r2, i1});
        REQUIRE_RECORD_ENTRY(record, newKey(i), Vec1Real{p_r1, p_r2, Real(p_i1)});
        if (a[0] == '[')
        {
            REQUIRE_RECORD_ENTRY(record,
                                 newKey(i),
                                 Vec2Real{
                                     {p_r1, p_r2, p_r1, p_r2},
                                     {p_r2, p_r1, p_r2},
                                     {p_r1, p_r2, p_r2, p_r1, p_r1}
            });
            REQUIRE_RECORD_ENTRY(record,
                                 newKey(i),
                                 Vec2Int{
                                     {p_i1, p_i2, p_i1, p_i2},
                                     {p_i2, p_i1, p_i2},
                                     {p_i1, p_i2, p_i2, p_i1, p_i1}
            });
            REQUIRE_RECORD_ENTRY(record,
                                 newKey(i),
                                 Vec2String{
                                     {bs, r2, i1, bs},
                                     {i2, bs, i2},
                                     {i1, r2, r2, i1, bs}
            });
            REQUIRE_RECORD_ENTRY(record,
                                 newKey(i),
                                 Vec2String{
                                     {bs, r2, i1, bs},
                                     {i2, p_qs1, i2},
                                     {i1, r2, p_qs2, i1, bs}
            });
        }
    }
}

BOOST_FIXTURE_TEST_CASE(ParseCommentLineVariousInvalid_CorrectExceptionThrown,
                        FixturePrimitiveTypes)
{
    Record record;
    Int i = 0;
    String line;

    line = " " + newKey(i) + " = ";

    for (const String& a : arrayMarkers)
    {
        String s = a.substr(2);
        REQUIRE_EXCEPTION_WHAT(
            parseCommentLine(record, line + a[0] + r1 + s + b1 + s + i1 + a[1]),
            std::runtime_error,
            "unsupported array containing a mixture of booleans and integers and/or "
            "floating-point numbers");
    }
}

BOOST_AUTO_TEST_CASE(ParseCommentLineStrings_CorrectItemRetrieved)
{
    Record record;
    String line = "string1 = \" ST\\\"RI!\\\\N\\nG1\" string2 = STRING2 string3 = \" valid INVA=LID \"";

    parseCommentLine(record, line);

    REQUIRE_RECORD_ENTRY(record, "string1", String(" ST\"RI!\\N\nG1")); 
    REQUIRE_RECORD_ENTRY(record, "string2", String("STRING2")); 
    REQUIRE_RECORD_ENTRY(record, "string3", String(" valid INVA=LID ")); 
}

BOOST_AUTO_TEST_CASE(ParseCommentLineRealWorldExample_CorrectItemRetrieved)
{
    Record record;
    String line = "config_type=isolated_atom gap_energy=-157.72725320 dft_virial=\"0.00000000      "
                  " 0.00000000       0.00000000       0.00000000       0.00000000       0.00000000 "
                  "      0.00000000       0.00000000       0.00000000\" dft_energy=-158.54496821 "
                  "nneightol=1.20000000 pbc=\"T T T\" Lattice=\"20.00000000       0.00000000       "
                  "0.00000000       0.00000000      20.00000000       0.00000000       0.00000000  "
                  "     0.00000000      20.00000000\" "
                  "Properties=species:S:1:pos:R:3:Z:I:1:map_shift:I:3:n_neighb:I:1:gap_force:R:3:"
                  "dft_force:R:3";

    parseCommentLine(record, line);

    REQUIRE_RECORD_ENTRY(record, "config_type", String("isolated_atom")); 
    REQUIRE_RECORD_ENTRY(record, "gap_energy", -157.72725320); 
    REQUIRE_RECORD_ENTRY(record, "dft_virial", Vec1Real{0, 0, 0, 0, 0, 0, 0, 0, 0}); 
    REQUIRE_RECORD_ENTRY(record, "dft_energy", -158.54496821); 
    REQUIRE_RECORD_ENTRY(record, "nneightol", 1.20000000); 
    REQUIRE_RECORD_ENTRY(record, "pbc", Vec1String{"T", "T", "T"});
    REQUIRE_RECORD_ENTRY(record, "Lattice", Vec1Real{20, 0, 0, 0, 20, 0, 0, 0, 20});
    REQUIRE_RECORD_ENTRY(
        record,
        "Properties",
        String("species:S:1:pos:R:3:Z:I:1:map_shift:I:3:n_neighb:I:1:gap_force:R:3:dft_force:R:3"));
}

BOOST_AUTO_TEST_SUITE_END()
