#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE text

#include "boost_helpers.hpp"

#include "Item.hpp"
#include "BoolFlavor.hpp"
#include "text.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace vaspml;
using BFJ = text::BoolFlavorJson;
using BFE = text::BoolFlavorExtXyz;
using BFI = text::BoolFlavorIncar;

const Real tol = 10.0 * std::numeric_limits<Real>::epsilon();

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(DetectCorrectlyWrittenInt_IntRecognized)
{
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("0"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("10"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("010"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("-10"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("-010"));

    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("1234567890"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("+1234567890"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("-1234567890"));

    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("-2147483648"));
    BOOST_REQUIRE_EQUAL(ItemIndex::INT, text::detectType<BFE>("2147483647"));
}

BOOST_AUTO_TEST_CASE(DetectMalformattedInt_FallbackRecognized)
{
    REQUIRE_EQUAL_MESSAGE(ItemIndex::REAL, text::detectType<BFE>("-2147483649"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::REAL, text::detectType<BFE>("2147483648"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::REAL, text::detectType<BFE>("0x12345"), "Hexadecimal base");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::REAL, text::detectType<BFE>("0x123.ae45P-02"), "Hexadecimal base");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("123abc"), "Trailing characters");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("123d567"), "Mixed in letter");
}

BOOST_AUTO_TEST_CASE(DetectCorrectlyWrittenReal_RealRecognized)
{
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 0."));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 0.0"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 0.1"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 2.1"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+0."));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+0.0"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+0.1"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+2.1"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-0."));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-0.0"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-0.1"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-2.1"));

    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" .12345678901234567890"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+.12345678901234567890"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-.12345678901234567890"));

    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 0.12345E-012"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+0.12345E-012"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-0.12345E-012"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 0.12345e+012"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+0.12345e+012"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-0.12345e+012"));

    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 2.226E-0308"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+2.226E-0308"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-2.226E-0308"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 1.796e+0308"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+1.796e+0308"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-1.796e+0308"));

    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" 12345678901234567890.12345678901234567890E123"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+12345678901234567890.12345678901234567890E+123"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-12345678901234567890.12345678901234567890E-123"));

#ifndef __INTEL_LLVM_COMPILER
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" inf"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" INF"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" infinity"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" INFINITY"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" nan"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" NaN"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>(" NAN"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+inf"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+INF"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+infinity"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+INFINITY"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+nan"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+NaN"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("+NAN"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-inf"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-INF"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-infinity"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-INFINITY"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-nan"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-NaN"));
    BOOST_REQUIRE_EQUAL(ItemIndex::REAL, text::detectType<BFE>("-NAN"));
#endif
}

BOOST_AUTO_TEST_CASE(DetectMalformattedReal_FallbackRecognized)
{
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("+."), "No numbers");

    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(" 2.224E-0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("+2.224E-0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("-2.224E-0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(" 1.798e+0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("+1.798e+0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("-1.798e+0308"), "Out-of-bounds");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("-1.798e+0308"), "Out-of-bounds");

    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("infty"), "Incorrect spelling");

    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(" 0.12345D-012"), "d/D exponent");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("+0.12345D-012"), "d/D exponent");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("-0.12345D-012"), "d/D exponent");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(" 0.12345d+012"), "d/D exponent");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("+0.12345d+012"), "d/D exponent");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("-0.12345d+012"), "d/D exponent");
}

BOOST_AUTO_TEST_CASE(DetectCorrectlyWrittenBoolJson_BoolRecognized)
{
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFJ>("true"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFJ>("false"));
}

BOOST_AUTO_TEST_CASE(DetectMalformattedBoolJson_FallbackRecognized)
{
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("T"), "only letters");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("F"), "only letters");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("True"), "Uppercase T");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("False"), "Uppercase F");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("TRUE"), "All uppercase");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("FALSE"), "All uppercase");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".true."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".false."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".TRUE."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".FALSE."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".TrUe."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>(".fALsE."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("t"), "Lowercase t");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("f"), "Lowercase f");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("tRue"), "Spelling mistake");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFJ>("faLse"), "Spelling mistake");
}

BOOST_AUTO_TEST_CASE(DetectCorrectlyWrittenBoolExtXyz_BoolRecognized)
{
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("T"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("F"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("true"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("false"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("True"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("False"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("TRUE"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFE>("FALSE"));
}

BOOST_AUTO_TEST_CASE(DetectMalformattedBoolExtXyz_FallbackRecognized)
{
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("t"), "Lowercase t");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("f"), "Lowercase f");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("tRue"), "Spelling mistake");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>("faLse"), "Spelling mistake");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(".true."), "Fortran style");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFE>(".false."), "Fortran style");
}

BOOST_AUTO_TEST_CASE(DetectCorrectlyWrittenBoolIncar_BoolRecognized)
{
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>("T"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>("F"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>("Totally"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>("Finally"));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".true."));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".false."));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".TRUE."));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".FALSE."));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".TrUe."));
    BOOST_REQUIRE_EQUAL(ItemIndex::BOOL, text::detectType<BFI>(".fALsE."));
}

BOOST_AUTO_TEST_CASE(DetectMalformattedBoolIncar_FallbackRecognized)
{
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFI>("."), "Only a dot");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFI>("drue"), "Spelling mistake");
    REQUIRE_EQUAL_MESSAGE(ItemIndex::STRING, text::detectType<BFI>("valse"), "Spelling mistake");
}

BOOST_AUTO_TEST_CASE(ConvertCorrectlyWrittenIntegers_CorrectValues)
{
    BOOST_REQUIRE_EQUAL(0, text::convert<Int>("0"));
    BOOST_REQUIRE_EQUAL(10, text::convert<Int>("10"));
    BOOST_REQUIRE_EQUAL(10, text::convert<Int>("010"));
    BOOST_REQUIRE_EQUAL(-10, text::convert<Int>("-10"));
    BOOST_REQUIRE_EQUAL(-10, text::convert<Int>("-010"));

    BOOST_REQUIRE_EQUAL(1234567890, text::convert<Int>("1234567890"));
    BOOST_REQUIRE_EQUAL(1234567890, text::convert<Int>("+1234567890"));
    BOOST_REQUIRE_EQUAL(-1234567890, text::convert<Int>("-1234567890"));

    BOOST_REQUIRE_EQUAL(-2147483648, text::convert<Int>("-2147483648"));
    BOOST_REQUIRE_EQUAL(2147483647, text::convert<Int>("2147483647"));
}

BOOST_AUTO_TEST_CASE(ConvertMalformattedIntegers_Throws)
{
    BOOST_REQUIRE_THROW(text::convert<Int>("-2147483649"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Int>("2147483648"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Int>("0x12345"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Int>("0x123.ae45P-02"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Int>("123abc"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Int>("123d567"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConvertCorrectlyWrittenReal_CorrectValues)
{
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>(" 0."), tol);
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>(" 0.0"), tol);
    BOOST_REQUIRE_SMALL( 0.1 - text::convert<Real>(" 0.1"), tol);
    BOOST_REQUIRE_SMALL( 2.1 - text::convert<Real>(" 2.1"), tol);
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>("+0."), tol);
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>("+0.0"), tol);
    BOOST_REQUIRE_SMALL( 0.1 - text::convert<Real>("+0.1"), tol);
    BOOST_REQUIRE_SMALL( 2.1 - text::convert<Real>("+2.1"), tol);
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>("-0."), tol);
    BOOST_REQUIRE_SMALL( 0.0 - text::convert<Real>("-0.0"), tol);
    BOOST_REQUIRE_SMALL(-0.1 - text::convert<Real>("-0.1"), tol);
    BOOST_REQUIRE_SMALL(-2.1 - text::convert<Real>("-2.1"), tol);

    BOOST_REQUIRE_SMALL( 0.12345678901234567890 - text::convert<Real>(" .12345678901234567890"), tol);
    BOOST_REQUIRE_SMALL( 0.12345678901234567890 - text::convert<Real>("+.12345678901234567890"), tol);
    BOOST_REQUIRE_SMALL(-0.12345678901234567890 - text::convert<Real>("-.12345678901234567890"), tol);

    BOOST_REQUIRE_SMALL( 0.12345E-012 - text::convert<Real>(" 0.12345E-012"), tol);
    BOOST_REQUIRE_SMALL(+0.12345E-012 - text::convert<Real>("+0.12345E-012"), tol);
    BOOST_REQUIRE_SMALL(-0.12345E-012 - text::convert<Real>("-0.12345E-012"), tol);
    BOOST_REQUIRE_SMALL( 0.12345e+012 - text::convert<Real>(" 0.12345e+012"), tol);
    BOOST_REQUIRE_SMALL(+0.12345e+012 - text::convert<Real>("+0.12345e+012"), tol);
    BOOST_REQUIRE_SMALL(-0.12345e+012 - text::convert<Real>("-0.12345e+012"), tol);

    BOOST_REQUIRE_SMALL( 2.226E-0308 - text::convert<Real>(" 2.226E-0308"), tol);
    BOOST_REQUIRE_SMALL(+2.226E-0308 - text::convert<Real>("+2.226E-0308"), tol);
    BOOST_REQUIRE_SMALL(-2.226E-0308 - text::convert<Real>("-2.226E-0308"), tol);
    BOOST_REQUIRE_SMALL( 1.796e+0308 - text::convert<Real>(" 1.796e+0308"), tol);
    BOOST_REQUIRE_SMALL(+1.796e+0308 - text::convert<Real>("+1.796e+0308"), tol);
    BOOST_REQUIRE_SMALL(-1.796e+0308 - text::convert<Real>("-1.796e+0308"), tol);

    BOOST_REQUIRE_SMALL( 12345678901234567890.12345678901234567890E123 - text::convert<Real>(" 12345678901234567890.12345678901234567890E123"), tol);
    BOOST_REQUIRE_SMALL(+12345678901234567890.12345678901234567890E+123 - text::convert<Real>("+12345678901234567890.12345678901234567890E+123"), tol);
    BOOST_REQUIRE_SMALL(-12345678901234567890.12345678901234567890E-123 - text::convert<Real>("-12345678901234567890.12345678901234567890E-123"), tol);

#ifndef __INTEL_LLVM_COMPILER
    BOOST_REQUIRE(std::isinf(text::convert<Real>(" inf")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>(" INF")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>(" infinity")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>(" INFINITY")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>(" nan")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>(" NaN")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>(" NAN")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("+inf")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("+INF")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("+infinity")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("+INFINITY")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("+nan")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("+NaN")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("+NAN")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("-inf")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("-INF")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("-infinity")));
    BOOST_REQUIRE(std::isinf(text::convert<Real>("-INFINITY")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("-nan")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("-NaN")));
    BOOST_REQUIRE(std::isnan(text::convert<Real>("-NAN")));
#endif
}

BOOST_AUTO_TEST_CASE(ConvertMalformattedReal_Throws)
{
    BOOST_REQUIRE_THROW(text::convert<Real>(" 2.224E-0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>("+2.224E-0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>("-2.224E-0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>(" 1.798e+0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>("+1.798e+0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>("-1.798e+0308"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<Real>("-1.798e+0308"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConvertCorrectlyWrittenBoolJson_CorrectValues)
{
    BOOST_REQUIRE_EQUAL(true, text::convert<BFJ>("true"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFJ>("false"));
}

BOOST_AUTO_TEST_CASE(ConvertMalformattedBoolJson_Throws)
{
    BOOST_REQUIRE_THROW(text::convert<BFJ>("T"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("F"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("True"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("False"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("TRUE"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("FALSE"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".true."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".false."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".TRUE."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".FALSE."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".TrUe."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>(".fALsE."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("t"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("f"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("tRue"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFJ>("faLse"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConvertCorrectlyWrittenBoolExtXyz_CorrectValues)
{
    BOOST_REQUIRE_EQUAL(true, text::convert<BFE>("T"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFE>("F"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFE>("true"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFE>("false"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFE>("True"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFE>("False"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFE>("TRUE"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFE>("FALSE"));
}

BOOST_AUTO_TEST_CASE(ConvertMalformattedBoolExtXyz_Throws)
{
    BOOST_REQUIRE_THROW(text::convert<BFE>("t"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFE>("f"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFE>("tRue"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFE>("faLse"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFE>(".TrUe"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFE>("fALsE."), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConvertCorrectlyWrittenBoolIncar_CorrectValues)
{
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>("T"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>("F"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>("Totally"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>("Finally"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>("true"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>("false"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>("True"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>("False"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>("TRUE"));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>("FALSE"));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>(".true."));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>(".false."));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>(".TRUE."));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>(".FALSE."));
    BOOST_REQUIRE_EQUAL(true, text::convert<BFI>(".TrUe."));
    BOOST_REQUIRE_EQUAL(false, text::convert<BFI>(".fALsE."));
}

BOOST_AUTO_TEST_CASE(ConvertMalformattedBoolIncar_Throws)
{
    BOOST_REQUIRE_THROW(text::convert<BFI>("."), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFI>("drue"), std::runtime_error);
    BOOST_REQUIRE_THROW(text::convert<BFI>("valse"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
