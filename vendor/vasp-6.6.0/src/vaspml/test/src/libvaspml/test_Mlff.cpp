#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Mlff

#include "test_Mlff.hpp"

#include "boost_helpers.hpp"
#include "FixtureRunInDirectory.hpp"

#include "io.hpp"
#include "rec.hpp"
#include "Record.hpp"
#include "RecordDiff.hpp"
#include "constants.hpp"

#include "boost_unit_test_case.hpp"

#include <fstream>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseMlff> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE_F(FixtureRunInDirectory,
                       MLFFConvertCycleCompareBuffer_SameContent,
                       bdata::make(container.testCases),
                       testCase)
{
    setup(testNameFull() + "_" + testCase.name);

    Record mlffIn ;
    std::fstream strm;

    io::open(strm, testCase.sample->path, "r");
    io::processMlff(mlffIn, strm, false);
    io::close(strm);

    io::open(strm, "ML_FF.out", "w");
    io::processMlff(mlffIn, strm, true);
    io::close(strm);

    Buffer binOrig;
    io::open(strm, testCase.sample->path, "rb");
    io::readBuffer(strm, binOrig);
    io::close(strm);

    Buffer binCopy;
    io::open(strm, "ML_FF.out", "rb");
    io::readBuffer(strm, binCopy);
    io::close(strm);

    REQUIRE_EQUAL_MESSAGE(binOrig.size(), binCopy.size(), "Binary ML_FF comparison: size");

    // Exclude header from byte comparison, there may be small formatting differences (spaces,...).
    Size i = 0;
    for (i = testCase.hasHeader ? constants::mlffHeaderSize : 0; i < binOrig.size(); ++i)
    {
        if (binOrig[i] != binCopy[i]) break;
    }
    REQUIRE_EQUAL_MESSAGE(binOrig.size(), i, "Binary ML_FF comparison: first different byte");
}

BOOST_DATA_TEST_CASE_F(FixtureRunInDirectory,
                       MLFFConvertCycle_EmptyDiff,
                       bdata::make(container.testCases),
                       testCase)
{
    setup(testNameFull() + "_" + testCase.name);

    Record mlffIn ;
    std::fstream strm;

    io::open(strm, testCase.sample->path, "r");
    io::processMlff(mlffIn, strm, false);
    io::close(strm);

    Record mlffOut = mlffIn;
    io::open(strm, "ML_FF.out", "w");
    io::updateMlffHeaderDate(mlffOut);
    io::processMlff(mlffOut, strm, true);
    io::close(strm);

    Record mlffInAfterWrite;
    io::open(strm, "ML_FF.out", "r");
    io::processMlff(mlffInAfterWrite, strm, false);
    io::close(strm);

    Record diff = rec::diff(mlffIn, mlffInAfterWrite);

    BOOST_REQUIRE_MESSAGE(rec::detail::emptyRecordDiffLocal(diff),
                          "No differences in top level of ML_FFs");
    if (testCase.hasHeader)
    {
        const Record& header = *diff.cget<ShRec>("common-sub")->cget<ShRec>("header");
        BOOST_REQUIRE_MESSAGE(header.cget<ShRec>("only-in-1")->empty(),
                              "No items only in first dict");
        BOOST_REQUIRE_MESSAGE(header.cget<ShRec>("only-in-2")->empty(),
                              "No items only in second dict");
        const Record& localDiff = *header.cget<ShRec>("local-diff");
        REQUIRE_EQUAL_MESSAGE(localDiff.cget<Vec1String>("KEYS").size(),
                              1,
                              "Only single entry expected");
        REQUIRE_EQUAL_MESSAGE(localDiff.cget<Vec1String>("KEYS")[0],
                              "date",
                              "Date should be the only different item in header dict");
        BOOST_REQUIRE_MESSAGE(localDiff.contains("date-1"),
                              "Diff should contain date from first dict");
        BOOST_REQUIRE_MESSAGE(localDiff.contains("date-2"),
                              "Diff should contain date from second dict");
    }
}

BOOST_AUTO_TEST_SUITE_END()
