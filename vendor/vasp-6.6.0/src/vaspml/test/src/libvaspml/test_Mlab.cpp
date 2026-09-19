#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Mlab

#include "test_Mlab.hpp"

#include "boost_helpers.hpp"
#include "FixtureRunInDirectory.hpp"

#include "Record.hpp"
#include "RecordDiff.hpp"
#include "io.hpp"
#include "rec.hpp"

#include "boost_unit_test_case.hpp"

#include <fstream>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseMlab> container;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE_F(FixtureRunInDirectory,
                       MLABConvertCycle_EmptyDiff,
                       bdata::make(container.testCases),
                       testCase)
{
    setup(testNameFull() + "_" + testCase.name);

    Record mlabIn ;
    std::fstream strm;

    io::open(strm, testCase.sample->path, "r");
    io::readMlab(mlabIn, strm);
    io::close(strm);

    Record mlabOut = mlabIn;
    io::open(strm, "ML_AB.out", "w");
    io::writeMlab(mlabOut, strm);
    io::close(strm);

    Record mlabInAfterWrite;
    io::open(strm, "ML_AB.out", "r");
    io::readMlab(mlabInAfterWrite, strm);
    io::close(strm);

    Record diff = rec::diff(mlabIn, mlabInAfterWrite);
    BOOST_REQUIRE_MESSAGE(rec::detail::emptyRecordDiffRecursive(diff), "No differences in ML_ABs");
}

BOOST_AUTO_TEST_SUITE_END()
