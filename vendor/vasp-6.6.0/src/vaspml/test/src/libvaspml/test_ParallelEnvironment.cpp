#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE ParallelEnvironment

#include "ParallelEnvironment.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(SetUp_CorrectAlgos)
{
    ParallelEnvironment parallel;

    BOOST_REQUIRE_MESSAGE(parallel.off(), "Default algo: off");
    BOOST_REQUIRE_EQUAL(parallel.selected(), "off");

    parallel.init("serial");
    BOOST_REQUIRE_MESSAGE(!parallel.off(), "Algo serial != off should always be selectable");
    BOOST_REQUIRE_EQUAL(parallel.selected(), "serial");

#ifdef VASPML_PALGO_THREADS
    parallel.init("threads");
    BOOST_REQUIRE(!parallel.off());
    BOOST_REQUIRE(!parallel.gpu());
    BOOST_REQUIRE_EQUAL(parallel.selected(), "threads");
#endif

#ifdef VASPML_PALGO_GPU
    parallel.init("gpu");
    BOOST_REQUIRE(!parallel.off());
    BOOST_REQUIRE(parallel.gpu());
    BOOST_REQUIRE_EQUAL(parallel.selected(), "gpu");
#endif
}

BOOST_AUTO_TEST_CASE(ListSupported_CorrectString)
{
    ParallelEnvironment parallel;

    String expected{"off,serial"};
#ifdef VASPML_PALGO_THREADS
    expected += ",threads";
#endif
#ifdef VASPML_PALGO_GPU
    expected += ",gpu";
#endif

    BOOST_REQUIRE_EQUAL(expected, parallel.supported());
}

BOOST_AUTO_TEST_SUITE_END()

