#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Tutor

#include "Tutor.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <stdexcept>

using namespace vaspml;

// unit test for error
BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(TestErrorHandling_VisualOutputOnly)
{
    Tutor tutor;
    BOOST_REQUIRE_THROW(tutor.error("Testing error message!"), std::runtime_error);
}

// unit test for warning
BOOST_AUTO_TEST_CASE(TestWarningHandling_VisualOutputOnly)
{
    Tutor tutor;
    tutor.warning("Testing warning message!");
    // Dummy test.
    BOOST_REQUIRE_EQUAL(1, 1);
}

// unit test for bug
BOOST_AUTO_TEST_CASE(TestBugHandling_VisualOutputOnly)
{
    Tutor tutor;
    BOOST_REQUIRE_THROW(tutor.bug("Testing bug message!"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
