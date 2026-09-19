#include "Fixture.hpp"

#include "boost_unit_test.hpp"

using namespace vaspml;

Fixture::Fixture(String contentDescription)
{
    BOOST_TEST_MESSAGE("Loading fixture \"" + contentDescription + "\".");
}
