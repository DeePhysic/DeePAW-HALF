#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE io

#include "boost_helpers.hpp"

#include "Record.hpp"
#include "io.hpp"
#include "types.hpp"

#include "boost_unit_test.hpp"

#include <filesystem>
#include <sstream>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(CreateDirectory_Success)
{
    namespace fs = std::filesystem;
    String dirName = "VASPML_RUN_" + testNameFull();

    bool success = io::createDirectory(dirName);

    BOOST_REQUIRE_MESSAGE(success, "function created new directory");
    BOOST_REQUIRE_MESSAGE(fs::exists(fs::path(dirName)), "directory exists");

    fs::remove_all(dirName);
}

BOOST_AUTO_TEST_CASE(RemoveDirectory_Success)
{
    namespace fs = std::filesystem;
    String dirName = "VASPML_RUN_" + testNameFull();
    fs::create_directories(fs::path(dirName));

    Int deleted = io::removeRecursively(dirName);

    REQUIRE_EQUAL_MESSAGE(deleted, 1, "function returns one item deleted");
    BOOST_REQUIRE_MESSAGE(!fs::exists(fs::path(dirName)), "directory is gone");
}

BOOST_AUTO_TEST_CASE(IncarPassRealExpectsInt_WarningIssued)
{
    String rawIncar =
R"raw(ML_LMAX2 = 2.5
)raw";
    std::istringstream strm(rawIncar);

    Record incar;
    io::readIncar(incar, strm);

    Record setup;
    setup.put<bool>("LMLABEXIST", true);
    setup.put<bool>("LMLFFEXIST", true);
    setup.put<Int>("NTYP", 3);
    setup.put<Int>("NSW", 0);
    setup.put<Int>("IBRION", 0);
    io::setupFromIncar(incar, setup);

    //rec::toJson(incar, std::cout);
    //rec::toJson(setup, std::cout);

    BOOST_REQUIRE_EQUAL(setup.cget<Int>("ML_LMAX2"), 2);
}

BOOST_AUTO_TEST_SUITE_END()
