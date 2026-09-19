#include "FixtureRunInDirectory.hpp"

#include "boost_unit_test.hpp"

#include <regex>

using namespace vaspml;
namespace fs = std::filesystem;

const String prefix = "VASPML_RUN_TEST_";

FixtureRunInDirectory::FixtureRunInDirectory() :
    Fixture("Run in " + prefix + "... subdirectory"),
    setupCalled(false),
    teardownCalled(false)
{
    bool inside = inTestDirectory();
    BOOST_REQUIRE_MESSAGE(!inside,
                          "Attempt to construct \"Run in subdirectory\" fixture inside existing "
                          "test run directory \""
                              + String(fs::current_path())
                              + "\"). Maybe previous tests did not reset current path correctly?");
    topLevelDirectory = fs::current_path();
}

FixtureRunInDirectory::~FixtureRunInDirectory()
{
    fs::current_path(topLevelDirectory);
    if (setupCalled)
    {
        BOOST_WARN_MESSAGE(
            teardownCalled,
            "Fixture teardown() was not called, assuming a failed test. Test run directory \""
                + runDirectoryName + "\" will be kept for further inspection.");
    }
    else
    {
        BOOST_WARN_MESSAGE(
            teardownCalled,
            "Destructing \"Run in subdirectory\" fixture without ever calling setup().");
    }
}

void FixtureRunInDirectory::setup(String testName)
{
    setupCalled = true;
    runDirectoryName = prefix + sanitizeTestName(testName);
    BOOST_REQUIRE_MESSAGE(!inTestDirectory(),
                          "Attempt to call \"Run in subdirectory\" fixture setup() inside existing "
                          "test run directory \"" + String(fs::current_path()) + "\").");
    BOOST_WARN_MESSAGE(fs::create_directories(fs::path(runDirectoryName)),
                       "Unable to create test directory \"" + runDirectoryName
                           + "\", directory already exists. Test will be continued but may be "
                             "executed in non-empty directory.");
    fs::current_path(runDirectoryName);

    return;
}

void FixtureRunInDirectory::teardown()
{
    teardownCalled = true;
    BOOST_WARN_MESSAGE(setupCalled,
                       "Attempt to call \"Run in subdirectory\" fixture teardown() without "
                       "calling setup() first.");
    if (!setupCalled) return;

    fs::current_path(topLevelDirectory);
    auto deleted = fs::remove_all(fs::path(runDirectoryName));
    BOOST_WARN_MESSAGE(deleted > 0,
                       "Attempt to delete test directory \"" + runDirectoryName + "\" failed.");

    return;
}

bool FixtureRunInDirectory::inTestDirectory()
{
    fs::path currentPath = fs::current_path();

    if (currentPath.begin() == currentPath.end()) return false;

    String dir = (--currentPath.end())->string();
    if (dir.find(prefix) != dir.npos) return true;

    return false;
}

String FixtureRunInDirectory::sanitizeTestName(String testName)
{
    String result = testName;

    result = std::regex_replace(result, std::regex("/"), ".");
    result = std::regex_replace(result, std::regex("\\s+"), "_");

    return result;
}
