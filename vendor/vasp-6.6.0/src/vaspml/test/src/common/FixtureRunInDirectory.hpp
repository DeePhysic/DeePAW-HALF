#ifndef FIXTURERUNINDIRECTORY_HPP
#define FIXTURERUNINDIRECTORY_HPP

#include "Fixture.hpp"

#include "types.hpp"

#include <filesystem>

namespace vaspml
{

struct FixtureRunInDirectory : public Fixture
{
    FixtureRunInDirectory();
    ~FixtureRunInDirectory();
    void   setup(String testName);
    void   teardown();
    bool   inTestDirectory();
    String sanitizeTestName(String);

    bool                  setupCalled;
    bool                  teardownCalled;
    String                runDirectoryName;
    std::filesystem::path topLevelDirectory;
};

}

#endif
