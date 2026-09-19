#ifndef FILE_HELPERS_HPP
#define FILE_HELPERS_HPP

#include "types.hpp"

#include "boost_unit_test.hpp"

#include <filesystem>

namespace vaspml
{

inline void requireFileExists(String path)
{
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(path),
                          "Cannot find file for this test case (relative path \"" + path
                              + "\") from  current directory \""
                              + String(std::filesystem::current_path()) + "\".");
}

}; // namespace vaspml

#endif

