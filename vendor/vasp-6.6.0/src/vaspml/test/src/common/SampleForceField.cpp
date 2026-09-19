#include "file_helpers.hpp"
#include "SampleForceField.hpp"

#include "Record.hpp"
#include "utils.hpp"
#include "io.hpp"

#include <fstream>

using namespace vaspml;

SampleForceField::SampleForceField(String id) :
    SampleFile("ff", id),
    mode(string_tools::split(id, "/")[1]),
    fast(mode == "refit" ? true : false)
{}

Record& SampleForceField::load()
{
    ff = std::make_shared<Record>(); 

    requireFileExists(path);

    std::fstream strm;
    io::open(strm, path, "r");
    io::processMlff(*ff, strm, false);
    io::close(strm);

    return *ff;
}
