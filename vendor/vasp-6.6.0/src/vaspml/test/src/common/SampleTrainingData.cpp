#include "SampleTrainingData.hpp"
#include "file_helpers.hpp"

#include "Record.hpp"
#include "utils.hpp"
#include "io.hpp"

#include <fstream>

using namespace vaspml;

SampleTrainingData::SampleTrainingData(String id) :
    SampleFile("set", id),
    mode(string_tools::split(id, "/")[1])
{}

Record& SampleTrainingData::load()
{
    set = std::make_shared<Record>(); 

    requireFileExists(path);

    std::fstream strm;
    io::open(strm, path, "r");
    io::readMlab(*set, strm);
    io::close(strm);

    return *set;
}
