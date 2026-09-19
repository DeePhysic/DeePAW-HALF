#include "SampleFile.hpp"
#include "utils.hpp"

using namespace vaspml;

SampleFile::SampleFile(String group, String id) :
    pathPrefix("../../../data/"),
    group(group),
    id(id),
    path(pathPrefix + group + (group.back() == '/' ? "" : "/") + id),
    version(SemanticVersion(string_tools::split(id, "/").front()))
{}
