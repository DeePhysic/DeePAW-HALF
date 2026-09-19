#include "Record.hpp"
#include "io.hpp"
#include "rec.hpp"
#include "types.hpp"

#include <fstream>
#include <iostream>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    Record       record;
    std::fstream f;

    io::open(f, "ML_FF", "rb");
    io::processMlff(record, f, false);
    io::close(f);

    record.put<Vec1String>("_sort", {"versionMajor", "versionMinor", "versionPatch", "header"});
    std::cout << rec::printMemoryUsage(rec::computeMemory(record));

    io::open(f, "ML_FF.json", "w");
    rec::toJson(record, f);
    io::close(f);

    io::open(f, "ML_FF.out", "wb");
    io::processMlff(record, f, true);
    io::close(f);

    return 0;
}
