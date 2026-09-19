#include "Record.hpp"
#include "io.hpp"
#include "rec.hpp"
#include "types.hpp"

#include <fstream>
#include <iostream>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    std::ifstream fileIn;
    io::open(fileIn, "ML_AB");
    Record dataset;
    io::readMlab(dataset, fileIn);
    io::close(fileIn);

    std::cout << rec::printMemoryUsage(rec::computeMemory(dataset));

    std::ofstream fileOut;
    io::open(fileOut, "ML_AB.json");
    rec::toJson(dataset, fileOut);
    io::close(fileOut);

    return 0;
}
