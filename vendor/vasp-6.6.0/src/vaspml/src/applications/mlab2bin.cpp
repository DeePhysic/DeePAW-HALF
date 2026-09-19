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

    Buffer buffer;
    rec::toBuffer(dataset, buffer);

    std::cout << rec::printMemoryUsage(rec::computeMemory(dataset));
    std::cout << "Size of binary representation of ML_AB: "
              << buffer.size() * sizeof(Buffer::value_type) / 1000. << " kB\n\n";

    std::ofstream fileOut;
    io::open(fileOut, "ML_AB.bin", "wb");
    io::writeBuffer(buffer, fileOut);
    io::close(fileOut);

    Buffer bufferIn;
    io::open(fileIn, "ML_AB.bin", "rb");
    io::readBuffer(fileIn, bufferIn);
    io::close(fileIn);

    Record recFromBuffer;
    rec::fromBuffer(bufferIn, recFromBuffer);

    io::open(fileOut, "ML_AB.write");
    io::writeMlab(recFromBuffer, fileOut);
    io::close(fileOut);

    return 0;
}
