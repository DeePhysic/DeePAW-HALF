#include "Record.hpp"
#include "dataset.hpp"
#include "io.hpp"

#include <fstream>

using namespace vaspml;

int main(int /*argc*/, char** /*argv*/)
{
    std::ifstream fileIn;
    fileIn.open("ML_AB.MAPbI3", std::ifstream::in);
    Record dataset1;
    io::readMlab(dataset1, fileIn);
    fileIn.close();

    fileIn.open("ML_AB.CsPbBr3", std::ifstream::in);
    Record dataset2;
    io::readMlab(dataset2, fileIn);
    fileIn.close();

    dataset::concat(dataset1, dataset2);

    std::ofstream fileOut;
    fileOut.open("ML_AB.concat");
    io::writeMlab(dataset1, fileOut);
    fileOut.close();

    return 0;
}
