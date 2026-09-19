#include "Record.hpp"
#include "io.hpp"
#include "rec.hpp"
#include "types.hpp"

#include "MemoryCounter.hpp"

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

    std::cout << rec::detail::printRecordByteSizes();
    std::cout << rec::printMemoryUsage(rec::computeMemory(dataset));

    std::ofstream fileOutTmp;
    io::open(fileOutTmp, "ML_AB.json");
    rec::toJson(dataset, fileOutTmp, 4, false);
    io::close(fileOutTmp);

    std::ofstream fileOutTmp2;
    io::open(fileOutTmp2, "ML_AB.out");
    io::writeMlab(dataset, fileOutTmp2);
    io::close(fileOutTmp2);

    //std::ofstream fileOut;
    //io::open(fileOut, "ML_AB.xyz");
    //io::writeXyz(dataset, fileOut);
    //io::close(fileOut);

    /*============================================================================================+
     | Read back in XYZ file.
     +============================================================================================*/
    io::open(fileIn, "ML_AB.xyz");
    Record datasetXyz;
    io::readXyz(datasetXyz, fileIn);
    io::close(fileIn);

    io::open(fileOutTmp, "ML_AB.xyz.json");
    rec::toJson(datasetXyz, fileOutTmp, 4, false);
    io::close(fileOutTmp);

    //io::open(fileOutTmp2, "ML_AB.xyz.out");
    //io::writeMlab(datasetXyz, fileOutTmp2);
    //io::close(fileOutTmp2);

    return 0;
}
