#include "BoolFlavor.hpp"
#include "Record.hpp"
#include "io.hpp"
#include "rec.hpp"
#include "text.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace vaspml;

int main(int argc, char** argv)
{
    if (argc < 5)
    {
        std::cout << str("USAGE: %s <ML_FF1> <ML_FF2> <log_vector_diffs> <float_tolerance>\n",
                         argv[0]);
        throw std::runtime_error("ERROR: Wrong number of arguments.");
    }

    String file1 = argv[1];
    String file2 = argv[2];
    bool   logVectorDiffs = text::convert<text::BoolFlavorIncar>(String(argv[3]));
    Real   tolerance = text::convert<Real>(String(argv[4]));

    Record mlff1;
    Record mlff2;

    std::fstream f;

    io::open(f, file1, "rb");
    io::processMlff(mlff1, f, false);
    io::close(f);

    io::open(f, file2, "rb");
    io::processMlff(mlff2, f, false);
    io::close(f);

    Record diff = rec::diff(mlff1, mlff2, logVectorDiffs, tolerance);

    io::open(f, "ML_FF-diff.json", "w");
    rec::toJson(diff, f);
    io::close(f);

    return 0;
}
