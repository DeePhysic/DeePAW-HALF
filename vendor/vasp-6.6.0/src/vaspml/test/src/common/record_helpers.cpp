#include "record_helpers.hpp"

#include <stdexcept>

using namespace vaspml;

Record vaspml::generateSampleRecord(String sampleName)
{
    Record record;

    const Real       energy = 1.2345;
    const Real       cutoff = 5.0123;
    const Int        numAtoms = 100;
    const Int        numTypes = 4;
    const String     system = "Test system";
    const String     subsystem = "Subsystem";
    const bool       typeSort = false;
    const bool       distSort = true;
    const Vec1Real   dist = {1.0, 2.0, 3.0, 4.0, 5.0};
    const Vec1Int    atoms = {1, 4, 6, 8};
    const Vec1String types = {"Fe", "Si", "H", "B", "Cr"};
    const Vec2Real   desc = {
        {-1.0, -2.0, -3.0, -4.0, -5.0},
        {3.4, 1.2, 5.6},
        {5.5, 4.4, 3.3, 2.2}
    };
    const Vec2Int nlm = {
        {-1, -4, -6, -8},
        {1, 3, 4, 7, 9, 11},
        {0, 0, 3},
        {1, 1, 1, 2}
    };
    const Vec2String info = {
        {"ML_RCUT1", "5.0123", "from_incar"},
        {"close", "bracket ] closes", "brace ] closes"},
        {"ENERGY", "1.2345", "ab_initio", "TOTEN"},
        {"types", "Fe Si H B Cr"},
        {"open", "bracket [ opens", "brace { opens"},
        {"descriptions", "This is a longer description of a test system"}
    };

    /*============================================================================================+
     | Nesting in this example:
     |
     |   record
     |   |
     |   |-> ShRec sub
     |   |   |-> justData
     |   |   |-> Vec1ShRec subsub
     |   |       [
     |   |        -> 0 -> ShRec -> justData
     |   |        -> 1 -> ShRec -> justData
     |   |       ]
     |   |
     |   |-> ShRec sub2
     |   |   |-> justData
     |   |   |-> Vec1ShRec subsub
     |   |       [
     |   |        -> 0 -> ShRec -> justData
     |   |        -> 1 -> ShRec -> justData
     |   |        -> 2 -> ShRec -> justData
     |   |       ]
     |   |
     |   |-> ShRec sub3
     |   |   |-> justData
     |   |
     |   |-> Vec1ShRec vsub
     |   |   [
     |   |    -> 0 -> ShRec -> justData
     |   |    -> 1 -> ShRec -> justData
     |   |   ]
     |   |
     |   |-> Vec1ShRec vsub2
     |       [
     |       ]
     |
     +============================================================================================*/
    if (sampleName == "Nested 1")
    {
        record["energy"] = energy;
        record["cutoff"] = cutoff;
        record["numAtoms"] = numAtoms;
        record["numTypes"] = numTypes;
        record["system"] = system;
        record["subsystem"] = subsystem;
        record["typeSort"] = typeSort;
        record["distSort"] = distSort;
        record["dist"] = dist;
        record["atoms"] = atoms;
        record["types"] = types;
        record["desc"] = desc;
        record["nlm"] = nlm;
        record["info"] = info;
        Record justData = record;

        // Create sub-record.
        record.add("sub", "ShRec");
        Record& sub = record.dget<ShRec>("sub");
        // Copy assignment, should have same data now in sub-record.
        sub = justData;
        sub.add("subsub", "Vec1ShRec");
        // Create vector of sub-sub-record and add same content twice (uses copy constructor).
        sub.get<Vec1ShRec>("subsub").push_back(std::make_shared<Record>(justData));
        sub.get<Vec1ShRec>("subsub").push_back(std::make_shared<Record>(justData));

        // Second sub-record.
        record.add("sub2", "ShRec");
        Record& sub2 = record.dget<ShRec>("sub2");
        sub2 = justData;
        // In second sub-record create also vector of sub-sub-record and add content trice.
        sub2.add("sub2sub", "Vec1ShRec");
        sub2.get<Vec1ShRec>("sub2sub").push_back(std::make_shared<Record>(justData));
        sub2.get<Vec1ShRec>("sub2sub").push_back(std::make_shared<Record>(justData));
        sub2.get<Vec1ShRec>("sub2sub").push_back(std::make_shared<Record>(justData));

        // Third sub-record.
        record.add("sub3", "ShRec");
        Record& sub3 = record.dget<ShRec>("sub3");
        sub3 = justData;

        // Vector of sub-records.
        record.add("vsub", "Vec1ShRec");
        record.get<Vec1ShRec>("vsub").push_back(std::make_shared<Record>(justData));
        record.get<Vec1ShRec>("vsub").push_back(std::make_shared<Record>(justData));

        // Empty vector of sub-records.
        record.add("vsub2", "Vec1ShRec");
    }
    else
    {
        throw std::runtime_error("Unknown Record sample name \"" + sampleName + "\".");
    }

    return record;
}
