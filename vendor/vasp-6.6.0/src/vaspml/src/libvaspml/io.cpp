#include "io.hpp"

#include "Item.hpp"
#include "Record.hpp"
#include "SemanticVersion.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "buffer.hpp"
#include "constants.hpp"
#include "dataset.hpp"
#include "debug.hpp"
#include "io_detail.hpp"
#include "io_detail_extxyz.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm>  // for max, copy
#include <exception>  // for exception
#include <filesystem> // for path, create_directories, remove_all
#include <iostream>   // for cout
#include <iterator>   // for istream_iterator, operator!=
#include <stdexcept>  // for logic_error, runtime_error

using namespace vaspml;
using namespace vaspml::io::detail;

template<>
SmartEnum<io::TagState>::EnumMap SmartEnum<io::TagState>::mapEnumsNames()
{
    return SmartEnum<io::TagState>::EnumMap{
        {io::TagState::UNSET,     "(?)"},
        {io::TagState::DEFAULT,   "   "},
        {io::TagState::INCAR,     "(I)"},
        {io::TagState::INCAR_ALT, "(i)"},
        {io::TagState::OVERRIDE,  "(!)"}
    };
}

void io::open(std::fstream& strm, const String& fileName, const String& mode)
{
    std::ios_base::openmode imode = std::ios_base::in;

    if (mode == "r") imode = std::ios_base::in;
    else if (mode == "rb") imode = std::ios_base::in | std::ios_base::binary;
    else if (mode == "w") imode = std::ios_base::out;
    else if (mode == "wb") imode = std::ios_base::out | std::ios_base::binary;
    else if (mode == "a") imode = std::ios_base::app;
    else if (mode == "ab") imode = std::ios_base::app | std::ios_base::binary;
    else global::tutor.bug("Unknown file open mode \"" + mode + "\".");

    strm.open(fileName, imode);

    if (!strm.is_open()) { global::tutor.error("Could not open file \"" + fileName + "\"."); }

    return;
}

void io::open(std::ifstream& strm, const String& fileName, const String& mode)
{
    std::ios_base::openmode imode = std::ios_base::in;

    if (mode == "r") imode = std::ios_base::in;
    else if (mode == "rb") imode = std::ios_base::in | std::ios_base::binary;
    else global::tutor.bug("Unknown file open mode \"" + mode + "\".");

    strm.open(fileName, imode);

    if (!strm.is_open()) { global::tutor.error("Could not open file \"" + fileName + "\"."); }

    return;
}

void io::open(std::ofstream& strm, const String& fileName, const String& mode)
{
    std::ios_base::openmode imode = std::ios_base::in;

    if (mode == "w") imode = std::ios_base::out;
    else if (mode == "wb") imode = std::ios_base::out | std::ios_base::binary;
    else if (mode == "a") imode = std::ios_base::app;
    else if (mode == "ab") imode = std::ios_base::app | std::ios_base::binary;
    else global::tutor.bug("Unknown file open mode \"" + mode + "\".");

    strm.open(fileName, imode);

    if (!strm.is_open()) { global::tutor.error("Could not open file \"" + fileName + "\"."); }

    return;
}

void io::close(std::fstream& strm)
{
    if (!strm.is_open()) { global::tutor.bug("Invalid attempt to close already closed fstream"); }
    strm.close();

    return;
}

void io::close(std::ifstream& strm)
{
    if (!strm.is_open()) { global::tutor.bug("Invalid attempt to close already closed ifstream"); }
    strm.close();

    return;
}

void io::close(std::ofstream& strm)
{
    if (!strm.is_open()) { global::tutor.bug("Invalid attempt to close already closed ofstream"); }
    strm.close();

    return;
}

bool io::createDirectory(String directoryPath)
{
    namespace fs = std::filesystem;
    bool success = false;

    try
    {
        success = fs::create_directories(fs::path(directoryPath));
    }
    catch (std::exception& e)
    {
        global::tutor.error("Cannot create directory \"" + directoryPath + "\" (" + e.what()
                            + ").");
    }

    return success;
}

Int io::removeRecursively(String fileOrDirectory)
{
    namespace fs = std::filesystem;
    Int deleted = 0;

    try
    {
        deleted = fs::remove_all(fs::path(fileOrDirectory));
    }
    catch (std::exception& e)
    {
        global::tutor.error("Cannot remove \"" + fileOrDirectory + "\" (" + e.what() + ").");
    }

    return deleted;
}

void io::writeBuffer(const Buffer& buffer, std::ostream& strm)
{
    strm.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());

    return;
}

void io::readBuffer(std::istream& strm, Buffer& buffer)
{
    // Read file into buffer.
    buffer.resize(bytesLeftInStream(strm));
    if (buffer.size() == 0) { global::tutor.error("Cannot process stream because it is empty."); }
    // Avoid skipping whitespaces!
    std::noskipws(strm);
    std::copy(std::istream_iterator<char>(strm),
              std::istream_iterator<char>(),
              reinterpret_cast<char*>(buffer.data()));

    return;
}

String io::toString(std::istream& strm)
{
    constexpr Size bufferSize = 4096;
    String         result;

    String buffer = std::string(bufferSize, '\0');
    while (strm.read(&buffer[0], bufferSize)) { result.append(buffer, 0, strm.gcount()); }
    result.append(buffer, 0, strm.gcount());

    return result;
}

void io::readMlab(Record& record, std::istream& input)
{
    if (!record.empty()) global::tutor.error("Cannot read ML_AB file, record is not empty.");

    using namespace string_tools;
    String line;
    Int    n = 0;

    // Read first line and check for version string.
    n += readLine(input, line);
    if (!checkForSubstring(line, "Version"))
    {
        std::runtime_error("ERROR: Invalid header found in input file in ML_AB format.");
    }
    // Store version in record.
    Vec1String      strSplit = split(trim(line), " ");
    SemanticVersion version(strSplit.at(0));
    record["version"] = version.toString();

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The number of configurations", n););
    n += readLine(input, line, 1);
    record["numStructures"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atom type", n););
    n += readLine(input, line, 1);
    record["maxTypes"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The atom types in the data file", n););
    n += readLine(input, line, 0);
    record["typeNames"] = readMlabVector<Vec1String>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atoms per system", n););
    n += readLine(input, line, 1);
    record["maxAtoms"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "The maximum number of atoms per atom type", n););
    n += readLine(input, line, 1);
    record["maxAtomsPerType"] = readMlabScalar<Int>(line);

    n += readLine(input, line, 1);
    VASPML_DEBUG_L2(checkLine(line, "Reference atomic energy (eV)", n););
    n += readLine(input, line, 0);
    record["atomicRefEnergy"] = readMlabVector<Vec1Real>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "Atomic mass", n););
    n += readLine(input, line, 0);
    record["atomicMass"] = readMlabVector<Vec1Real>(input, n);

    n += readLine(input, line, 0);
    VASPML_DEBUG_L2(checkLine(line, "The numbers of basis sets per atom type", n););
    n += readLine(input, line, 0);
    record["numLrc"] = readMlabVector<Vec1Int>(input, n);
    const Int maxTypes = record.cget<Int>("maxTypes");
    record["lrcStructure"] = Vec2Int(maxTypes);
    record["lrcAtom"] = Vec2Int(maxTypes);

    for (Int i = 0; i < maxTypes; ++i)
    {
        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Basis set for", n););
        n += readLine(input, line, 0);
        const Int numLrc = record.cget<Vec1Int>("numLrc").at(i);
        Vec1Int   lrcStructure;
        Vec1Int   lrcAtom;
        for (Int j = 0; j < numLrc; ++j)
        {
            n += readLine(input, line, 0);
            Vec1String values = split(line, " ");
            lrcStructure.push_back(readMlabScalar<Int>(values.at(0)) - 1);
            lrcAtom.push_back(readMlabScalar<Int>(values.at(1)) - 1);
        }
        record.get<Vec2Int>("lrcStructure").at(i) = lrcStructure;
        record.get<Vec2Int>("lrcAtom").at(i) = lrcAtom;
        n += readLine(input, line, 0);
    }

    record["structures"] = Vec1ShRec();
    Vec1ShRec& structures = record.get<Vec1ShRec>("structures");
    const Int  numStructures = record.cget<Int>("numStructures");
    bool       hasCharges = false;
    for (Int i = 0; i < numStructures; ++i)
    {
        structures.push_back(std::make_shared<Record>());
        Record& structure = *structures.back();
        // If charges are not present at the end the line has already been read.
        if (i == 0 || hasCharges) n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Configuration num.", n););
        structure["id"] = (Int)(i + 1);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "System name", n););
        n += readLine(input, line, 1);
        structure["system"] = readMlabScalar<String>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "The number of atom types", n););
        n += readLine(input, line, 1);
        structure["numTypes"] = readMlabScalar<Int>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "The number of atoms", n););
        n += readLine(input, line, 1);
        structure["numAtoms"] = readMlabScalar<Int>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "Atom types and atom numbers", n););
        n += readLine(input, line, 0);
        const Int  numTypes = structure.cget<Int>("numTypes");
        Vec1String types;
        Vec1Int    numAtomsPerType;
        Vec1Int    atomTypes;
        for (Int j = 0; j < numTypes; ++j)
        {
            n += readLine(input, line, 0);
            Vec1String values = split(line, " ");
            types.push_back(readMlabScalar<String>(values.at(0)));
            Int m = readMlabScalar<Int>(values.at(1));
            numAtomsPerType.push_back(m);
            atomTypes.resize(atomTypes.size() + m, j);
        }
        // Contains type names, e.g. [Cs Pb Br].
        structure["typeNames"] = types;
        // Contains number of atoms per type, e.g. [8 24 8].
        structure["numAtomsPerType"] = numAtomsPerType;
        // Contains per-structure type index for each atom, e.g. [0 0 0... 1 1 1 1... 2 2 2..].
        structure["atomTypes"] = atomTypes;

        n += readLine(input, line, 1);
        String tmp = string_tools::trim(readMlabScalar<String>(line));
        if (tmp == "CTIFOR")
        {
            n += readLine(input, line, 1);
            structure["CTIFOR"] = readMlabScalar<Real>(line);
            n += readLine(input, line, 1);
        }
        VASPML_DEBUG_L2(checkLine(line, "Primitive lattice vectors (ang.)", n););
        n += readLine(input, line, 0);
        structure["lattice"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Atomic positions (ang.)", n););
        n += readLine(input, line, 0);
        structure["positions"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Total energy (eV)", n););
        n += readLine(input, line, 1);
        structure["energy"] = readMlabScalar<Real>(line);

        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "Forces (eV ang.^-1)", n););
        n += readLine(input, line, 0);
        structure["forces"] = readMlabVector<Vec1Real>(input, n);

        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "Stress (kbar)", n););
        n += readLine(input, line, 1);
        VASPML_DEBUG_L2(checkLine(line, "XX YY ZZ", n););
        n += readLine(input, line, 0);
        Vec1Real stress = readMlabVector<Vec1Real>(input, n);
        n += readLine(input, line, 0);
        VASPML_DEBUG_L2(checkLine(line, "XY YZ ZX", n););
        n += readLine(input, line, 0);
        Vec1Real stress2 = readMlabVector<Vec1Real>(input, n, (i == numStructures - 1));
        stress.insert(stress.end(), stress2.begin(), stress2.end());
        structure["stress"] = stress;

        Int nextra = readLine(input, line, 0, (i == numStructures - 1));
        if (nextra > 0)
        {
            n += nextra;
            tmp = string_tools::trim(readMlabScalar<String>(line));
            if (tmp == "Charges (e)")
            {
                hasCharges = true;
                n += readLine(input, line, 0);
                structure["charges"] = readMlabVector<Vec1Real>(input, n, (i == numStructures - 1));
            }
            else hasCharges = false;
        }
    }

    dataset::setPerStructureLrc(record.cget<Vec2Int>("lrcStructure"),
                                record.cget<Vec2Int>("lrcAtom"),
                                record);

    return;
}

void io::readMlab(Record& record, String fileName)
{
    std::ifstream file;

    io::open(file, fileName, "r");
    io::readMlab(record, file);
    io::close(file);

    return;
}

void io::writeMlab(const Record& record, std::ostream& output)
{
    const String separatorStar = String(50, '*') + '\n';
    const String separatorDash = String(50, '-') + '\n';
    const String separatorEqual = String(50, '=') + '\n';

    const SemanticVersion version(record.cget<String>("version"));
    output << version.major() << "." << version.minor() << " Version\n";

    output << mlabItemString("The number of configurations",
                             record.cget<Int>("numStructures"),
                             '*',
                             '-');

    output << mlabItemString("The maximum number of atom type",
                             record.cget<Int>("maxTypes"),
                             '*',
                             '-');

    output << mlabItemString("The atom types in the data file",
                             record.cget<Vec1String>("typeNames"),
                             '*',
                             '-');

    output << mlabItemString("The maximum number of atoms per system",
                             record.cget<Int>("maxAtoms"),
                             '*',
                             '-');

    output << mlabItemString("The maximum number of atoms per atom type",
                             record.cget<Int>("maxAtomsPerType"),
                             '*',
                             '-');

    output << mlabItemString("Reference atomic energy (eV)",
                             record.cget<Vec1Real>("atomicRefEnergy"),
                             '*',
                             '-');

    output << mlabItemString("Atomic mass", record.cget<Vec1Real>("atomicMass"), '*', '-');

    output << mlabItemString("The numbers of basis sets per atom type",
                             record.cget<Vec1Int>("numLrc"),
                             '*',
                             '-');

    for (Size i = 0; i < record.cget<Vec1String>("typeNames").size(); ++i)
    {
        output << separatorStar;
        output << "Basis set for " + record.cget<Vec1String>("typeNames")[i] + "\n";
        output << separatorDash;
        const Vec1Int& iStructure = record.cget<Vec2Int>("lrcStructure")[i];
        const Vec1Int& iAtom = record.cget<Vec2Int>("lrcAtom")[i];
        for (Int j = 0; j < record.cget<Vec1Int>("numLrc")[i]; ++j)
        {
            output << str("%d %d\n", iStructure[j] + 1, iAtom[j] + 1);
        }
    }

    for (Int i = 0; i < record.cget<Int>("numStructures"); ++i)
    {
        const Record& s = *record.vcget<Vec1ShRec>("structures")[i];

        output << separatorStar;
        output << str("Configuration num. %d\n", i + 1);

        output << mlabItemString("System name", s.cget<String>("system"), '=', '-');

        output << mlabItemString("The number of atom types", s.cget<Int>("numTypes"), '=', '-');

        output << mlabItemString("The number of atoms", s.cget<Int>("numAtoms"), '=', '-');

        output << separatorStar;
        output << "Atom types and atom numbers\n";
        output << separatorDash;
        const Vec1String& types = s.cget<Vec1String>("typeNames");
        const Vec1Int&    numAtomsPerType = s.cget<Vec1Int>("numAtomsPerType");
        for (Size j = 0; j < types.size(); ++j)
        {
            output << str("%s %d\n", types[j].c_str(), numAtomsPerType[j]);
        }

        if (s.contains("CTIFOR"))
        {
            output << mlabItemString("CTIFOR", s.cget<Real>("CTIFOR"), '=', '-');
        }

        output << mlabItemString("Primitive lattice vectors (ang.)",
                                 s.cget<Vec1Real>("lattice"),
                                 '=',
                                 '-');

        output << mlabItemString("Atomic positions (ang.)",
                                 s.cget<Vec1Real>("positions"),
                                 '=',
                                 '-');

        output << mlabItemString("Total energy (eV)", s.cget<Real>("energy"), '=', '-');

        output << mlabItemString("Forces (eV ang.^-1)", s.cget<Vec1Real>("forces"), '=', '-');

        output << separatorEqual;
        output << "Stress (kbar)\n";
        const Vec1Real& stress = s.cget<Vec1Real>("stress");
        Vec1Real        stressPart(stress.begin(), stress.begin() + 3);
        output << mlabItemString("XX YY ZZ", stressPart, '-', '-');
        stressPart.clear();
        stressPart.insert(stressPart.end(), stress.begin() + 3, stress.end());
        output << mlabItemString("XY YZ ZX", stressPart, '-', '-');

        if (s.contains("charges"))
        {
            output << mlabItemString("Charges (e)", s.cget<Vec1Real>("charges"), '=', '-', 1);
        }
    }

    return;
}

void io::writeMlab(const Record& record, String fileName)
{
    std::ofstream file;

    io::open(file, fileName, "w");
    io::writeMlab(record, file);
    io::close(file);

    return;
}

void io::readXyz(Record& record, std::istream& input)
{
    if (!record.empty()) global::tutor.error("Cannot read XYZ file, record is not empty.");

    String line;
    Int    n = 0;

    record["structures"] = Vec1ShRec();
    Vec1ShRec& structures = record.get<Vec1ShRec>("structures");

    while (true)
    {
        Int m = readLine(input, line, 0, true);
        std::cout << "line1: " << line << "\n";
        // If EOF or empty line, stop reading.
        if (m == 0 || line.empty()) break;
        n += m;
        Int numAtoms = extxyz::readScalar<Int>(line);
        structures.push_back(std::make_shared<Record>());
        Record& s = *structures.back();
        s.put("numAtoms", numAtoms);

        n += readLine(input, line);
        std::cout << "line2: " << line << "\n";
        s.put("exyz", line);

        extxyz::parseCommentLine(s, line);
        //break;

        n += readLine(input, line, numAtoms - 1);
        std::cout << "line3: " << line << " " << n << "\n";
        break;
    }

    return;
}

void io::writeXyz(const Record& record, std::ostream& output)
{
    using namespace extxyz;

    const Vec1Real& atomicMass = record.cget<Vec1Real>("atomicMass");
    const Vec1Real& atomicRefEnergy = record.cget<Vec1Real>("atomicRefEnergy");

    for (Int i = 0; i < record.cget<Int>("numStructures"); ++i)
    {
        const Record& s = *record.vcget<Vec1ShRec>("structures")[i];

        Int numAtoms = s.cget<Int>("numAtoms");
        output << valueToString(numAtoms, "%d") << "\n";
        output << "Lattice=\"" << valueToString(s.cget<Vec1Real>("lattice"), "%13.6E", true)
               << "\" ";
        output << "Energy=" << valueToString(s.cget<Real>("energy"), "%24.16E", true) << " ";
        output << "Stress=\"" << valueToString(s.cget<Vec1Real>("stress"), "%13.6E", true) << "\" ";
        output << "Mass=\"" << valueToString(atomicMass, "%13.6E", true) << "\" ";
        output << "RefEnergy=\"" << valueToString(atomicRefEnergy, "%13.6E", true) << "\" ";
        if (s.contains("CTIFOR"))
        {
            output << "CTIFOR=" << valueToString(s.cget<Real>("CTIFOR"), "%13.6E", true) << " ";
        }
        Vec1Int lrcplus1 = s.cget<Vec1Int>("lrc");
        for (auto& ilrc : lrcplus1) ilrc++;
        output << "LocalRef=\"" << valueToString(lrcplus1, "%d", true) << "\" ";
        output << "Properties=species:S:1:pos:R:3:forces:R:3\n";
        const Vec1Real&   positions = s.cget<Vec1Real>("positions");
        const Vec1Real&   forces = s.cget<Vec1Real>("forces");
        const Vec1String& types = s.cget<Vec1String>("typeNames");
        const Size        maxLengthTypes = string_tools::maxLength(types);
        const String      fmtTypes(str("%%-%zus", maxLengthTypes));
        const Vec1Int&    numAtomsPerType = s.cget<Vec1Int>("numAtomsPerType");
        Size              ai = 0;
        for (Size j = 0; j < types.size(); ++j)
        {
            for (Int k = 0; k < numAtomsPerType[j]; ++k)
            {
                output << valueToString(types[j], fmtTypes) << " ";
                output << valueToString(positions[3 * ai], "%24.16E") << " "
                       << valueToString(positions[3 * ai + 1], "%24.16E") << " "
                       << valueToString(positions[3 * ai + 2], "%24.16E") << " ";
                output << valueToString(forces[3 * ai], "%24.16E") << " "
                       << valueToString(forces[3 * ai + 1], "%24.16E") << " "
                       << valueToString(forces[3 * ai + 2], "%24.16E") << "\n";
                ai++;
            }
        }
    }

    return;
}

void io::updateMlffHeaderDate(Record& rec)
{
    if (!rec.contains("header")) return;
    Record& hRec = rec.dget<ShRec>("header");
    if (hRec.contains("date")) hRec.erase("date");
    hRec.put<String>("date", now());

    return;
}

/*================================================================================================+
 | Guidelines for adding new items in ML_FF record:
 |
 | * Consider cases where the existence of some items (e.g. ML_LMAX2_RED) depend on the value of
 |   another (e.g. ML_DESC_TYPE): older M_FF version may not have the "deciding" item at all
 |   (ML_DESC_TYPE) but we should still add it with the default value (ML_DESC_TYPE = 0) to allow
 |   writing the old force field with the newest file format. Examples:
 |   - ML_DESC_TYPE
 |   - ML_BASIS_TYPE2
 |
 | * Similarly, some items were only added explicitly in a certain version. Before that, the value
 |   was maybe deduced from some other item. These items should also be added even when reading in
 |   an old file, e.g.:
 |   - ML_LFAST
 |   - WRITECMAT
 +================================================================================================*/
void io::processMlff(Record& rec, std::fstream& strm, bool doWrite)
{
    Buffer buf;
    Size   pos = 0;

    if (!doWrite && !rec.empty())
    {
        std::logic_error("Attempt to read in ML_FF file but input record is not empty.");
    }
    if (doWrite && rec.empty())
    {
        std::logic_error("Attempt to write ML_FF file but input record is empty.");
    }

    if (!doWrite) readBuffer(strm, buf);

    processMlffHeader(buf, pos, doWrite, rec);
    processMlffEntry<Int>(buf, pos, doWrite, rec, "versionMajor");
    processMlffEntry<Int>(buf, pos, doWrite, rec, "versionMinor");
    processMlffEntry<Int>(buf, pos, doWrite, rec, "versionPatch");
    SemanticVersion versionIn(rec.cget<Int>("versionMajor"),
                              rec.cget<Int>("versionMinor"),
                              rec.cget<Int>("versionPatch"));

    if (!doWrite) rec.put<Int>("ML_DESC_TYPE", 0);
    if (versionIn >= SemanticVersion(0, 2, 1))
    {
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_DESC_TYPE", -1, true);
    }
    const Int descriptorType = rec.cget<Int>("ML_DESC_TYPE");

    processMlffEntry<Int>(buf, pos, doWrite, rec, "numTypes");
    const Int numTypes = rec.cget<Int>("numTypes");

    if (!doWrite) rec.put<Vec1String>("types", Vec1String(numTypes, "??"));
    processMlffEntry<Vec1String>(buf, pos, doWrite, rec, "types", 2, true);
    processMlffEntry<Vec1Real>(buf, pos, doWrite, rec, "ML_EATOM_REF", 0, false, numTypes);
    processMlffEntry<Vec1Real>(buf, pos, doWrite, rec, "atomicMass", 0, false, numTypes);
    processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IWEIGHT");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_WTOTEN");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_WTIFOR");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_WTSIF");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "rmseEnergy");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "rmseForces");
    processMlffEntry<Real>(buf, pos, doWrite, rec, "rmseStress");

    processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LMLMB");

    if (rec.cget<bool>("ML_LMLMB")) // Should always be true!
    {
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LSIC");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LSUPERVEC");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_W1");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_W2");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_ICUT1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_ICUT2");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_RCUT1");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_RCUT2");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IBROAD1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IBROAD2");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_SION1");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_SION2");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MRB1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MRB2");
        if (!doWrite) rec.put<Int>("ML_BASIS_TYPE2", 0);
        if (versionIn >= SemanticVersion(0, 2, 2))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_BASIS_TYPE2", -1, true);
            const Int basisType = rec.cget<Int>("ML_BASIS_TYPE2");
            if (basisType == 4 || basisType == 5)
            {
                processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MRB2_MAX");
            }
            if (descriptorType == 2 || descriptorType == 4 || descriptorType == 7
                || descriptorType == 9)
            {
                processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MRB2_RED");
                processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_BASIS_TYPE2_RED");
                const Int basisTypeRed = rec.cget<Int>("ML_BASIS_TYPE2_RED");
                if (basisTypeRed == 4 || basisTypeRed == 5)
                {
                    processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MRB2_MAX_RED");
                }
            }
        }
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NR1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NR2");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MSPL1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_MSPL2");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_LMAX1");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_LMAX2");
        if (descriptorType == 2 || descriptorType == 4 || descriptorType == 7
            || descriptorType == 9)
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_LMAX2_RED");
            processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_DESC_RATIO_DUAL");
            processMlffEntry<Real>(buf, pos, doWrite, rec, "RATIO_DESC_LOW_RES");
            processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_DESC_FACTORE_TESTE");
        }
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NHYP");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NHYP2");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LNORM1");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LNORM2");

        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LWINDOW1");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LWINDOW2");
        if (rec.cget<bool>("ML_LWINDOW1"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IWINDOW1");
        }
        if (rec.cget<bool>("ML_LWINDOW2"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IWINDOW2");
        }

        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LAFILT2");
        if (rec.cget<bool>("ML_LAFILT2"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_IAFILT2");
            if (rec.cget<Int>("ML_IAFILT2") == 2)
            {
                processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_AFILT2");
            }
        }

        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LMETRIC1");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LMETRIC2");
        if (rec.cget<bool>("ML_LMETRIC1"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NMETRIC1");
            processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_RMETRIC1");
        }
        if (rec.cget<bool>("ML_LMETRIC2"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NMETRIC2");
            processMlffEntry<Real>(buf, pos, doWrite, rec, "ML_RMETRIC2");
        }

        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LVARTRAN1");
        processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LVARTRAN2");
        if (rec.cget<bool>("ML_LVARTRAN1"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NVARTRAN1");
        }
        if (rec.cget<bool>("ML_LVARTRAN2"))
        {
            processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_NVARTRAN2");
        }

        if (rec.cget<Real>("ML_W2") > 0.0)
        {
            processMlffEntry<Vec1Int>(buf,
                                      pos,
                                      doWrite,
                                      rec,
                                      "SHS3-3-body-number-radial-basis",
                                      0,
                                      false,
                                      rec.cget<Int>("ML_LMAX2") + 1);
            if (descriptorType == 2 || descriptorType == 4 || descriptorType == 7
                || descriptorType == 9)
            {
                processMlffEntry<Vec1Int>(buf,
                                          pos,
                                          doWrite,
                                          rec,
                                          "SHS3-3-body-reduced-number-radial-basis",
                                          0,
                                          false,
                                          rec.cget<Int>("ML_LMAX2_RED") + 1);
            }
            if (!doWrite)
            {
                rec.put<Vec1Int>("SHS3-3-body-number-descriptors-per-type", Vec1Int(numTypes));
                rec.put<Vec2Int>("SHS3-3-body-descriptor-list", Vec2Int(numTypes));
                if (descriptorType == 2 || descriptorType == 4 || descriptorType == 7
                    || descriptorType == 9)
                {
                    rec.put<Vec1Int>("ndesc_per_type2_dual", Vec1Int(numTypes));
                    rec.put<Vec2Int>("descriptor_list2_dual", Vec2Int(numTypes));
                }
            }
            for (Int itype = 0; itype < numTypes; ++itype)
            {
                Int n;
                if (doWrite)
                {
                    n = rec.get<Vec1Int>("SHS3-3-body-number-descriptors-per-type")[itype];
                    write(n, buf, 0);
                }
                else
                {
                    read(buf, n, pos);
                    rec.get<Vec1Int>("SHS3-3-body-number-descriptors-per-type")[itype] = n;
                    rec.get<Vec2Int>("SHS3-3-body-descriptor-list")[itype].resize(n);
                }
                process(buf,
                        rec.get<Vec2Int>("SHS3-3-body-descriptor-list")[itype],
                        pos,
                        doWrite,
                        0);
                if (descriptorType == 2 || descriptorType == 4 || descriptorType == 7
                    || descriptorType == 9)
                {
                    if (doWrite)
                    {
                        n = rec.get<Vec1Int>("ndesc_per_type2_dual")[itype];
                        write(n, buf, 0);
                    }
                    else
                    {
                        read(buf, n, pos);
                        rec.get<Vec1Int>("ndesc_per_type2_dual")[itype] = n;
                        rec.get<Vec2Int>("descriptor_list2_dual")[itype].resize(n);
                    }
                    process(buf, rec.get<Vec2Int>("descriptor_list2_dual")[itype], pos, doWrite, 0);
                }
            }
            if (rec.cget<bool>("ML_LSIC"))
            {
                if (!doWrite)
                {
                    rec.put<Vec1Int>("num-sic-active-descriptors", Vec1Int(numTypes));
                    rec.put<Vec2Int>("sic-active-descriptors", Vec2Int(numTypes));
                }
                for (Int itype = 0; itype < numTypes; ++itype)
                {
                    Int n;
                    if (doWrite)
                    {
                        n = rec.get<Vec1Int>("num-sic-active-descriptors")[itype];
                        write(n, buf, 0);
                    }
                    else
                    {
                        read(buf, n, pos);
                        rec.get<Vec1Int>("num-sic-active-descriptors")[itype] = n;
                        rec.get<Vec2Int>("sic-active-descriptors")[itype].resize(n);
                    }
                    process(buf,
                            rec.get<Vec2Int>("sic-active-descriptors")[itype],
                            pos,
                            doWrite,
                            0);
                }
            }
        } // ML_W2

        processMlffEntry<Vec1Int>(buf, pos, doWrite, rec, "numLrc", 0, false, numTypes);
        const Vec1Int& numLrc = rec.cget<Vec1Int>("numLrc");
        processMlffEntry<Int>(buf, pos, doWrite, rec, "ML_ISCALE_TOTEN");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "average-energy-per-atom");

        if (!doWrite) rec.put<Vec2Real>("regression-coeff", Vec2Real(numTypes));
        for (Int itype = 0; itype < numTypes; ++itype)
        {
            Vec1Real& coeff = rec.get<Vec2Real>("regression-coeff")[itype];
            if (!doWrite) coeff.resize(numLrc[itype]);
            process(buf, coeff, pos, doWrite, 0);
        }

        if (rec.cget<Real>("ML_W1") > 0.0)
        {
            if (!doWrite)
            {
                rec.put<Int>("SHS2-2-body-number-descriptors-per-type",
                             rec.cget<Int>("ML_MRB1") * numTypes);
            }
            if (!doWrite) rec.put<Vec2Real>("SHS2-2-body-reference-configs", Vec2Real(numTypes));
            for (Int itype = 0; itype < numTypes; ++itype)
            {
                Vec1Real& refConf = rec.get<Vec2Real>("SHS2-2-body-reference-configs")[itype];
                if (!doWrite)
                {
                    refConf.resize(numLrc[itype]
                                   * rec.cget<Int>("SHS2-2-body-number-descriptors-per-type"));
                }
                process(buf, refConf, pos, doWrite, 0);
            }
        } // ML_W1

        if (rec.cget<Real>("ML_W2") > 0.0)
        {
            if (!doWrite) rec.put<Vec2Real>("SHS3-3-body-reference-configs", Vec2Real(numTypes));
            for (Int itype = 0; itype < numTypes; ++itype)
            {
                Vec1Real& refConf = rec.get<Vec2Real>("SHS3-3-body-reference-configs")[itype];
                if (!doWrite)
                {
                    refConf.resize(
                        numLrc[itype]
                        * rec.cget<Vec1Int>("SHS3-3-body-number-descriptors-per-type")[itype]);
                }
                process(buf, refConf, pos, doWrite, 0);
            }

            if (descriptorType == 2 || descriptorType == 4)
            {
                if (!doWrite) rec.put<Vec2Real>("ang_desc_dual", Vec2Real(numTypes));
                for (Int itype = 0; itype < numTypes; ++itype)
                {
                    Vec1Real& refConf = rec.get<Vec2Real>("ang_desc_dual")[itype];
                    if (!doWrite)
                    {
                        refConf.resize(numLrc[itype]
                                       * rec.cget<Vec1Int>("ndesc_per_type2_dual")[itype]);
                    }
                    process(buf, refConf, pos, doWrite, 0);
                }
            }
        } // ML_W2

        processMlffEntry<Real>(buf, pos, doWrite, rec, "sigv");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "sigw");
        processMlffEntry<Real>(buf, pos, doWrite, rec, "variance-training-energies");
        processMlffEntry<Vec1Real>(buf, pos, doWrite, rec, "variance-training-force", 0, false, 3);
        processMlffEntry<Vec1Real>(buf, pos, doWrite, rec, "variance-training-stress", 0, false, 6);

        if (!doWrite)
        {
            rec.put<bool>("ML_LFAST", false);
            if (versionIn < SemanticVersion(0, 2, 0)) { rec.get<bool>("ML_LFAST") = false; }
            else if (versionIn < SemanticVersion(0, 2, 4))
            {
                rec.get<bool>("ML_LFAST") = rec.dcget<ShRec>("header").cget<bool>("ML_LFAST");
            }
        }
        if (versionIn >= SemanticVersion(0, 2, 4))
        {
            processMlffEntry<bool>(buf, pos, doWrite, rec, "ML_LFAST", -1, true);
        }

        if (!doWrite)
        {
            rec.put<bool>("WRITECMAT", true);
            if (versionIn < SemanticVersion(0, 2, 3)) { rec.get<bool>("WRITECMAT") = true; }
            else if (versionIn < SemanticVersion(0, 2, 5))
            {
                rec.get<bool>("WRITECMAT") = !rec.cget<bool>("ML_LFAST");
            }
        }
        if (versionIn >= SemanticVersion(0, 2, 5))
        {
            processMlffEntry<bool>(buf, pos, doWrite, rec, "WRITECMAT", -1, true);
        }

        if (rec.cget<bool>("WRITECMAT"))
        {
            Size numLrcTotal = 0;
            for (Real ilrc : numLrc) numLrcTotal += ilrc;
            processMlffEntry<Vec1Real>(buf,
                                       pos,
                                       doWrite,
                                       rec,
                                       "covariance-matrix",
                                       0,
                                       false,
                                       numLrcTotal * numLrcTotal);
        }
        else
        {
            if (!doWrite) rec.put<Vec2Real>("inverse-kernel-matrix", Vec2Real(numTypes));
            for (Int itype = 0; itype < numTypes; ++itype)
            {
                Vec1Real& cov = rec.get<Vec2Real>("inverse-kernel-matrix")[itype];
                cov.resize(numLrc[itype] * numLrc[itype]);
                process(buf, cov, pos, doWrite, 0);
            }
        }
    } // ML_LMLMB

    if (doWrite)
    {
        strm.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(Byte));
    }

    return;
}

void io::readMlffAndConvertUnits(Record& record, String fileName)
{
    std::fstream file;

    io::open(file, fileName, "rb");
    io::processMlff(record, file, false);
    io::close(file);

    setup::rescaleMlffUnits(record,
                            1.0 / constants::EUNIT,
                            1.0 / constants::AUTOA,
                            1.0 / constants::FUNIT,
                            1.0 / constants::SUNIT);

    return;
}

void io::readIncar(Record& record, std::istream& strm)
{
    io::detail::setupKnownTags(record);

    String input = toString(strm);
    record.put<String>("raw", input);
    if (input.empty())
    {
        record.get<String>("message").append("INCAR is empty, no tags were read in.");
    }
    // Add newline at end of file if not yet present (avoid crashes if line continuation is used at
    // the end of INCAR.
    if (input.back() != '\n') input.push_back('\n');

    io::detail::readIncarToStrings(record, input);
    io::detail::processIncarTags(record);

    return;
}

void io::setupFromIncar(const Record& incar, Record& setup)
{
    const Int maxSeed0 = 900000000;
    const Int maxSeed1 = 1000000;

    // Defaults from ml_reader.F
    Int  default_ISTART = 0;
    bool default_LFAST = false;
    Int  default_ICRITERIA = 1;
    Real default_CDOUB = 2.0;
    Int  default_OUTBLOCK = 1;
    Int  default_IALGO_LINREG = 1;
    Real default_EPS_LOW = 1.0E-9;

    // From ml_ff_struct.F
    Int default_ML_MB = 1500;
    Int default_ML_MCONF = 1500;

    // Pre-existing keys in setup record.
    const bool& mlabPresent = setup.cget<bool>("LMLABEXIST");
    const bool& mlffPresent = setup.cget<bool>("LMLFFEXIST");
    Int&  nsw = setup.get<Int>("NSW");
    const Int&  ntyp = setup.cget<Int>("NTYP");
    const Int&  ibrion = setup.cget<Int>("IBRION");

    const Record& types = incar.dcget<ShRec>("types");
    Record&       states = setup.add("states", ItemIndex::SHREC).dget<ShRec>("states");
    for (String tag : types.keys()) states.put<Int>(tag, static_cast<Int>(TagState::UNSET));

    setTagValue<String>(setup, incar, "ML_MODE", "none");
    const String& mode = setup.cget<String>("ML_MODE");
    if (mode == "train")
    {
        if (mlabPresent) default_ISTART = 1;
        else default_ISTART = 0;
    }
    else if (mode == "refit")
    {
        default_ISTART = 4;
        nsw = 1;
        default_LFAST = true;
        default_IALGO_LINREG = 4;
        default_EPS_LOW = 1.0E-11;
    }
    else if (mode == "refitbayesian")
    {
        default_ISTART = 4;
        nsw = 1;
        default_EPS_LOW = 1.0E-11;
    }
    else if (mode == "select")
    {
        default_ISTART = 3;
        nsw = 1;
        default_ICRITERIA = 3;
        default_CDOUB = 100.0;
    }
    else if (mode == "run")
    {
        default_ISTART = 2;
        default_OUTBLOCK = 10;
    }
    else if (mode == "delta")
    {
        default_ISTART = 2;
    }
    else if (mode == "none")
    {
        // TODO: Maybe this should not be ISTART = 0 but -1?
        default_ISTART = 0;
    }
    //else if (mode == "ranse")
    //{
    //    default_ISTART = 3;
    //    nsw = 1;
    //    default_ICRITERIA = 3;
    //    default_CDOUB = 100.0;
    //}
    //else if (mode == "fpsse")
    //{
    //    default_ISTART = 3;
    //    nsw = 1;
    //    default_ICRITERIA = 3;
    //    default_CDOUB = 100.0;
    //}
    else
    {
        global::tutor.error("Option for ML_MODE not recognized. Please specify a valid option.");
    }

    setTagValue<String>(setup, incar, "ML_TYPE", "kernel");
    const String& type = setup.cget<String>("ML_TYPE");
    if (type == "vasp" || type == "kernel") setup.get<String>("ML_TYPE") = "kernel";
    else if (type == "grace")
    {
#ifndef VASPML_ENABLE_GRACE
        global::tutor.error(
            "Found ML_TYPE = grace in INCAR but VASPml was compiled without GRACE force field "
            "support. Please recompile VASP with GRACE support (-DVASPML_ENABLE_GRACE).");
#endif
    }
    setTagValue<String>(setup, incar, "ML_GRACE_MODEL", "./SavedModel");
    setTagValue<bool>(setup, incar, "ML_LIB", true);
    setTagValue<String>(setup, incar, "ML_PALGO", "off");
    setTagValue<Int>(setup, incar, "ML_NCSHMEM", -1);
    setTagValue<Int>(setup, incar, "ML_CALGO", 0);
    const Int& calgo = setup.cget<Int>("ML_CALGO");
    if (calgo < 0 || calgo > 1)
    {
        global::tutor.error(
            "ERROR: Option for ML_CALGO not recognized. Please specify a valid option.");
    }
    setTagValue<Real>(setup, incar, "ML_CDOUB", default_CDOUB);
    setTagValue<Real>(setup, incar, "ML_CSIG", 4E-01);
    setTagValue<Real>(setup, incar, "ML_CSLOPE", 2E-01);
    if (setup.cget<Int>("ML_CALGO") == 0) setTagValue<Real>(setup, incar, "ML_CTIFOR", 0.002);
    else setTagValue<Real>(setup, incar, "ML_CTIFOR", 0.02);
    setTagValue<Real>(setup, incar, "ML_SCLC_CTIFOR", 0.6);
    setTagValue<Vec1Real>(setup, incar, "ML_EATOM_REF", Vec1Real(ntyp, 0.0));
    setTagValue<Real>(setup, incar, "ML_EPS_LOW", default_EPS_LOW);
    setTagValue<Real>(setup, incar, "ML_EPS_REG", 1E-15);
    setTagValue<Int>(setup, incar, "ML_IMAT_SPARS", 1);
    setTagValue<Int>(setup, incar, "ML_ISTART", default_ISTART);
    const Int& istart = setup.cget<Int>("ML_ISTART");
    if (istart == 3)
    {
        global::tutor.warning(
            "ML_MODE = select was selected. Structure in POSCAR is ignored, hence no prediction of "
            "error, forces and stress will be available on that structure at the end of this "
            "calculation.");
    }
    if (istart == 2 && type == "kernel" && !mlffPresent)
    {
        global::tutor.error("ERROR: ML_FF is missing, exiting...");
    }
    if ((istart == 1 || istart == 3 || istart == 4) && !mlabPresent)
    {
        global::tutor.error("ERROR: ML_AB is missing, exiting...");
    }
    // Note: LMLONLY is not used here.
    setTagValue<Int>(setup, incar, "ML_ISVD", 1);
    setTagValue<Int>(setup, incar, "ML_IUPDATE_CRITERIA", 1);
    setTagValue<bool>(setup, incar, "ML_LBASIS_DISCARD", true);
    setTagValue<bool>(setup, incar, "ML_LCONF_DISCARD", false);
    setTagValue<bool>(setup, incar, "ML_LUSE_NAMES", false);
    if (calgo == 0) setTagValue<Int>(setup, incar, "ML_ICRITERIA", default_ICRITERIA);
    else if (calgo == 1) setTagValue<Int>(setup, incar, "ML_ICRITERIA", 0);
    if (calgo == 1 && setup.cget<Int>("ML_ICRITERIA") > 0)
    {
        global::tutor.error("ERROR: For ML_CALGO = 1, ML_ICRITERIA = 0 must be set.");
    }
    setTagValue<bool>(setup, incar, "ML_LMLMB", true);
    setTagValue<Int>(setup, incar, "ML_NDIM_SCALAPACK", 2);
    setTagValue<Int>(setup, incar, "ML_FF_NWRITE", 2);
    setTagValue<Int>(setup, incar, "ML_IWEIGHT", 3);
    setTagValue<bool>(setup, incar, "ML_LNMDINT_RANDOM", false);
    setTagValue<bool>(setup, incar, "ML_LTEST", false);
    setTagValue<bool>(setup, incar, "ML_LTRJ", false);
    setTagValue<bool>(setup, incar, "ML_LFAST", default_LFAST);
    const bool& lfast = setup.cget<bool>("ML_LFAST");
    if (lfast && istart == 1 && nsw > 1)
    {
        global::tutor.error("ERROR: ML_LFAST only available for ML_ISTART = 1 and NSW = 0,1.");
    }
    if (lfast && (istart == 0 || istart == 3))
    {
        global::tutor.error(
            "ERROR: ML_LFAST not available for ML_ISTART = " + std::to_string(istart) + ".");
    }
    setTagValue<bool>(setup, incar, "ML_LFORCESLOW", false);
    setTagValue<Int>(setup, incar, "ML_IALGO_LINREG", default_IALGO_LINREG);
    const Int& ialgolinreg = setup.cget<Int>("ML_IALGO_LINREG");
    if (lfast && ialgolinreg < 3 && istart != 2)
    {
        global::tutor.error("ERROR: ML_LFAST only available for ML_IALGO_LINREG > 2");
    }
    setTagValue<bool>(setup, incar, "ML_LTOTEN_SYSTEM", false);
    setTagValue<Int>(setup, incar, "ML_MCONF", default_ML_MCONF);
    setTagValue<Int>(setup, incar, "ML_MCONF_NEW", 5);
    setTagValue<Int>(setup, incar, "ML_MHIS", 10);

    if (istart == 3) setTagValue<Int>(setup, incar, "ML_NMDINT", 1);
    else setTagValue<Int>(setup, incar, "ML_NMDINT", 10);
    const bool& ltest = setup.cget<bool>("ML_LTEST");
    if (ltest) setTagValue<Int>(setup, incar, "ML_NTEST", 10);
    if (istart == 2) setTagValue<Int>(setup, incar, "ML_OUTPUT_MODE", 0);
    else setTagValue<Int>(setup, incar, "ML_OUTPUT_MODE", 1);
    setTagValue<Int>(setup, incar, "ML_OUTBLOCK", default_OUTBLOCK);
    Int& outblock = setup.get<Int>("ML_OUTBLOCK");
    if (outblock < 1)
    {
        global::tutor.error("ERROR: ML_OUTBLOCK must be larger than or equal to 1.");
    }
    if (outblock > std::max(nsw, 1))
    {
        if (states.cget<Int>("ML_OUTBLOCK") != static_cast<Int>(TagState::DEFAULT))
        {
            global::tutor.warning("ML_OUTBLOCK > max(NSW,1), and will be reset to ML_OUTBLOCK = "
                                  + std::to_string(std::max(nsw, 1)) + ".");
            states.get<Int>("ML_OUTBLOCK") = static_cast<Int>(TagState::OVERRIDE);
        }
        outblock = std::max(nsw, 1);
    }
    if ((istart != 2 || mode == "delta") && outblock != 1)
    {
        global::tutor.warning("ML_OUTBLOCK > 1 is only supported for ML_MODE = run, "
                              "defaulting back to ML_OUTBLOCK = 1 ...");
        outblock = 1;
    }
    if (outblock > 1 && ibrion != 0)
    {
        if (states.cget<Int>("ML_OUTBLOCK") != static_cast<Int>(TagState::DEFAULT))
        {
            global::tutor.warning("ML_OUTBLOCK > 1 is only supported for IBRION = 0, defaulting "
                                  "back to ML_OUTBLOCK = 1 ...");
        }
        outblock = 1;
    }
    if (istart == 2 && mode != "delta") setTagValue<Int>(setup, incar, "ML_ESTBLOCK", outblock);
    else setTagValue<Int>(setup, incar, "ML_ESTBLOCK", 1);
    const Int& estblock = setup.cget<Int>("ML_ESTBLOCK");
    if (istart != 2 && estblock != 1)
    {
        global::tutor.error("ERROR: ML_ESTBLOCK=1 must be set for ML_MODE/=run.");
    }
    if (istart == 2 && estblock > 0)
    {
        if (estblock > outblock)
        {
            global::tutor.error("ML_ESTBLOCK is smaller than ML_OUTBLOCK. ML_ESTBLOCK needs "
                                "to be set as an integer multiple of ML_OUTBLOCK.");
        }
        else if (estblock % outblock != 0)
        {
            global::tutor.error(
                "ERROR: ML_ESTBLOCK must be integer divisible by ML_OUTBLOCK. ML_ESTBLOCK needs to "
                "be set as an integer multiple of ML_OUTBLOCK.");
        }
    }
    if (setup.cget<bool>("ML_LFORCESLOW") && estblock > 0)
    {
        global::tutor.error("ERROR: ML_ESTBLOCK=0 must be set for ML_LFORCESLOW = .TRUE.");
    }
    setTagValue<Real>(setup, incar, "ML_SF_LIMIT", 0.9);
    if (setup.cget<Int>("ML_IWEIGHT") == 1)
    {
        setTagValue<Real>(setup, incar, "ML_WTIFOR", 5E-02);
        setTagValue<Real>(setup, incar, "ML_WTOTEN", 5E-03);
        setTagValue<Real>(setup, incar, "ML_WTSIF", 5.0);
    }
    else
    {
        setTagValue<Real>(setup, incar, "ML_WTIFOR", 1.0);
        setTagValue<Real>(setup, incar, "ML_WTOTEN", 1.0);
        setTagValue<Real>(setup, incar, "ML_WTSIF", 1.0);
    }
    setTagValue<Real>(setup, incar, "ML_CX", 0.0);
    if (!setup.cget<bool>("ML_LMLMB")) return;
    setup.put<bool>("ML_LSIC", !setup.cget<bool>("ML_LFAST"));
    states.get<Int>("ML_LSIC") = static_cast<Int>(TagState::DEFAULT);
    setTagValue<Real>(setup, incar, "ML_W1", 0.1);
    setup.put<Real>("ML_W2", 1.0 - setup.get<Real>("ML_W1"));
    const Real& w2 = setup.cget<Real>("ML_W2");
    setTagValue<Int>(setup, incar, "ML_DESC_TYPE", 0);
    const Int& desctype = setup.cget<Int>("ML_DESC_TYPE");
    if (desctype > 0) setup.get<bool>("ML_LSIC") = false;
    if (desctype > 0 && w2 <= 0.0)
    {
        global::tutor.error("ERROR: ML_DESC_TYPE must be 0 for ML_W2 = 0.0.");
    }
    if (desctype == 2 || desctype == 4)
    {
        if (mode == "train" || mode == "select")
        {
            global::tutor.error("ERROR: ML_DESC_TYPE = 2,4 not available for ML_MODE = " + mode
                                + ".");
        }
        if (istart == 0 || istart == 1 || istart == 3)
        {
            global::tutor.error("ERROR: ML_DESC_TYPE = 2,4 not available for ML_ISTART = "
                                + std::to_string(istart) + ".");
        }
        if (mode == "refit")
        {
            global::tutor.error("ERROR: ML_DESC_TYPE = 2,4 not available for ML_MODE = " + mode
                                + ". Please select ML_MODE = refitbayesian.");
        }
    }
    setTagValue<Real>(setup, incar, "ML_DESC_RATIO_DUAL", 0.5);
    setTagValue<Real>(setup, incar, "ML_DESC_FACTORE_TESTE", 1.0);
    setTagValue<bool>(setup, incar, "ML_LAFILT2", true);
    const bool& lafilt2 = setup.cget<bool>("ML_LAFILT2");
    if (lafilt2)
    {
        setTagValue<Int>(setup, incar, "ML_IAFILT2", 2);
        if (setup.cget<Int>("ML_IAFILT2") == 2)
        {
            setTagValue<Real>(setup, incar, "ML_AFILT2", 2.0E-03);
        }
    }
    setTagValue<Int>(setup, incar, "ML_IBROAD1", 2);
    setTagValue<Int>(setup, incar, "ML_IBROAD2", 2);
    setTagValue<Int>(setup, incar, "ML_BASIS_TYPE2", 0);
    setTagValue<Int>(setup, incar, "ML_BASIS_TYPE2_RED", 4);
    setTagValue<Int>(setup, incar, "ML_ICUT1", 1);
    setTagValue<Int>(setup, incar, "ML_ICUT2", 1);
    setTagValue<Int>(setup, incar, "ML_INVERSE_SOAP", 2);
    setTagValue<Int>(setup, incar, "ML_IREG", 2);
    setTagValue<bool>(setup, incar, "ML_SAVECMAT", true);
    setTagValue<Int>(setup, incar, "ML_ISCALE_TOTEN", 2);
    setTagValue<bool>(setup, incar, "ML_LSUPERVEC", true);
    setTagValue<bool>(setup, incar, "ML_LEATOM", false);
    setTagValue<bool>(setup, incar, "ML_LHEAT", false);
    setTagValue<bool>(setup, incar, "ML_LMETRIC1", false);
    setTagValue<bool>(setup, incar, "ML_LMETRIC2", false);
    setTagValue<bool>(setup, incar, "ML_LSPARSDES", false);
    if (setup.cget<bool>("ML_LSPARSDES") && (mode == "train" || mode == "select"))
    {
        global::tutor.error("ML_LSPARSDES = .TRUE. is not allowed together with ML_MODE = "
                            "train or ML_MODE = select.");
    }
    setTagValue<Int>(setup, incar, "ML_NRANK_SPARSDES", 5);
    setTagValue<Real>(setup, incar, "ML_RDES_SPARSDES", 0.5);
    if (istart != 1 && istart != 4
        && states.cget<Int>("ML_LSPARSDES") != static_cast<Int>(TagState::DEFAULT))
    {
        setup.get<bool>("ML_LSPARSDES") = false;
        states.get<Int>("ML_LSPARSDES") = static_cast<Int>(TagState::OVERRIDE);
    }
    setTagValue<bool>(setup, incar, "ML_LVARTRAN1", false);
    setTagValue<bool>(setup, incar, "ML_LVARTRAN2", false);
    if (setup.cget<bool>("ML_LMETRIC1")) setTagValue<Int>(setup, incar, "ML_NMETRIC1", 6);
    if (setup.cget<bool>("ML_LMETRIC2")) setTagValue<Int>(setup, incar, "ML_NMETRIC2", 6);
    if (setup.cget<bool>("ML_LVARTRAN1")) setTagValue<Int>(setup, incar, "ML_NVARTRAN1", 6);
    if (setup.cget<bool>("ML_LVARTRAN2")) setTagValue<Int>(setup, incar, "ML_NVARTRAN2", 6);
    setTagValue<bool>(setup, incar, "ML_LWINDOW1", false);
    setTagValue<bool>(setup, incar, "ML_LWINDOW2", false);
    if (setup.cget<bool>("ML_LWINDOW1")) setTagValue<Int>(setup, incar, "ML_IWINDOW1", 6);
    if (setup.cget<bool>("ML_LWINDOW2")) setTagValue<Int>(setup, incar, "ML_IWINDOW2", 6);
    if (lafilt2) setTagValue<Int>(setup, incar, "ML_LMAX2", 3);
    else setTagValue<Int>(setup, incar, "ML_LMAX2", 6);
    setTagValue<Int>(setup, incar, "ML_LMAX2_RED", 3);
    setTagValue<bool>(setup, incar, "ML_LNORM1", true);
    setTagValue<bool>(setup, incar, "ML_LNORM2", true);
    setTagValue<Int>(setup, incar, "ML_MB", default_ML_MB);
    setTagValue<Int>(setup, incar, "ML_MB_MIN", 3);
    if (setup.cget<Int>("ML_MB_MIN") < 2)
    {
        global::tutor.error(
            "The minimum number of local reference configurations (ML_MB_MIN) must be at least 2.");
    }
    setTagValue<Int>(setup, incar, "ML_MSPL1", 1000);
    setTagValue<Int>(setup, incar, "ML_MSPL2", setup.cget<Int>("ML_MSPL1"));
    const bool& lsupervec = setup.cget<bool>("ML_LSUPERVEC");
    if (lsupervec) setTagValue<Int>(setup, incar, "ML_NHYP", 4);
    else setTagValue<Int>(setup, incar, "ML_NHYP", 1);
    setTagValue<Int>(setup, incar, "ML_NHYP2", 4);
    if (lsupervec) setup.get<Int>("ML_NHYP2") = setup.cget<Int>("ML_NHYP");
    setTagValue<Int>(setup, incar, "ML_NR1", 1000);
    setTagValue<Int>(setup, incar, "ML_NR2", setup.cget<Int>("ML_NR1"));
    const Real& w1 = setup.cget<Real>("ML_W1");
    setTagValue<Real>(setup, incar, "ML_RCUT1", 8.0);
    setTagValue<Real>(setup, incar, "ML_RCUT2", 5.0);
    if (w1 <= 0.0) setup.get<Real>("ML_RCUT1") = setup.cget<Real>("ML_RCUT2");
    if (w2 <= 0.0) setup.get<Real>("ML_RCUT2") = setup.cget<Real>("ML_RCUT1");
    if (setup.cget<bool>("ML_LMETRIC1")) setTagValue<Real>(setup, incar, "ML_RMETRIC1", 1.0);
    if (setup.cget<bool>("ML_LMETRIC2")) setTagValue<Real>(setup, incar, "ML_RMETRIC2", 1.0);
    setTagValue<Real>(setup, incar, "ML_SIGV0", 1.0);
    if (ialgolinreg == 4) setTagValue<Real>(setup, incar, "ML_SIGW0", 1E-7);
    else setTagValue<Real>(setup, incar, "ML_SIGW0", 1.0);
    if (setup.cget<Int>("ML_IBROAD1") != 1)
    {
        setTagValue<Real>(setup, incar, "ML_SION1", 5.0E-01);
    }
    if (setup.cget<Int>("ML_IBROAD2") != 1)
    {
        setTagValue<Real>(setup, incar, "ML_SION2", setup.cget<Real>("ML_SION1"));
    }
    setTagValue<bool>(setup, incar, "ML_LEMPPOT", false);
    const bool& lemppot = setup.cget<bool>("ML_LEMPPOT");
    setTagValue<bool>(setup, incar, "ML_LCOUPLE", lemppot);
    const bool& lcouple = setup.cget<bool>("ML_LCOUPLE");
    if (lemppot && !lcouple)
    {
        global::tutor.error("ERROR: ML_LCOUPLE must be .TRUE. if ML_LEMPPOT=.TRUE.");
    }
    setTagValue<Int>(setup, incar, "ML_NATOM_COUPLED", 0);
    setTagValue<Vec1Int>(setup, incar, "ML_ICOUPLE", Vec1Int{});
    if (lcouple
        && setup.cget<Int>("ML_NATOM_COUPLED")
               != static_cast<Int>(setup.cget<Vec1Int>("ML_ICOUPLE").size()))
    {
        global::tutor.error("ERROR: Number of atoms given via ML_ICOUPLE is not equal to "
                            "value of ML_NATOM_COUPLED.");
    }
    if (lemppot) setTagValue<Real>(setup, incar, "ML_RCOUPLE", 0.0);
    else setTagValue<Real>(setup, incar, "ML_RCOUPLE", 1.0);
    setTagValue<Int>(setup, incar, "ML_MRB1", 12);
    setTagValue<Int>(setup, incar, "ML_MRB2", 8);
    setTagValue<Int>(setup, incar, "ML_MRB2_MAX", std::max(8, setup.cget<Int>("ML_MRB2")));
    if (setup.cget<Int>("ML_MRB2_MAX") < setup.cget<Int>("ML_MRB2"))
    {
        global::tutor.error(
            "ERROR: ML_MRB2_MAX is smaller than ML_MRB2. Please set ML_MRB2_MAX >= ML_MRB2.");
    }
    setTagValue<Int>(setup, incar, "ML_MRB2_RED", 2);
    setTagValue<Int>(setup, incar, "ML_MRB2_MAX_RED", std::max(8, setup.cget<Int>("ML_MRB2_RED")));
    if (setup.cget<Int>("ML_MRB2_MAX_RED") < setup.cget<Int>("ML_MRB2_RED"))
    {
        global::tutor.error("ERROR: ML_MRB2_MAX_RED is smaller than ML_MRB2_RED. Please set "
                            "ML_MRB2_MAX_RED>=ML_MRB2_RED.");
    }
    setTagValue<Vec1Int>(setup, incar, "ML_RANDOM_SEED", Vec1Int{});
    Vec1Int& seed = setup.get<Vec1Int>("ML_RANDOM_SEED");
    if (seed.empty())
    {
        Vec1Int values;
        now(values);
        // clang-format off
        Int idumlong = values[2] * 24 * 60 * 60 * 300 +
                       values[3] * 60 * 60 * 1000 +
                       values[4] * 60 * 1000 +
                       values[5] * 1000 +
                       values[6];
        // clang-format on
        seed.resize(3, 0);
        seed[0] = idumlong % maxSeed0;
    }
    else if (seed.size() != 3)
    {
        global::tutor.error("The random number requires three integer values.");
    }
    else if (seed[0] < 0 || seed[0] > maxSeed0 || seed[1] < 0 || seed[1] > maxSeed1 || seed[2] < 0)
    {
        global::tutor.error(
            "The Random Number Generator may be seeded as follows:\n 0 <= SEED1 < 900000000\n 0 <= "
            "SEED2 < 1000000\n 0 <= SEED3\nYour current RANDOM_SEED: "
            + std::to_string(seed[0]) + " " + std::to_string(seed[1]) + " "
            + std::to_string(seed[2]) + ".");
    }
    if (istart == 3)
    {
        setTagValue<bool>(setup, incar, "ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS", false);
    }

    setTagValue<Real>(setup, incar, "ML_SRPOT_B0", 1.0);
    setTagValue<Real>(setup, incar, "ML_SRPOT_S0", 2.0);
    setTagValue<Int>(setup, incar, "ML_SRPOT_N0", 10);
    setTagValue<Real>(setup, incar, "ML_EMPPOT_RCUT", 0.0);
    setTagValue<Int>(setup, incar, "ML_MOPOT_NM", 0);
    const Int& mopot_nm = setup.cget<Int>("ML_MOPOT_NM");
    setTagValue<Vec1Int>(setup, incar, "ML_MOPOT_IJM", Vec1Int(2 * mopot_nm));
    setTagValue<Vec1Real>(setup, incar, "ML_MOPOT_DM", Vec1Real(mopot_nm));
    setTagValue<Vec1Real>(setup, incar, "ML_MOPOT_RM", Vec1Real(mopot_nm));
    setTagValue<Vec1Real>(setup, incar, "ML_MOPOT_RKM", Vec1Real(mopot_nm));

    // TODO: Discuss prelimiary tag names!
    setTagValue<String>(setup, incar, "ML_RALGO", "design");
    setTagValue<String>(setup, incar, "ML_SALGO", "ranse");
    setTagValue<String>(setup, incar, "ML_SMETRIC", "pkernel4");
    setTagValue<Vec1Int>(setup, incar, "ML_SNCONF", {});
    setTagValue<Vec1Real>(setup, incar, "ML_STRESH", {});
    setTagValue<String>(setup, incar, "ML_SPAR", "mpi");

    /*============================================================================================+
     | POST-READER SANITY CHECKS AND WARNINGS
     | TODO: Move to separate function, see settings::checkFeatureSupport().
     +============================================================================================*/
    if (type == "grace")
    {
        if (istart != 2)
        {
            global::tutor.error(
                "ML_TYPE = grace is currently only supported in combination with ML_MODE = run "
                "(prediction-only mode). Please pick a valid selection in your INCAR file.");
        }

        String& palgo = setup.get<String>("ML_PALGO");
        if (palgo != "off")
        {
            global::tutor.warning(
                "ML_TYPE = grace is only compatible with ML_PALGO = off, INCAR value \""
                + palgo + "\" was automatically changed to \"off\".");
            palgo = "off";
            states.get<Int>("ML_PALGO") = static_cast<Int>(TagState::OVERRIDE);
        }
    }

    return;
}

void io::writeMllogfile(const Record& setup, std::ostream& output)
{
    output << detail::mllogfileSection('*', "MACHINE LEARNING SETTINGS");
    output << "This section lists the available machine-learning related settings with a short "
              "description, their\nselected values and the INCAR tags. The column between the "
              "value and the INCAR tag may contain a\n\"state indicator\" highlighting the origin "
              "of the value. Here is a list of possible indicators:\n";
    output << "\n";
    output << " * " << TagState::DEFAULT
           << " : (empty) Tag was not provided in the INCAR file, a default value was chosen "
              "automatically.\n";
    output << " * " << TagState::INCAR << " : Value was provided in the INCAR file.\n";
    output << " * " << TagState::INCAR_ALT
           << " : Value was provided in the INCAR file, deprecated tag.\n";
    output << " * " << TagState::OVERRIDE
           << " : A value found in the INCAR file was overwritten by the contents of the ML_FF or "
              "ML_AB file.\n";
    output
        << " * " << TagState::UNSET
        << " : The value for this tag was never set (please report this to the VASP developers).\n";
    output << "\n";
    output << "Tag values with associated units are given here in Angstrom/eV, if not specified "
              "otherwise.\n";
    output << "\n";
    output << "Please refer to the VASP online manual for a detailed description of available "
              "INCAR tags.\n";
    output << "\n";
    output << "\n";
    output << detail::mllogfileSection('-', "General settings");
    output << detail::mllogfileTag("ML_MODE", setup);
    output << detail::mllogfileTag("ML_ISTART", setup);
    output << detail::mllogfileTag("ML_TYPE", setup);
    output << detail::mllogfileTag("ML_OUTPUT_MODE", setup);
    output << detail::mllogfileTag("ML_OUTBLOCK", setup);

    const String& type = setup.cget<String>("ML_TYPE");
    if (type == "kernel") writeMllogfileKernel(setup, output);
    else if (type == "grace") writeMllogfileGrace(setup, output);

    output << "\n";
    output << setup.dcget<ShRec>("incar").cget<String>("message");
    output << "\n";
    output << detail::mllogfileSection('*');

    return;
}
