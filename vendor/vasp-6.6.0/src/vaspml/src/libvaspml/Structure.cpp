#include "Structure.hpp"

#include "Record.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <cmath>
#include <iostream>

using namespace vaspml;

MapString vaspml::dataMapStructure()
{
    return MapString{
        {"positions",           "Vec1Real"  },
        {"typeNames",           "Vec1String"},
        {"atomTypeNames",       "Vec1String"},
        {"typeIDs",             "Vec1Int"   },
        {"numAtomsPerType",     "Vec1Int"   },
        {"atomIndx",            "Vec1Int"   },
        {"latticeShRec",        "ShRec"     },
        {"inverseLatticeShRec", "ShRec"     }
    };
}

Structure::Structure(ShRec record, const bool& isDirect) :
    data(assignOrMakeRecord(record, dataMapStructure())),
    direct(isDirect),
    lattice(data->get<ShRec>("latticeShRec")),
    inverseLattice(data->get<ShRec>("inverseLatticeShRec")),
    positions(data->get<Vec1Real>("positions")),
    typeNames(data->get<Vec1String>("typeNames")),
    typeIDs(data->get<Vec1Int>("typeIDs")),
    atomTypeNames(data->get<Vec1String>("atomTypeNames")),
    numAtomsPerType(data->get<Vec1Int>("numAtomsPerType")),
    atomIndx(data->get<Vec1Int>("atomIndx"))
{
    if (data->contains("lattice")) set_lattice(data->get<Vec1Real>("lattice"));
    set_types(typeNames, numAtomsPerType);
}

void Structure::readPoscar(const String& fname)
{
    std::ifstream infile_stream = file_io::openFileI(fname);

    String dummy;
    // read header from POSCAR file
    getline(infile_stream, dummy);
    read_lattice(infile_stream);
    read_types(infile_stream);
    getline(infile_stream, dummy);

    direct = string_tools::checkForSubstring(string_tools::makeLowerCase(dummy), "d");
    read_positions(infile_stream);
    if (direct) { shift_atoms_primitive_cell(); }
    else
    {
        cartesianToDirect();
        shift_atoms_primitive_cell();
        directToCartesian();
    }

    build_types();
    infile_stream.close();

    return;
}

void Structure::set_lattice(const Vec1Real& lattice_in)
{
    lattice.set_components(lattice_in[0],
                           lattice_in[1],
                           lattice_in[2],
                           lattice_in[3],
                           lattice_in[4],
                           lattice_in[5],
                           lattice_in[6],
                           lattice_in[7],
                           lattice_in[8]);
    inverseLattice.set_components(lattice.get_components());
    inverseLattice.invert();

    return;
}

void Structure::set_positions(const Vec1Real& positions, const bool isDirect_in)
{
    direct = isDirect_in;

    this->positions = positions;
    if (direct) { shift_atoms_primitive_cell(); }
    else
    {
        cartesianToDirect();
        shift_atoms_primitive_cell();
        directToCartesian();
    }

    return;
}

void Structure::set_types(const Vec1String& atom_types_in, const Vec1Int& number_atoms_in)
{
    typeNames = atom_types_in;
    numAtomsPerType = number_atoms_in;
    numAtoms = math::sumVector(number_atoms_in);
    build_types();

    return;
}

void Structure::importPoscarData(const Vec1Real&   lattice_in,
                                 const Vec1String& atom_types_in,
                                 const Vec1Int&    number_atoms_in,
                                 const Vec1Real&   positions,
                                 const bool        isDirect_in)
{
    set_lattice(lattice_in);
    set_types(atom_types_in, number_atoms_in);
    set_positions(positions, isDirect_in);
    build_types();

    return;
}

void Structure::directToCartesian()
{
    if (direct) { rescale_coordinates(lattice); }
    else
    {
        std::cout << "Warning in Structure::position_direct_to_cartesian" << std::endl;
        std::cout << "You already have direct coordinates. Doing nothing" << std::endl;
    }
    direct = false;

    return;
}

void Structure::cartesianToDirect()
{
    if (!direct) { rescale_coordinates(inverseLattice); }
    else
    {
        std::cout << "Warning in Structure::position_cartesian_to_direct" << std::endl;
        std::cout << "You already have direct coordinates. Doing nothing" << std::endl;
    }
    direct = true;

    return;
}

const std::tuple<const Real&, const Real&, const Real&> Structure::getAtom(
    const Size atomIndx) const
{
    return std::tie(positions[3 * atomIndx],
                    positions[3 * atomIndx + 1],
                    positions[3 * atomIndx + 2]);
}

bool Structure::isDirect() const
{
    return direct;
}

const Vec1Real& Structure::get_positions() const
{
    return positions;
}

const Lattice& Structure::get_lattice() const
{
    return lattice;
}

const Lattice& Structure::get_inverseLattice() const
{
    return inverseLattice;
}

const Vec1Int& Structure::get_types() const
{
    return typeIDs;
}

const Vec1String& Structure::get_typeNames() const
{
    return typeNames;
}

const Vec1String& Structure::get_atomTypeNames() const
{
    return atomTypeNames;
}

const String& Structure::get_atomTypeNames(const Size& atomIdx) const
{
    return atomTypeNames[atomIdx];
}

const Vec1Int& Structure::get_numAtomsPerType() const
{
    return numAtomsPerType;
}

Size Structure::get_numAtoms() const
{
    return numAtoms;
}

void Structure::writeToScreen() const
{
    std::cout << "-------------------------------------------\n";
    std::cout << "WARNING: All output uses Bohr length units.\n";
    std::cout << "-------------------------------------------\n";
    std::cout << str("Scaling factor: %24.16E\n", scale);
    std::cout << "Lattice:\n";
    lattice.writeToScreen();
    std::cout << "Inverse lattice:\n";
    inverseLattice.writeToScreen();
    std::cout << "(lattice) * (inverse lattice):\n";
    for (Size i = 0; i < 3; ++i)
    {
        for (Size j = 0; j < 3; ++j)
        {
            Real mij = 0.0;
            for (Size k = 0; k < 3; ++k)
            {
                mij += lattice.component(i, k) * inverseLattice.component(k, j);
            }
            std::cout << str(" %10.2E", mij);
        }
        std::cout << "\n";
    }
    std::cout << str("Number of atoms: %zu\n", numAtoms);
    std::cout << "Atom types:";
    for (auto t : typeNames) std::cout << " " + t;
    std::cout << "\n";
    std::cout << "Number of atoms per type:";
    for (auto n : numAtomsPerType) std::cout << str(" %d", n);
    std::cout << "\n";
    if (direct) std::cout << "Direct\n";
    else std::cout << "Cartesian\n";

    std::cout << "positions " << positions.size() << std::endl;
    std::cout << "typeIDs " << typeIDs.size() << std::endl;
    std::cout << "typeNames " << typeNames.size() << std::endl;

    for (Size i = 0; i < numAtoms; ++i)
    {
        std::cout << str("%24.16E %24.16E %24.16E %d %s\n",
                         positions[i],
                         positions[i + 1],
                         positions[i + 2],
                         typeIDs[i],
                         typeNames[typeIDs[i]].c_str());
    }

    return;
}

//**************************************************************************************************
// Private member functions
//**************************************************************************************************

void Structure::read_lattice(std::ifstream& infile_stream)
{
    String line;
    getline(infile_stream, line);
    line = string_tools::trim(line);
    Vec1Real scaleVector = string_tools::extractData<Real>(line);
    if (scaleVector.size() != 1)
    {
        global::tutor.error("Unsupported scaling factor format encountered while reading "
                            "POSCAR, only a single value is supported.");
    }
    scale = scaleVector.at(0);
    Vec1Real temp_lattice;

    for (auto i = 0; i < 3; i++)
    {
        getline(infile_stream, line);
        Vec1Real tmp = string_tools::extractData<Real>(line);
        math::vectorTimesScalarNoCopy(tmp, scale);
        temp_lattice.insert(temp_lattice.end(), tmp.begin(), tmp.end());
    }

    lattice.set_components(temp_lattice);
    inverseLattice.set_components(lattice.get_components());
    inverseLattice.invert();

    return;
}

void Structure::read_positions(std::ifstream& infile_stream)
{
    numAtoms = math::sumVector(numAtomsPerType);
    positions.resize(3 * numAtoms);
    Size counter = 0;
    for (Size i = 0; i < numAtoms; i++)
    {
        String line;
        getline(infile_stream, line);
        Vec1Real atom = string_tools::extractData<Real>(line);
        if (!direct) math::vectorTimesScalarNoCopy(atom, scale);
        positions[counter] = atom[0];
        counter++;
        positions[counter] = atom[1];
        counter++;
        positions[counter] = atom[2];
        counter++;
    }

    return;
}

void Structure::read_types(std::ifstream& infile_stream)
{
    String line;
    getline(infile_stream, line);
    typeNames = string_tools::extractData<String>(line);
    getline(infile_stream, line);
    numAtomsPerType = string_tools::extractData<Int>(line);

    return;
}

void Structure::shift_atoms_primitive_cell()
{
    Real one = (Real)1;
    if (direct)
    {
        for (Size i = 0; i < positions.size(); i++)
        {
            positions[i] = std::modf(positions[i] + (Real)100, &one);
        }
    }
    else
    {
        std::cout << "Error in Structure::shift_atoms_primitive_cell " << std::endl;
        std::cout << "Function only works on direct coordinates" << std::endl;
    }

    return;
}

void Structure::rescale_coordinates(const Lattice& lattice)
{
    for (Size i = 0; i < 3 * numAtoms; i += 3)
    {
        lattice.timesVectorInPlace(positions[i], positions[i + 1], positions[i + 2]);
    }

    return;
}

void Structure::build_types()
{
    typeIDs.resize(numAtoms);
    atomTypeNames.resize(numAtoms);
    Size indx = 0;
    for (Size i = 0; i < numAtomsPerType.size(); i++)
    {
        for (Size j = 0; j < (Size)numAtomsPerType[i]; j++)
        {
            typeIDs[indx] = i;
            atomTypeNames[indx] = typeNames[i];
            indx++;
        }
    }

    return;
}

void Structure::set_indexOffset(const Size& offset)
{
    indexOffset = offset;
}
Size Structure::get_indexOffset() const
{
    return indexOffset;
}

void Structure::set_typeMap(const Vec1String superTypes)
{
    if (typeMap == nullptr) typeMap = std::make_shared<TypeMap>();
    typeMap->update(superTypes, typeNames);

    return;
}

std::shared_ptr<TypeMap> Structure::get_typeMap() const
{
    return typeMap;
}

void Structure::convertUnits(const Real& factor)
{
    // no rescaling needed if coords are direct. Then units are contained in the lattice.
    if (!direct) math::vectorTimesScalarNoCopy(positions, factor);
    lattice.convertUnits(factor);
    inverseLattice.convertUnits(1 / factor);
}

void Structure::set_atomIndx(const Vec1Int& atomIndx)
{
    this->atomIndx = atomIndx;
}

const Vec1Int& Structure::get_atomIndx() const
{
    return atomIndx;
}

Vec1Int Structure::get_typesLocal() const
{
    return vector_tools::filterVector(typeIDs, atomIndx);
}

//Vec1String Structure::get_atomTypeNamesLocal( void ) const
//{
//    return vector_tools::filterVector( typeNames, atomIndx );
//}
