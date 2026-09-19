#include "Inference.hpp"

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "Lattice.hpp"
#include "Record.hpp"
#include "Structure.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iostream>

using namespace vaspml;
using namespace vaspml::setup;

Inference::Inference(ShRec data) :
    data((data == nullptr) ? std::make_shared<Record>() : data),
    frame(data),
    polyKernel(data),
    predictor(data)
{
    this->data->add("predictions", "Vec1ShRec");
}

void Inference::init(const Record& mlff, const std::shared_ptr<MlMPI> mpiIn)
{
    basisFunctions = setup::makeBasisFunctions(mlff);
    frame.init(mlff, basisFunctions);
    polyKernel.init(mlff, mpiIn);
    predictor.init(mlff, mpiIn);
    mpi = mpiIn;
}

void Inference::predict(const String& poscarName, const Real& distScale)
{
    frame.update(poscarName, distScale);
    predict();
    addPredictedData(poscarName);
}

void Inference::predict(const Vec1String& poscars, const Real& distScale)
{
    std::for_each(poscars.cbegin(),
                  poscars.end(),
                  [&](const String& poscarName) { predict(poscarName, distScale); });
}

void Inference::predict(std::vector<Structure>& structure)
{
    std::for_each(structure.begin(), structure.end(), [&](Structure& struc) { predict(struc); });
}

void Inference::predict(Structure& structure)
{
    frame.update(structure);

    //const Vec1Real pos = structure.get_positions();
    //Int i = 0;
    //for ( const auto& x : pos )
    //{
    //    std::cout << x << "   ";
    //    i++;
    //    if ( i%3==0 ) std::cout << std::endl;
    //}
    predict();
    addPredictedData("sample system");
}

void Inference::predict()
{
    polyKernel.updatePolynomialKernel(frame);

    Real                  volume = 0.0;
    [[maybe_unused]] bool descFound = false;
    for (const String& key : constants::descriptorKeyList)
    {
        if (polyKernel.get_descriptorCollection().getDescriptor(key).get_weight() > 0)
        {
            descFound = true;
            volume = polyKernel.get_descriptorCollection()
                         .getDescriptor(key)
                         .get_neighborList_ptr()
                         ->get_latticeVolume();
            break;
        }
    }
    VASPML_DEBUG_L1(
        if (!descFound)
        {
            global::tutor.error("ERROR: " + flf(VASPML_FLF) + " no valid descriptor found!");
        }
    );
    predictor.update(polyKernel);
    predictor.compute_atomicForces(polyKernel.get_descriptorCollection());
    predictor.compute_stressTensor(polyKernel.get_descriptorCollection(), volume);
}

void Inference::writeCurrentResultToScreen(const Real& energyScale,
                                           const Real& forceScale,
                                           const Real& stressScale) const
{
    Size numAtoms = polyKernel.get_descriptorCollection()
                        .getDescriptor("SHS2-2-body")
                        .get_neighborList_ptr()
                        ->get_nAtoms();
    std::cout << str("NATOM %zu\n", numAtoms);
    std::cout << str("ENERGY %24.16E\n", predictor.get_totalEnergy() * energyScale);
    Real pressure = 0.0;
    for (Size i = 0; i < 3; i++)
    {
        std::cout << "STRESS ";
        for (Size j = 0; j < 3; j++)
        {
            std::cout << str(" %24.16E", predictor.get_totalStressTensor(i, j) * stressScale);
        }
        pressure += predictor.get_totalStressTensor(i, i);
        std::cout << "\n";
    }
    pressure *= stressScale / 3.0;
    std::cout << str("PRESSURE %24.16E\n", pressure);

    const Vec1Real fxyz = predictor.get_atomicForces();
    for (Size atom = 0; atom < numAtoms; atom++)
    {
        std::cout << str("FORCE %-6zu %24.16E %24.16E %24.16E\n",
                         atom + 1,
                         fxyz[3 * atom] * forceScale,
                         fxyz[3 * atom + 1] * forceScale,
                         fxyz[3 * atom + 2] * forceScale);
    }
}

void Inference::addPredictedData(const String& key)
{
    Vec1ShRec& predictions = data->get<Vec1ShRec>("predictions");
    predictions.push_back(std::make_shared<Record>());
    Record& structure = *predictions.back();
    structure.put<String>("system", key);
    const Vec1Int& numAtomsPerType = frame.get_structure()->get_numAtomsPerType();
    structure.put<Int>("numTypes", numAtomsPerType.size());
    structure.put<Vec1Int>("numAtomsPerType", numAtomsPerType);
    structure.put<Int>("numAtoms", frame.get_structure()->get_numAtoms());
    structure.put<Vec1String>("typeNames", frame.get_structure()->get_typeNames());
    structure.put<Vec1Real>("lattice", frame.get_structure()->get_lattice().get_components());
    structure.put<Vec1Real>("positions", frame.get_structure()->get_positions());
    structure.put<Real>("energy", predictor.get_totalEnergy());
    structure.put<Vec1Real>("forces", predictor.get_atomicForces());
    const Vec1Real& stress = predictor.get_totalStressTensor();
    // xx xy xz, yx, yz, zx, zy, zz
    // 0   1  2, 3, 4, 5, 6, 7, 8
    structure.put<Vec1Real>("stress",
                            {stress[0], stress[1], stress[2], stress[4], stress[5], stress[8]});
}

void Inference::writePredictionsToJson([[maybe_unused]] const String& fname) const
{
    //io::open(f, fname, "w");
    //rec::toJson(record, f);
    //io::close(f);
    //TODO
    //we would need a ShVec1Rec to json then this would be sound
}

ShRec Inference::getPredictions() const
{
    return data;
}
