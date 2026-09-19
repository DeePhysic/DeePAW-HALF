#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Predictor

#include "test_Predictor.hpp"

#include "boost_helpers.hpp"

#include "BasisFunctions.hpp"
#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "ParallelEnvironment.hpp"
#include "Predictor.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include "boost_unit_test_case.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCase_Predictor> container;

void prediction(String                    palgo,
                const TestCase_Predictor& testCase,
                Real&                     energy,
                Vec1Real&                 stress,
                Vec1Real&                 forces)
{
    global::parallel.init(palgo);

    const Record& ff = testCase.ff;

    BasisFunctionMap basisFunctions = setup::makeBasisFunctions(ff);

    Frame frame;
    frame.init(ff, basisFunctions);
    frame.update(testCase.structure, 1.0 / constants::AUTOA);

    KernelPolynomial kernel(ff);
    kernel.updatePolynomialKernel(frame);

    Predictor predictor(ff);
    predictor.update(kernel);
    Real volume = kernel.get_descriptorCollection()
                      .getDescriptor("SHS2-2-body")
                      .get_neighborList_ptr()
                      ->get_latticeVolume();
    predictor.compute_atomicForces(kernel.get_descriptorCollection());
    predictor.compute_stressTensor(kernel.get_descriptorCollection(), volume);

    energy = predictor.get_totalEnergy() * constants::EUNIT;
    std::cout << str("        Absolute difference ENERGY = %16.8E (limit: %16.8E)\n",
                     std::abs(energy - testCase.energy),
                     testCase.toleranceEnergy);

    stress = predictor.get_totalStressTensor();
    for (Real& s : stress) s *= constants::SUNIT;
    Real maxAbsDiffStress = 0.0;
    for (Size i = 0; i < stress.size(); ++i)
    {
        maxAbsDiffStress = std::max(maxAbsDiffStress, std::abs(stress[i] - testCase.stress[i]));
    }
    std::cout << str("Maximum absolute difference STRESS = %16.8E (limit: %16.8E)\n",
                     maxAbsDiffStress,
                     testCase.toleranceStress);

    forces = predictor.get_atomicForces();
    for (Real& f : forces) f *= constants::FUNIT;
    Real maxAbsDiffForces = 0.0;
    for (Size i = 0; i < forces.size(); ++i)
    {
        maxAbsDiffForces = std::max(maxAbsDiffForces, std::abs(forces[i] - testCase.forces[i]));
    }
    std::cout << str("Maximum absolute difference FORCES = %16.8E (limit: %16.8E)\n",
                     maxAbsDiffForces,
                     testCase.toleranceForces);
    return;
}

BOOST_AUTO_TEST_SUITE(IntegrationTests)

BOOST_DATA_TEST_CASE(PredictPalgoOff_CorrectEnergiesForcesStress,
                     bdata::make(container.testCases),
                     testCase)
{
    Real energy;
    Vec1Real stress;
    Vec1Real forces;
    prediction("off", testCase, energy, stress, forces);

    BOOST_REQUIRE_SMALL(energy - testCase.energy, testCase.toleranceEnergy);
    REQUIRE_CLOSE_COLLECTIONS(stress, testCase.stress, testCase.toleranceStress, "stress");
    REQUIRE_CLOSE_COLLECTIONS(forces, testCase.forces, testCase.toleranceForces, "forces");
}

BOOST_DATA_TEST_CASE(PredictPalgoSerial_CorrectEnergiesForcesStress,
                     bdata::make(container.testCases),
                     testCase)
{
    if (testCase.name.find("DESC_TYPE_1") == String::npos)
    {
        Real energy;
        Vec1Real stress;
        Vec1Real forces;
        prediction("serial", testCase, energy, stress, forces);

        BOOST_REQUIRE_SMALL(energy - testCase.energy, testCase.toleranceEnergy);
        REQUIRE_CLOSE_COLLECTIONS(stress, testCase.stress, testCase.toleranceStress, "stress");
        REQUIRE_CLOSE_COLLECTIONS(forces, testCase.forces, testCase.toleranceForces, "forces");
    }
    else
    {
        BOOST_TEST_MESSAGE("Force fields with DESC_TYPE = 1 are not yet supported with "
                           "PALGO = serial.");
    }
}

BOOST_AUTO_TEST_SUITE_END()
