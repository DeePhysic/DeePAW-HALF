#include "setup.hpp"

#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "DescriptorSHS3ReducedLinElem.hpp"
#include "Record.hpp"
#include "Selector.hpp"
#include "ShmemArray.hpp"
#include "Structure.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "math.hpp"

#include <algorithm> // for for_each

using namespace vaspml;

BasisFunctionMap setup::makeBasisFunctions(const Record& record)
{
    using CT = vaspml::math::CutoffType;

    BasisFunctionMap basisFunctions;
    basisFunctions["2-body"] =
        std::make_shared<BasisFunctionsRadialSpline>(static_cast<CT>(record.cget<Int>("ML_ICUT1")),
                                                     record.cget<Real>("ML_RCUT1"),
                                                     record.cget<Real>("ML_SION1"),
                                                     record.cget<Int>("ML_NR1"),
                                                     0,
                                                     record.cget<Int>("ML_MRB1"),
                                                     BasisFunctionType::bodyOrder2);
    basisFunctions["3-body"] =
        std::make_shared<BasisFunctionsAngular>(static_cast<CT>(record.cget<Int>("ML_ICUT2")),
                                                record.cget<Real>("ML_RCUT2"),
                                                record.cget<Real>("ML_SION2"),
                                                record.cget<Int>("ML_NR2"),
                                                record.cget<Int>("ML_LMAX2"),
                                                record.cget<Int>("ML_MRB2"),
                                                BasisFunctionType::bodyOrder3);

    return basisFunctions;
}

std::vector<Structure> setup::makeStructureStorage(const ShRec& mlab)
{
    Int                    nStrucs = mlab->cget<Int>("numStructures");
    std::vector<Structure> structures;
    Vec1ShRec&             inputStructures = mlab->get<Vec1ShRec>("structures");
    for (Int i = 0; i < nStrucs; i++)
    {
        structures.push_back(Structure(inputStructures[i], false));
    }

    return structures;
}

NeighborListMap setup::createNeighborLists(const std::map<String, Real>& cutoff, ShRec data)
{
    NeighborListMap neighborLists;
    neighborLists["2-body"] =
        std::make_shared<NearestNeighborNSquare>(cutoff.at("2-body"),
                                                 true,
                                                 false,
                                                 data->get<ShRec>("2-body-neighbor-list"));

    neighborLists["3-body"] =
        std::make_shared<NearestNeighborNSquare>(cutoff.at("3-body"),
                                                 true,
                                                 false,
                                                 data->get<ShRec>("3-body-neighbor-list"));

    return neighborLists;
}

DescriptorMap setup::createDescriptors(const Record&              incar,
                                       const Size&                nTypes,
                                       const BasisFunctions&      basis3Body,
                                       const Descriptor3BodyType& descriptor3BodyType,
                                       ShRec                      data)
{
    DescriptorMap descriptors;

    descriptors["SHS2-2-body"] = std::make_shared<DescriptorSHS2>(incar.cget<Real>("ML_W1"),
                                                                  incar.cget<bool>("ML_LNORM1"),
                                                                  DescriptorType::bodyOrder2,
                                                                  data->get<ShRec>("2-body-SHS2"));

    descriptors["SHS2-3-body"] = std::make_shared<DescriptorSHS2>(incar.cget<Real>("ML_W2"),
                                                                  incar.cget<bool>("ML_LNORM2"),
                                                                  DescriptorType::bodyOrder2,
                                                                  data->get<ShRec>("3-body-SHS2"));
    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        descriptors["SHS3-3-body"] =
            std::make_shared<DescriptorSHS3>(nTypes,
                                             incar.cget<bool>("ML_LAFILT2"),
                                             incar.cget<Int>("ML_IAFILT2"),
                                             incar.cget<Real>("ML_AFILT2"),
                                             basis3Body.get_nValidRoots(),
                                             incar.cget<Int>("ML_LMAX2"),
                                             incar.cget<Real>("ML_W2"),
                                             incar.cget<bool>("ML_LNORM2"),
                                             DescriptorType::bodyOrder3,
                                             data->get<ShRec>("3-body-SHS3"));
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        descriptors["SHS3-3-body"] =
            std::make_shared<DescriptorSHS3ReducedLinElem>(nTypes,
                                                           incar.cget<bool>("ML_LAFILT2"),
                                                           incar.cget<Int>("ML_IAFILT2"),
                                                           incar.cget<Real>("ML_AFILT2"),
                                                           basis3Body.get_nValidRoots(),
                                                           incar.cget<Int>("ML_LMAX2"),
                                                           incar.cget<Real>("ML_W2"),
                                                           incar.cget<bool>("ML_LNORM2"),
                                                           DescriptorType::bodyOrder3,
                                                           data->get<ShRec>("3-body-SHS3"));
    }

    return descriptors;
}

DescriptorMap setup::createDescriptors(const Record&              setup,
                                       const BasisFunctions&      basis3Body,
                                       const Descriptor3BodyType& descriptor3BodyType,
                                       ShRec                      data)
{
    DescriptorMap descriptors;

    descriptors["SHS2-2-body"] = std::make_shared<DescriptorSHS2>(setup.cget<Real>("ML_W1"),
                                                                  setup.cget<bool>("ML_LNORM1"),
                                                                  DescriptorType::bodyOrder2,
                                                                  data->get<ShRec>("2-body-SHS2"));

    descriptors["SHS2-3-body"] = std::make_shared<DescriptorSHS2>(setup.cget<Real>("ML_W2"),
                                                                  setup.cget<bool>("ML_LNORM2"),
                                                                  DescriptorType::bodyOrder2,
                                                                  data->get<ShRec>("3-body-SHS2"));
    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        if (setup.cget<Real>("ML_W2") > 0)
        {
            descriptors["SHS3-3-body"] =
                std::make_shared<DescriptorSHS3>(setup.cget<Vec2Int>("SHS3-3-body-descriptor-list"),
                                                 setup.cget<bool>("ML_LAFILT2"),
                                                 setup.cget<Int>("ML_IAFILT2"),
                                                 setup.cget<Real>("ML_AFILT2"),
                                                 basis3Body.get_nValidRoots(),
                                                 setup.cget<Int>("ML_LMAX2"),
                                                 setup.cget<Real>("ML_W2"),
                                                 setup.cget<bool>("ML_LNORM2"),
                                                 DescriptorType::bodyOrder3,
                                                 data->get<ShRec>("3-body-SHS3"));
        }
        else descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3>();
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        if (setup.cget<Real>("ML_W2") > 0)
        {
            descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3ReducedLinElem>(
                setup.cget<Vec2Int>("SHS3-3-body-descriptor-list"),
                setup.cget<bool>("ML_LAFILT2"),
                setup.cget<Int>("ML_IAFILT2"),
                setup.cget<Real>("ML_AFILT2"),
                basis3Body.get_nValidRoots(),
                setup.cget<Int>("ML_LMAX2"),
                setup.cget<Real>("ML_W2"),
                setup.cget<bool>("ML_LNORM2"),
                DescriptorType::bodyOrder3LinearElement,
                data->get<ShRec>("3-body-SHS3LinElement"));
        }
        else descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3ReducedLinElem>();
    }

    return descriptors;
}

std::map<String, ShVec1Int> setup::make_numberDescriptorsMap(const Record& setup)
{
    std::map<String, ShVec1Int> numberDescriptors;

    if (setup.cget<Real>("ML_W2") > 0)
        numberDescriptors["SHS3-3-body"] = std::make_shared<Vec1Int>(
            setup.cget<Vec1Int>("SHS3-3-body-number-descriptors-per-type"));
    else numberDescriptors["SHS3-3-body"] = std::make_shared<Vec1Int>();

    if (setup.cget<Real>("ML_W1") > 0)
    {
        Int     twoBody = setup.cget<Int>("SHS2-2-body-number-descriptors-per-type");
        Vec1Int twoBodyVector(setup.cget<Int>("numTypes"));
        for (Size i = 0; i < twoBodyVector.size(); i++) { twoBodyVector[i] = twoBody; }
        numberDescriptors["SHS2-2-body"] = std::make_shared<Vec1Int>(twoBodyVector);
    }
    else numberDescriptors["SHS2-2-body"] = std::make_shared<Vec1Int>();

    return numberDescriptors;
}

std::map<String, ShVec1Int> setup::generate_featureSpaceMap(const Record& setup)
{
    std::map<String, ShVec1Int> featureSpaceSize;
    const Vec1Int&              locRef = setup.cget<Vec1Int>("numLrc");

    // TODO make this to a getif with the record iff possible
    //for (const String& key : constants::descriptorKeyList)
    //{
    //    String tag = key + "-number-descriptors-per-type";
    //    if (std::get_if<Int>(&data.at(translator.get_InputTag(tag))))
    //    {
    //        featureSpaceSize[key] = std::make_shared<Vec1Int>(
    //            math::vectorTimesScalar(locRef,
    //            data.at(translator.get_InputTag(tag)).cget<Int>()));
    //    }
    //    else if (std::get_if<ShVec1Int>(&data.at(translator.get_InputTag(tag))))
    //    {
    //        featureSpaceSize[key] = std::make_shared<Vec1Int>(
    //            math::elementwiseProduct(locRef,
    //                                     data.at(translator.get_InputTag(tag)).dcget<ShVec1Int>()));
    //    }
    //}

    if (setup.cget<Real>("ML_W1") > 0)
    {
        featureSpaceSize["SHS2-2-body"] = std::make_shared<Vec1Int>(
            math::vectorTimesScalar(locRef,
                                    setup.cget<Int>("SHS2-2-body-number-descriptors-per-type")));
    }
    else featureSpaceSize["SHS2-2-body"] = std::make_shared<Vec1Int>();
    if (setup.cget<Real>("ML_W2") > 0)
        featureSpaceSize["SHS3-3-body"] = std::make_shared<Vec1Int>(math::elementwiseProduct(
            locRef,
            setup.cget<Vec1Int>("SHS3-3-body-number-descriptors-per-type")));
    else featureSpaceSize["SHS3-3-body"] = std::make_shared<Vec1Int>();

    return featureSpaceSize;
}

std::shared_ptr<ShmemArray2DVariableLen<Real>> setup::make_regressionCoefficients(
    const Record&                 inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn)
{
    const Vec2Real& regCoeff = inputParameters.cget<Vec2Real>("regression-coeff");
    const Vec1Int&  sizes = inputParameters.cget<Vec1Int>("numLrc");
    std::shared_ptr<ShmemArray2DVariableLen<Real>> shmemPtr =
        std::make_shared<ShmemArray2DVariableLen<Real>>(sizes, mpiIn);
    shmemPtr->set_value(regCoeff);

    return shmemPtr;
}

std::map<String, Real> setup::make_weightsMap(const Record& inputParameters)
{
    std::map<String, Real> weights;
    weights["SHS2-2-body"] = inputParameters.cget<Real>("ML_W1");
    weights["SHS3-3-body"] = inputParameters.cget<Real>("ML_W2");

    return weights;
}

Selector setup::makeSelector(const ShRec&                  incar,
                             const ShRec&                  mlab,
                             const ShRec&                  data,
                             const std::shared_ptr<MlMPI>& mlmpi)
{
    Selector selector(incar, mlab, data, mlmpi);
    //selector.init( incar, mlab, data, mlmpi );

    return selector;
}

Vec1Int setup::makeFarthestPointNumberSamples(const ShRec& incar, const ShRec& mlab)
{
    Vec1Int   numSamples = incar->get<Vec1Int>("ML_SNCONF");
    const Int numberTypes = mlab->cget<Int>("maxTypes");
    if (numSamples.size() != (Size)numberTypes)
    {
        global::tutor.error(
            "ERROR: " + flf(VASPML_FLF)
            + " the number of local reference configuration selections does"
              " not match the number of atom types in the ML_AB file. Please adjust the number of"
            + " values supplied to ML_SNCONF.");
    }

    return numSamples;
}

Vec1Real setup::makeFarthestPointTresholds(const ShRec& incar, const ShRec& mlab)
{
    Vec1Real  tresh = incar->get<Vec1Real>("ML_STRESH");
    const Int numberTypes = mlab->cget<Int>("maxTypes");
    if (tresh.size() != (Size)numberTypes)
    {
        global::tutor.error(
            "ERROR: " + flf(VASPML_FLF)
            + " the number of local reference configuration selections does"
              " not match the number of atom types in the ML_AB file. Please adjust the number of"
            + " values supplied to ML_STRESH.");
    }

    return tresh;
}

std::variant<Vec1Int, Vec1Real> setup::makeFarhtestPointCriterion(const ShRec& incar,
                                                                  const ShRec& mlab)
{
    if (incar->cget<String>("ML_SALGO") == "fpsn" || incar->cget<String>("ML_SALGO") == "rann"
        || incar->cget<String>("ML_SALGO") == "slhc")
    {
        return makeFarthestPointNumberSamples(incar, mlab);
    }
    else return makeFarthestPointTresholds(incar, mlab);
}

void setup::rescaleIncarUnits(Record& incar, const Real distFactor)
{
    incar.get<Real>("ML_RCUT1") *= distFactor;
    incar.get<Real>("ML_RCUT2") *= distFactor;
    incar.get<Real>("ML_SION1") *= distFactor;
    incar.get<Real>("ML_SION2") *= distFactor;

    return;
}

void setup::rescaleMlabUnits(Record&    mlab,
                             const Real energyFactor,
                             const Real distFactor,
                             const Real forceFactor,
                             const Real stressFactor)
{
    math::vectorTimesScalarNoCopy(mlab.get<Vec1Real>("atomicRefEnergy"), energyFactor);
    Vec1ShRec& structures = mlab.get<Vec1ShRec>("structures");
    std::for_each(structures.begin(),
                  structures.end(),
                  [&](ShRec& struc)
                  {
                      math::vectorTimesScalarNoCopy(struc->get<Vec1Real>("positions"), distFactor);
                      struc->get<Real>("energy") *= energyFactor;
                      math::vectorTimesScalarNoCopy(struc->get<Vec1Real>("forces"), forceFactor);
                      math::vectorTimesScalarNoCopy(struc->get<Vec1Real>("stress"), stressFactor);
                      math::vectorTimesScalarNoCopy(struc->get<Vec1Real>("lattice"), distFactor);
                  });

    return;
}

void setup::rescaleMlffUnits(Record&    mlff,
                             const Real energyFactor,
                             const Real distFactor,
                             const Real forceFactor,
                             const Real stressFactor)
{
    math::vectorTimesScalarNoCopy(mlff.get<Vec1Real>("ML_EATOM_REF"), energyFactor);
    mlff.get<Real>("ML_WTOTEN") *= energyFactor;
    mlff.get<Real>("ML_WTIFOR") *= forceFactor;
    mlff.get<Real>("ML_WTSIF") *= stressFactor;
    mlff.get<Real>("rmseEnergy") *= energyFactor;
    mlff.get<Real>("rmseForces") *= forceFactor;
    mlff.get<Real>("rmseStress") *= stressFactor;
    mlff.get<Real>("ML_RCUT1") *= distFactor;
    mlff.get<Real>("ML_RCUT2") *= distFactor;
    mlff.get<Real>("ML_SION1") *= distFactor;
    mlff.get<Real>("ML_SION2") *= distFactor;
    mlff.get<Real>("average-energy-per-atom") *= energyFactor;
    math::vectorTimesScalarNoCopy(mlff.get<Vec2Real>("regression-coeff"), energyFactor);
    mlff.get<Real>("variance-training-energies") *= energyFactor;
    math::vectorTimesScalarNoCopy(mlff.get<Vec1Real>("variance-training-force"), forceFactor);
    math::vectorTimesScalarNoCopy(mlff.get<Vec1Real>("variance-training-stress"), stressFactor);

    return;
}
