#include "MlMPI.hpp"

#include "Selector.hpp"

#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "FarthestPointSampling.hpp"
#include "Frame.hpp"
#include "MetricFunctions.hpp"
#include "Record.hpp"
#include "ShmemArray.hpp"
#include "SequentialLikelihoodCoverage.hpp"
#include "Timer.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "dataset.hpp"
#include "debug.hpp"
#include "io.hpp"
#include "math.hpp"
#include "rec_mpi.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm>                         // for for_each
#include <iostream>                          // for basic_ostream, operator<<
#include <limits>                            // for numeric_limits
#include <numeric>                           // for iota
#include <variant>                           // for get, variant

using namespace vaspml;

SelectionMethod vaspml::chooseSelectionMethod(const ShRec& incar)
{
    SelectionMethod method = SelectionMethod::None;
    if (incar->cget<String>("ML_SALGO") == "fpsn")
    {
        method = SelectionMethod::FarthestPointSamplingN;
    }
    else if (incar->cget<String>("ML_SALGO") == "rann")
    {
        method = SelectionMethod::RandomSelection;
    }
    else if (incar->cget<String>("ML_SALGO") == "fpst")
    {
        method = SelectionMethod::FarthestPointSamplingThreshold;
    }
    else if (incar->cget<String>("ML_SALGO") == "slhc")
    {
        method = SelectionMethod::SequentialLikelihoodCover;
    }
    else
    {
        global::tutor.error("ERROR:" + flf(VASPML_FLF) + " selection method not recognized.\n");
    }

    return method;
}

SelectionMetric vaspml::chooseMetric(const ShRec& incar)
{
    SelectionMetric metric = SelectionMetric::None;
    if (incar->cget<String>("ML_SMETRIC") == "l1norm") { metric = SelectionMetric::L1Norm; }
    else if (incar->cget<String>("ML_SMETRIC") == "l2norm") { metric = SelectionMetric::L2Norm; }
    else if (incar->cget<String>("ML_SMETRIC") == "lpnorm") { metric = SelectionMetric::LpNorm; }
    else if (incar->cget<String>("ML_SMETRIC") == "pkernel4")
    {
        metric = SelectionMetric::PolyKernelNorm;
    }
    else if (incar->cget<String>("ML_SMETRIC") == "pcorr")
    {
        metric = SelectionMetric::antiPearsonNorm2;
    }
    else
    {
        global::tutor.error("ERROR:" + flf(VASPML_FLF) + " selection metric not recognized.\n");
    }

    return metric;
}

void vaspml::mpiRankCheck(const Vec1Int& numberAtoms, MlMPI& mlmpi)
{
    Int minAtoms = vector_tools::get_min(numberAtoms);
    if (minAtoms < mlmpi.get_numberRanks())
    {
        global::tutor.error("Error:" + flf(VASPML_FLF)
                            + "You are using more MPI ranks than there are atoms in your ML_AB.\n"
                            + "Please reduce the total number of cores, this is harmful for the"
                            + "environment and won't reduce computation time.\n");
    }

    return;
}

ParallelMode vaspml::get_parallelMode(const Record& incar)
{
    if (incar.cget<String>("ML_SPAR") == "mpi") { return ParallelMode::MPI; }
    else if (incar.cget<String>("ML_SPAR") == "openmp") { return ParallelMode::OpenMP; }
    else
    {
        global::tutor.error("ERROR: " + flf(VASPML_FLF)
                            + "Trying to use Selector with non existent ParallelMode.");
        throw;
    }
}

Selector::Selector(const ShRec& data) :
    incar(std::make_shared<Record>()),
    mlab(std::make_shared<Record>()),
    data(data == nullptr ? std::make_shared<Record>() : data),
    parallelMode(ParallelMode::None),
    randomSeed({0, 0, 0}),
    randomUnique(randomSeed),
    randomUniform(randomSeed),
    doWrite(false),
    descNormalized(false)
{
    selectionMethod = SelectionMethod::None;
    selectionMetric = SelectionMetric::None;
}

Selector::Selector(const ShRec&                  input,
                   const ShRec&                  mlab,
                   const ShRec&                  data,
                   const std::shared_ptr<MlMPI>& mlmpi) :
    incar(input),
    mlab(mlab),
    mlmpi(mlmpi),
    data(data == nullptr ? std::make_shared<Record>() : data),
    basisFunctions(setup::makeBasisFunctions(*incar)),
    parallelMode(get_parallelMode(*incar)),
    randomSeed(input->cget<Vec1Int>("ML_RANDOM_SEED")),
    randomUnique(randomSeed),
    randomUniform(randomSeed),
    structures(dataset::prepareStructures(*mlab)),
    selectionMethod(chooseSelectionMethod(input)),
    selectionMetric(chooseMetric(input)),
    dataDist(DataDist::roundRobin),
    doWrite(!mlmpi or (mlmpi and mlmpi->get_rank() == 0) ? true : false),
    descNormalized(false)
{
    const auto variant = setup::makeFarhtestPointCriterion(input, mlab);
    if (variant.index() == 0) nLocConf = std::get<Vec1Int>(variant);
    else if (variant.index() == 1) typeTresholds = std::get<Vec1Real>(variant);
    this->data->add("descriptors", "Vec1ShRec");

    if (mlmpi != nullptr) mpiRankCheck(dataset::extractNumberAtoms(*mlab), *mlmpi);
}

void Selector::init(const ShRec&                  input,
                    const ShRec&                  mlab,
                    const ShRec&                  data,
                    const std::shared_ptr<MlMPI>& mlmpi)
{
    incar = input;
    this->mlab = mlab;
    this->data = data == nullptr ? std::make_shared<Record>() : data;
    this->mlmpi = mlmpi;
    basisFunctions = setup::makeBasisFunctions(*incar);
    this->parallelMode = get_parallelMode(*incar);
    randomSeed = input->cget<Vec1Int>("ML_RANDOM_SEED");
    randomUnique.init(randomSeed);
    randomUniform.init(randomSeed);
    doWrite = !mlmpi or (mlmpi and mlmpi->get_rank() == 0) ? true : false;
    descNormalized = false;
    structures = dataset::prepareStructures(*mlab);
    const auto variant = setup::makeFarhtestPointCriterion(input, mlab);
    if (variant.index() == 0) nLocConf = std::get<Vec1Int>(variant);
    else if (variant.index() == 1) typeTresholds = std::get<Vec1Real>(variant);
    this->data->add("descriptors", "Vec1ShRec");
    // chose selection method
    selectionMethod = chooseSelectionMethod(input);
    selectionMetric = chooseMetric(input);
    dataDist = DataDist::roundRobin;
    if (mlmpi != nullptr) mpiRankCheck(dataset::extractNumberAtoms(*mlab), *mlmpi);

    return;
}

void Selector::select()
{
    if (selectionMethod == SelectionMethod::None)
        global::tutor.error("Error Selector::select(). No selection method chosen");
    if (doWrite) std::cout << "Start select" << std::endl;
    setUpDataDistribution();
    makeAtomLists();
    if (selectionMethod != SelectionMethod::RandomSelection)
    {
        computeDescriptors();
        collectDistributedData();
    }
    for (Int type = 0; type < mlab->cget<Int>("maxTypes"); type++)
    {
        if (selectionMethod == SelectionMethod::RandomSelection)
        {
            if (doWrite) randomSelection((Size)type);
        }
        else if (selectionMethod == SelectionMethod::FarthestPointSamplingN)
        {
            farthestPointSamplingNPoints((Size)type);
        }
        else if (selectionMethod == SelectionMethod::FarthestPointSamplingThreshold)
        {
            farthestPointSamplingThreshold((Size)type);
        }
        else if (selectionMethod == SelectionMethod::SequentialLikelihoodCover)
        {
            sequentialLikelihoodCover((Size)type);
        }
    }
    if (doWrite) std::cout << "Done select" << std::endl;

    return;
}

void Selector::computeDescriptors()
{
    if (doWrite) std::cout << " -Start computeDescriptors" << std::endl;
    VASPML_PROFILING_START("Compute Descriptors");
    const Vec1String& typeOrder = mlab->get<Vec1String>("typeNames");
    Frame             frame;
    frame.initRefit(*incar, typeOrder, basisFunctions);
    data->get<Vec1ShRec>("descriptors").resize(structures.size());
    DescriptorCollector collector(1, std::map<String, ShVec2Real>());
    centralIndexes.resize(structures.size());
    for (Size struc = 0; struc < structures.size(); struc++)
    {
        centralIndexes[struc] = choseDataDistribution(structures[struc].get_numAtoms());

        VASPML_PROFILING_START("Frame.update");
        frame.update(structures[struc], centralIndexes[struc]);
        VASPML_PROFILING_STOP("Frame.update");

        const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors =
            frame.get_descriptors();
        const std::shared_ptr<const TypeMap>& typeMap = frame.get_typeMap();
        VASPML_PROFILING_START("collector.updateCollector");
        collector.setDescriptorMap(structureDescriptors);
        collector.updateCollector();
        collector.rearrangeSHS2Body(*typeMap);
        VASPML_PROFILING_STOP("collector.updateCollector");
        VASPML_PROFILING_START("collector.get_descriptorsNormalized");
        const auto& descNormalized = collector.get_descriptorsNormalized();
        data->get<Vec1ShRec>("descriptors")[struc] = std::make_shared<Record>();
        for (const auto& [key, desc] : descNormalized)
        {
            if (structureDescriptors.at(key)->get_weight() > 0)
                data->get<Vec1ShRec>("descriptors")[struc]->put<Vec2Real>(key, *desc);
        }
        VASPML_PROFILING_STOP("collector.get_descriptorsNormalized");
        structures[struc].directToCartesian();
    }
    if (doWrite) std::cout << " -Done computeDescriptors" << std::endl;
    // compute array sizes
    const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors =
        frame.get_descriptors();
    for (const String& key : constants::descriptorKeyList)
    {
        if (structureDescriptors.at(key)->get_weight() > 0)
        {
            descSizes[key] = structureDescriptors.at(key)->get_sizeDescriptor(0);
            descOffset[key] = 0;
            activeDescriptors.push_back(key);
            if (structureDescriptors.at(key)->get_isNormalized()) descNormalized = true;
            else descNormalized = false;
        }
    }
    VASPML_PROFILING_STOP("Compute Descriptors");

    return;
}

void Selector::resizeTypeOrderDesc(const Size typeIndex)
{
    Size nAtomsType = 0;
    for (Size i = 0; i < structures.size(); i++)
    {
        const Structure&                      struc = structures[i];
        const std::shared_ptr<const TypeMap>& typeMap = struc.get_typeMap();
        const Int                             strucType = typeMap->toSubType(typeIndex);
        if (strucType < 0) continue;
        nAtomsType += struc.get_numAtomsPerType()[strucType];
    }
    typeOrderDesc.resize(nAtomsType);
    Size nDesc = 0;
    for (Size i = 0; i < activeDescriptors.size(); ++i)
    {
        nDesc += descSizes.at(activeDescriptors[i]);
    }
    std::for_each(typeOrderDesc.begin(),
                  typeOrderDesc.end(),
                  [&](Vec1Real& v) { v.resize(nDesc); });

    return;
}

void Selector::reorderDescriptors(const Size typeIndex)
{
    if (doWrite) std::cout << " -Start reorderDescriptors " << typeIndex << std::endl;
    resizeTypeOrderDesc(typeIndex);
    Size totAtom = 0;
    for (Size struc = 0; struc < structures.size(); struc++)
    {
        const Vec1Int& atomsPerType = structures[struc].get_numAtomsPerType();
        const TypeMap& typeMap = *structures[struc].get_typeMap();
        const Int      subType = typeMap.toSubType(typeIndex);
        if (subType < 0) continue;
        const Record& descriptor = *data->get<Vec1ShRec>("descriptors")[struc];
        const Vec1Int offset = math::partialSum0(atomsPerType);
        for (Int atom = 0; atom < atomsPerType[subType]; atom++)
        {
            Size nDesc = 0;
            for (const String& key : descriptor.keys())
            {
                const Vec2Real& desc = descriptor.cget<Vec2Real>(key);
                for (Int i = 0; i < descSizes.at(key); i++)
                {
                    typeOrderDesc[totAtom][nDesc] = desc[offset[subType] + atom][i];
                    nDesc++;
                }
            }
            totAtom++;
        }
    }
    if (doWrite) std::cout << " -Done reorderDescriptors" << std::endl;

    return;
}

void Selector::randomSelection(const Size type)
{
    makeRandomSelection(type);
    addListToMlab(type);

    return;
}

void Selector::makeRandomSelection(const Size type)
{
    std::cout << " -Start makeRandomSelection: Selecting " << nLocConf[type] << " samples "
              << std::endl;
    VASPML_PROFILING_START("makeRandomSelection");
    randomUnique.reinitializeRange(0, strucIndex[type].size());
    indexSelection = randomUnique.draw(nLocConf[type]);
    VASPML_PROFILING_STOP("makeRandomSelection");
    std::cout << " -Done makeRandomSelection " << std::endl;

    return;
}

void Selector::farthestPointSamplingNPoints(const Size type)
{
    if (parallelMode == ParallelMode::OpenMP)
    {
        if (mlmpi == nullptr or (mlmpi != nullptr and mlmpi->get_rank() == 0))
        {
            reorderDescriptors(type);
            makeFarthestPointSamplingNPointsOMP(type);
        }
    }
    else if (parallelMode == ParallelMode::MPI)
    {
        if (mlmpi) makeFarthestPointSamplingNPointsMPI(type);
        else
            global::tutor.bug(
                "ERROR: " + flf(VASPML_FLF)
                + " You are trying to use farthestPointSamplingKPointsMPI without \n"
                  "giving MPI communicator to your Selector class. Try using OMP version instead.");
    }
    addListToMlab(type);

    return;
}

void Selector::makeFarthestPointSamplingNPointsOMP(const Size type)
{
    if (doWrite)
        std::cout << " -Start makeFarthestPointSamplingOMP: Selecting " << nLocConf[type]
                  << " samples " << std::endl;
    VASPML_PROFILING_START("farthestPointSamplingNPointsOMP");
    const Size seed = randomUniform.uniformInt(1, std::numeric_limits<Int>::max());
    if (selectionMetric == SelectionMetric::L1Norm)
        indexSelection = fps::farthestPointSamplingKPointsOMP<Real>(typeOrderDesc,
                                                                    nLocConf[type],
                                                                    MetricFunctions<Real>::l1Norm,
                                                                    seed);
    else if (selectionMetric == SelectionMetric::L2Norm)
        if (descNormalized)
            indexSelection =
                fps::farthestPointSamplingKPointsOMP<Real>(typeOrderDesc,
                                                           nLocConf[type],
                                                           MetricFunctions<Real>::l2NormUnity,
                                                           seed);
        else
            indexSelection =
                fps::farthestPointSamplingKPointsOMP<Real>(typeOrderDesc,
                                                           nLocConf[type],
                                                           MetricFunctions<Real>::l2Norm,
                                                           seed);
    else if (selectionMetric == SelectionMetric::PolyKernelNorm)
    {
        const Int power = incar->cget<Int>("ML_NHYP");
        if (power == 4)
            indexSelection =
                fps::farthestPointSamplingKPointsOMP<Real>(typeOrderDesc,
                                                           nLocConf[type],
                                                           MetricFunctions<Real>::polyKernelNorm<4>,
                                                           seed);
        else global::tutor.error("ERROR: " + flf(VASPML_FLF) + " only available for ML_NHYP=4.");
    }
    else if (selectionMetric == SelectionMetric::antiPearsonNorm2)
        indexSelection =
            fps::farthestPointSamplingKPointsOMP<Real>(typeOrderDesc,
                                                       nLocConf[type],
                                                       MetricFunctions<Real>::antiPearsonNorm2,
                                                       seed);
    VASPML_PROFILING_STOP("farthestPointSamplingNPointsOMP");
    if (doWrite) std::cout << " -Done makeFarthestPointSampling " << std::endl;

    return;
}

void Selector::makeFarthestPointSamplingNPointsMPI(const Size type)
{
    if (doWrite)
        std::cout << " -Start makeFarthestPointSampling: Selecting " << nLocConf[type]
                  << " samples " << std::endl;
    VASPML_PROFILING_START("farthestPointSamplingNPointsMPI");
    const Size seed = randomUniform.uniformInt(1, std::numeric_limits<Int>::max());
    if (selectionMetric == SelectionMetric::L1Norm)
        indexSelection =
            fps::farthestPointSamplingKPointsMPI<Real>(*descriptorsShmem[type],
                                                       nLocConf[type],
                                                       MetricFunctions<Real>::l1NormPtr,
                                                       mlmpi,
                                                       seed);
    else if (selectionMetric == SelectionMetric::L2Norm)
        if (descNormalized)
            indexSelection =
                fps::farthestPointSamplingKPointsMPI<Real>(*descriptorsShmem[type],
                                                           nLocConf[type],
                                                           MetricFunctions<Real>::l2NormPtrUnity,
                                                           mlmpi,
                                                           seed);
        else
            indexSelection =
                fps::farthestPointSamplingKPointsMPI<Real>(*descriptorsShmem[type],
                                                           nLocConf[type],
                                                           MetricFunctions<Real>::l2NormPtr,
                                                           mlmpi,
                                                           seed);
    else if (selectionMetric == SelectionMetric::PolyKernelNorm)
    {
        const Int power = incar->cget<Int>("ML_NHYP");
        if (power == 4)
            indexSelection = fps::farthestPointSamplingKPointsMPI<Real>(
                *descriptorsShmem[type],
                nLocConf[type],
                MetricFunctions<Real>::polyKernelNormPtr<4>,
                mlmpi,
                seed);
        else global::tutor.error("ERROR: " + flf(VASPML_FLF) + " only available for ML_NHYP=4.");
    }
    else if (selectionMetric == SelectionMetric::antiPearsonNorm2)
        indexSelection =
            fps::farthestPointSamplingKPointsMPI<Real>(*descriptorsShmem[type],
                                                       nLocConf[type],
                                                       MetricFunctions<Real>::antiPearsonNorm2Ptr,
                                                       mlmpi,
                                                       seed);
    VASPML_PROFILING_STOP("farthestPointSamplingNPointsMPI");
    if (doWrite) std::cout << " -Done makeFarthestPointSampling " << std::endl;

    return;
}

void Selector::farthestPointSamplingThreshold(const Size type)
{
    if (parallelMode == ParallelMode::OpenMP)
    {
        if (mlmpi == nullptr or (mlmpi != nullptr and mlmpi->get_rank() == 0))
        {
            reorderDescriptors(type);
            makeFarthestPointSamplingThresholdOMP(type);
        }
    }
    else if (parallelMode == ParallelMode::MPI) { makeFarthestPointSamplingThresholdMPI(type); }
    addListToMlab(type);

    return;
}

void Selector::makeFarthestPointSamplingThresholdOMP(const Size type)
{
    const Real tresh = typeTresholds[type];
    if (doWrite)
        std::cout << " -Start makeFarthestPointSamplingTreshold: Selecting with treshold " << tresh
                  << std::endl;
    const Size seed = randomUniform.uniformInt(1, std::numeric_limits<Int>::max());
    VASPML_PROFILING_START("makeFarthestPointSamplingThresholdOMP");
    if (selectionMetric == SelectionMetric::L1Norm)
        indexSelection = fps::farthestPointSamplingDistTreshOMP<Real>(typeOrderDesc,
                                                                      tresh,
                                                                      MetricFunctions<Real>::l1Norm,
                                                                      seed);
    else if (selectionMetric == SelectionMetric::L2Norm)
        if (descNormalized)
            indexSelection =
                fps::farthestPointSamplingDistTreshOMP<Real>(typeOrderDesc,
                                                             tresh,
                                                             MetricFunctions<Real>::l2NormUnity,
                                                             seed);
        else
            indexSelection =
                fps::farthestPointSamplingDistTreshOMP<Real>(typeOrderDesc,
                                                             tresh,
                                                             MetricFunctions<Real>::l2Norm,
                                                             seed);
    else if (selectionMetric == SelectionMetric::PolyKernelNorm)
    {
        const Int power = incar->cget<Int>("ML_NHYP");
        if (power == 4)
            indexSelection = fps::farthestPointSamplingDistTreshOMP<Real>(
                typeOrderDesc,
                tresh,
                MetricFunctions<Real>::polyKernelNorm<4>,
                seed);
        else global::tutor.error("ERROR: " + flf(VASPML_FLF) + " only available for ML_NHYP=4.");
    }
    else if (selectionMetric == SelectionMetric::antiPearsonNorm2)
        indexSelection =
            fps::farthestPointSamplingDistTreshOMP<Real>(typeOrderDesc,
                                                         tresh,
                                                         MetricFunctions<Real>::antiPearsonNorm2,
                                                         seed);
    VASPML_PROFILING_STOP("makeFarthestPointSamplingThresholdOMP");
    if (doWrite) std::cout << " -Done makeFarthestPointSamplingTreshold " << std::endl;

    return;
}

void Selector::makeFarthestPointSamplingThresholdMPI(const Size type)
{
    const Real tresh = typeTresholds[type];
    if (doWrite)
        std::cout << " -Start makeFarthestPointSamplingTresholdMPI: Selecting with treshold "
                  << tresh << std::endl;
    const Size seed = randomUniform.uniformInt(1, std::numeric_limits<Int>::max());
    VASPML_PROFILING_START("makeFarthestPointSamplingThresholdMPI");
    if (selectionMetric == SelectionMetric::L1Norm)
        indexSelection =
            fps::farthestPointSamplingDistTreshMPI<Real>(*descriptorsShmem[type],
                                                         tresh,
                                                         MetricFunctions<Real>::l1NormPtr,
                                                         mlmpi,
                                                         seed);
    else if (selectionMetric == SelectionMetric::L2Norm)
        if (descNormalized)
            indexSelection =
                fps::farthestPointSamplingDistTreshMPI<Real>(*descriptorsShmem[type],
                                                             tresh,
                                                             MetricFunctions<Real>::l2NormPtrUnity,
                                                             mlmpi,
                                                             seed);
        else
            indexSelection =
                fps::farthestPointSamplingDistTreshMPI<Real>(*descriptorsShmem[type],
                                                             tresh,
                                                             MetricFunctions<Real>::l2NormPtr,
                                                             mlmpi,
                                                             seed);
    else if (selectionMetric == SelectionMetric::PolyKernelNorm)
    {
        const Int power = incar->cget<Int>("ML_NHYP");
        if (power == 4)
            indexSelection = fps::farthestPointSamplingDistTreshMPI<Real>(
                *descriptorsShmem[type],
                tresh,
                MetricFunctions<Real>::polyKernelNormPtr<4>,
                mlmpi,
                seed);
        else global::tutor.error("ERROR: " + flf(VASPML_FLF) + " only available for ML_NHYP=4.");
    }
    else if (selectionMetric == SelectionMetric::antiPearsonNorm2)
        indexSelection =
            fps::farthestPointSamplingDistTreshMPI<Real>(*descriptorsShmem[type],
                                                         tresh,
                                                         MetricFunctions<Real>::antiPearsonNorm2Ptr,
                                                         mlmpi,
                                                         seed);
    VASPML_PROFILING_STOP("makeFarthestPointSamplingThresholdMPI");
    if (doWrite) std::cout << " -Done makeFarthestPointSamplingTresholdMPI " << std::endl;

    return;
}

void Selector::sequentialLikelihoodCover(const Size type)
{
    if (parallelMode == ParallelMode::OpenMP)
    {
        if (mlmpi == nullptr or (mlmpi != nullptr and mlmpi->get_rank() == 0))
        {
            reorderDescriptors(type);
            makeSequentialLikelihoodCover(type);
        }
    }
    else if (parallelMode == ParallelMode::MPI)
    {
        std::cout << "NOT IMPLEMENTED YET !!!!!" << std::endl;
        throw;
    }
    addListToMlab(type);

    return;
}

void Selector::makeSequentialLikelihoodCover(const Size type)
{
    if (doWrite)
        std::cout << " -Start makeSequentialLikelihoodCover with maximum number of samples "
                  << nLocConf[type] << std::endl;
    VASPML_PROFILING_START("makeSequentialLikelihoodCoverOMP");
    if (selectionMetric == SelectionMetric::L1Norm)
    {
        indexSelection =
            bayes_select::sequentialLikelihoodCoverageSelectOMP<Real>(typeOrderDesc,
                                                                      MetricFunctions<Real>::l1Norm,
                                                                      nLocConf[type]);
    }
    else if (selectionMetric == SelectionMetric::L2Norm)
    {
        if (descNormalized)
        {
            indexSelection = bayes_select::sequentialLikelihoodCoverageSelect<Real>(
                typeOrderDesc,
                MetricFunctions<Real>::l2NormUnity,
                nLocConf[type],
                0.0);
        }
        else
        {
            indexSelection = bayes_select::sequentialLikelihoodCoverageSelect<Real>(
                typeOrderDesc,
                MetricFunctions<Real>::l2Norm,
                nLocConf[type]);
        }
    }
    else if (selectionMetric == SelectionMetric::PolyKernelNorm)
    {
        const Int power = incar->cget<Int>("ML_NHYP");
        if (power == 4)
            indexSelection = bayes_select::sequentialLikelihoodCoverageSelect<Real>(
                typeOrderDesc,
                MetricFunctions<Real>::polyKernelNorm<4>,
                nLocConf[type]);
        else global::tutor.error("ERROR: " + flf(VASPML_FLF) + " only available for ML_NHYP=4.");
    }
    else if (selectionMetric == SelectionMetric::antiPearsonNorm2)
    {
        indexSelection = bayes_select::sequentialLikelihoodCoverageSelect<Real>(
            typeOrderDesc,
            MetricFunctions<Real>::antiPearsonNorm2,
            nLocConf[type]);
    }
    VASPML_PROFILING_STOP("makeSequentialLikelihoodCoverOMP");
    if (doWrite) std::cout << " -Done makeSequentialLikelihoodCover " << std::endl;

    return;
}

void Selector::resizeIndexArrays()
{
    Int nTypes = mlab->cget<Int>("maxTypes");
    strucIndex.resize(nTypes);
    atomIndexInStruc.resize(nTypes);

    return;
}

void Selector::makeAtomLists()
{
    resizeIndexArrays();
    for (Size struc = 0; struc < structures.size(); struc++)
    {
        const Vec1Int& atomsPerType = structures[struc].get_numAtomsPerType();
        const TypeMap& typeMap = *structures[struc].get_typeMap();
        Int            atomOffset = 0;
        for (Size type = 0; type < atomsPerType.size(); type++)
        {
            Int superType = typeMap.toType(type);
            for (Int atom = 0; atom < atomsPerType[type]; atom++)
            {
                strucIndex[superType].push_back(struc);
                atomIndexInStruc[superType].push_back(atom + atomOffset);
            }
            atomOffset += atomsPerType[type];
        }
    }

    return;
}

void Selector::addListToMlab(const Size typeIndex)
{
    if (parallelMode == ParallelMode::OpenMP and mlmpi and mlmpi->get_rank() != 0) return;
    VASPML_PROFILING_START("addListToMlab");
    mlab->get<Vec1Int>("numLrc")[typeIndex] = indexSelection.size();
    Vec1Int& lrcStructure = mlab->get<Vec2Int>("lrcStructure")[typeIndex];
    Vec1Int& lrcAtom = mlab->get<Vec2Int>("lrcAtom")[typeIndex];
    lrcStructure.resize(indexSelection.size());
    lrcAtom.resize(indexSelection.size());
    for (Size i = 0; i < indexSelection.size(); i++)
    {
        const Int idx = indexSelection[i];
        lrcStructure[i] = strucIndex[typeIndex][idx] + 1;
        lrcAtom[i] = atomIndexInStruc[typeIndex][idx] + 1;
    }
    VASPML_PROFILING_STOP("addListToMlab");

    return;
}

void Selector::writeMlabFile(const String& filename) const
{
    if (mlmpi->get_rank() == 0) io::writeMlab(*mlab, filename);

    return;
}

Vec1Int Selector::choseDataDistribution(const Int& numberElements)
{
    if (mlmpi and dataDist == DataDist::roundRobin)
        return mlmpi::getRoundRobinIndexes(numberElements,
                                           mlmpi->get_numberRanks(),
                                           mlmpi->get_rank());
    else if (mlmpi and dataDist == DataDist::block)
        return mlmpi::getBlockDistributionIndexes(numberElements,
                                                  mlmpi->get_numberRanks(),
                                                  mlmpi->get_rank());
    else
    {
        Vec1Int index(numberElements);
        std::iota(index.begin(), index.end(), 0);
        return index;
    }
}

void Selector::setUpDataDistribution()
{
    centralIndexes.resize(structures.size());
    for (Size struc = 0; struc < structures.size(); struc++)
    {
        centralIndexes[struc] = choseDataDistribution(structures[struc].get_numAtoms());
    }

    return;
}

void Selector::collectDistributedData()
{
    if (parallelMode == ParallelMode::OpenMP) collectDistributedDataOMP();
    else if (parallelMode == ParallelMode::MPI) collectDistributedDataShmemArray();

    return;
}

void Selector::collectDistributedDataOMP()
{
    VASPML_PROFILING_START("collectDistributedData");
    if (mlmpi != nullptr)
    {
        for (Size struc = 0; struc < structures.size(); struc++)
        {
            ShRec& descShRec = data->get<Vec1ShRec>("descriptors")[struc];
            ShRec  collectIndex = std::make_shared<Record>();
            collectIndex->put("index", centralIndexes[struc]);
            collectIndex = rec::mpi::gatherV(*collectIndex, *mlmpi);
            ShRec descShRecTot = rec::mpi::gatherV(*descShRec, *mlmpi);
            if (mlmpi->get_rank() == 0)
            {
                descShRec = descShRecTot;
                Vec1Int newOrder = vector_tools::invertIndex(collectIndex->cget<Vec1Int>("index"));
                for (const String& key : descShRec->keys())
                {
                    // undo round robin
                    vector_tools::inplaceReorderVector(descShRec->get<Vec2Real>(key),
                                                       newOrder); // only Vec2Real
                }
            }
        }
    }
    VASPML_PROFILING_STOP("collectDistributedData");

    return;
}

void Selector::collectDistributedDataShmemArray()
{
    VASPML_PROFILING_START("collectDistributedDataShmemArray");
    if (mlmpi != nullptr)
    {
        for (Size struc = 0; struc < structures.size(); struc++)
        {
            ShRec& descShRec = data->get<Vec1ShRec>("descriptors")[struc];
            ShRec  collectIndex = std::make_shared<Record>();
            collectIndex->put("index", centralIndexes[struc]);
#ifdef use_shmem
            collectIndex = rec::mpi::gatherV(*collectIndex, *mlmpi);
            ShRec descShRecTot = rec::mpi::gatherV(*descShRec, *mlmpi);
#else
            collectIndex = rec::mpi::allGatherV(*collectIndex, *mlmpi);
            ShRec descShRecTot = rec::mpi::allGatherV(*descShRec, *mlmpi);
#endif
#ifdef use_shmem
            if (mlmpi->get_rank() == 0)
            {
#endif
                descShRec = descShRecTot;
                Vec1Int newOrder = vector_tools::invertIndex(collectIndex->cget<Vec1Int>("index"));
                for (const String& key : descShRec->keys())
                {
                    // undo round robin
                    vector_tools::inplaceReorderVector(descShRec->get<Vec2Real>(key),
                                                       newOrder); // only Vec2Real
                }
#ifdef use_shmem
            }
#endif
        }
        allocate_descriptorsShmem();
        fill_descriptorsShmem();
    }
    VASPML_PROFILING_STOP("collectDistributedDataShmemArray");

    return;
}

void Selector::allocate_descriptorsShmem()
{
    if (doWrite) std::cout << "-Start allocate_descriptorsShmem" << std::endl;

    Vec1Int totalNumberAtoms = dataset::countTotalAtomsPerType(*mlab);
    descriptorsShmem.resize(totalNumberAtoms.size());
    Size nDesc = 0;
    for (Size i = 0; i < activeDescriptors.size(); ++i)
    {
        nDesc += descSizes.at(activeDescriptors[i]);
    }
    for (Size type = 0; type < totalNumberAtoms.size(); type++)
    {
        descriptorsShmem[type] =
            std::make_shared<ShmemArray2D<Real>>(totalNumberAtoms[type], nDesc, mlmpi);
    }
    if (doWrite) std::cout << "-Done allocate_descriptorsShmem" << std::endl;

    return;
}

void Selector::fill_descriptorsShmem()
{
    if (doWrite) std::cout << "-Start fill_descriptorsShmem" << std::endl;
    Vec1Size totAtom(mlab->get<Int>("maxTypes"), 0);
    for (Size struc = 0; struc < structures.size(); struc++)
    {
        ShRec& descShRec = data->get<Vec1ShRec>("descriptors")[struc];
#ifdef use_shmem
        if (mlmpi and mlmpi->get_rank() == 0)
        {
#endif
            const Vec1Int& atomsPerType = structures[struc].get_numAtomsPerType();
            const TypeMap& typeMap = *structures[struc].get_typeMap();
            const Record&  descriptor = *descShRec;
            for (Size subType = 0; subType < typeMap.countStructureTypes(); subType++)
            {
                const Size    typeIndex = typeMap.toType(subType);
                const Vec1Int offset = math::partialSum0(atomsPerType);
                for (Int atom = 0; atom < atomsPerType[subType]; atom++)
                {
                    Size nDesc = 0;
                    for (const String& key : descriptor.keys())
                    {
                        const Vec2Real& desc = descriptor.cget<Vec2Real>(key);
                        for (Int i = 0; i < descSizes.at(key); i++)
                        {
                            descriptorsShmem[typeIndex]->set_value(totAtom[typeIndex],
                                                                   nDesc,
                                                                   desc[offset[subType] + atom][i]);
                            nDesc++;
                        }
                    }
                    totAtom[typeIndex]++;
                }
#ifdef use_shmem
            }
#endif
        }
        descShRec.reset();
    }
    if (doWrite) std::cout << "-Done fill_descriptorsShmem" << std::endl;
#ifdef use_shmem
    if (mlmpi) mlmpi->barrier();
#endif

    return;
}

const ShCRec Selector::get_mlab() const
{
    return mlab;
}
