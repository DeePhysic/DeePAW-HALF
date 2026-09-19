#include "Refit.hpp"

#include "BatchTools.hpp"
#include "Descriptor.hpp"
#include "DescriptorSHS3.hpp"
#include "Linalg.hpp"
#include "Random.hpp"
#include "Record.hpp"
#include "Tutor.hpp"
#include "dataset.hpp"
#include "debug.hpp"
#include "io.hpp"
#include "rec_shmem.hpp"
#include "setup.hpp"
#include "utils.hpp"

#include <algorithm> // for fill, for_each
#include <cmath>     // for ceil
#include <fstream>   // for basic_ofstream, basic_fstream
#include <iostream>  // for basic_ostream, operator<<, cout, endl

namespace vaspml
{
class DescriptorCollector;
}

using namespace vaspml;

void vaspml::fillRegCoefficients(Vec1Real& regCoeff, const Int n)
{
    Int m = 1;
    regCoeff.resize(n);
    if (m == 0) { std::fill(regCoeff.begin(), regCoeff.end(), 0); }
    else if (m == 1) { std::fill(regCoeff.begin(), regCoeff.end(), 1); }
    else if (m == 2)
    {
        Random<Real> rand(12091810);
        std::for_each(regCoeff.begin(),
                      regCoeff.end(),
                      [&](Real& x) { x = 10 * (rand.uniformFloat() + 1); });
    }

    return;
}

std::vector<Frame> vaspml::makeBatchFrameStorage(const Int&              nFrames,
                                                 const Record&           incar,
                                                 const Record&           mlab,
                                                 const BasisFunctionMap& basisFunctions)
{
    std::vector<Frame> frames(nFrames);
    std::for_each(frames.begin(),
                  frames.end(),
                  [&](Frame& frame)
                  { frame.initRefit(incar, mlab.cget<Vec1String>("typeNames"), basisFunctions); });

    return frames;
}

std::vector<AtomBatchMap> vaspml::initAtomBatchMaps(const std::vector<Structure>& structures,
                                                    const Size                    nFrames,
                                                    const Size                    batchSize)
{
    std::vector<AtomBatchMap> atomBatchMap(nFrames);
    Size                      nStrucs = structures.size();
    for (Size batchNum = 0; batchNum < nFrames; batchNum++)
    {
        Size start = batchNum * batchSize;
        Size end = batchSize * (batchNum + 1);
        if (nStrucs < end) { end = nStrucs; }
        atomBatchMap[batchNum].makeMap(structures.begin() + start, structures.begin() + end);
    }

    return atomBatchMap;
}

Refit::Refit(const ShRec& incar, const ShRec& mlab, std::shared_ptr<MlMPI> mpi) :
    mlab(mlab),
    incar(incar),
    mpi(mpi),
    mlff(std::make_shared<Record>()),
    basisFunctions(setup::makeBasisFunctions(*incar)),
    nFrames(std::ceil((Real)mlab->cget<Int>("numStructures") / (Real)batchSize)),
    structures(dataset::prepareStructures(*mlab)),
    frames(makeBatchFrameStorage(nFrames, *incar, *mlab, basisFunctions)),
    nStrucs(structures.size()),
    nRefConfs(vector_tools::sum(mlab->cget<Vec1Int>("numLrc"))),
    nLocRefConfs(mlab->cget<Vec1Int>("numLrc")),
    atomBatchMap(initAtomBatchMaps(structures, nFrames, batchSize)),
    strucToBatch(batchSize),
    kernelPolynomial(*incar, *mlab, nullptr),
    fitMatrix(dataset::extractNumberAtoms(*mlab), nRefConfs),
    energyShift(
        dataset::computeEnergyShiftIsolatedAtom(*mlab, incar->cget<Vec1Real>("ML_EATOM_REF"))),
    fitAlgorithm(incar->cget<String>("ML_RALGO") == "design" ? FitAlgorithm::DesignMatrix
                                                             : FitAlgorithm::NormalEquation)
{
    //dgemmTests();
    switch (fitAlgorithm)
    {
    case FitAlgorithm::DesignMatrix:
        fitMatrix.initDesignMatrixMode();
        std::cout << "Using design matrix approach for fitting" << std::endl;
        break;
    case FitAlgorithm::NormalEquation:
        fitMatrix.initNormalEqMode();
        std::cout << "Using normal equation approach for fitting" << std::endl;
        break;
    default:
        global::tutor.error(
            "Unknown fitting algorithm specified! Please choose between 'design' and 'normal'\n");
    }
}

void Refit::prepareFit()
{
    nStrucs = mlab->cget<Int>("numStructures");
    computeFrames();
    prepare_numberDescriptors();
    prepare_kernelPolynomial();
    computeFitMatrix();

    return;
}

void Refit::computeFrames()
{
    // @NOTE remove DEBUG
    //nFrames = 20;
    for (Int struc = 0; struc < nFrames; struc++)
    {
        frames[struc].update(structures[struc]);
    }

    return;
}

void Refit::prepare_numberDescriptors()
{
    // are the same for every batchedFrame. Because the number of
    // descriptors only depends on the number of basis functions and
    // the total number of types.
    const auto& descriptors = frames[0].get_descriptors();
    const Int&  nTypes = mlab->cget<Int>("maxTypes");
    for (const auto& [key, desc] : descriptors)
    {
        if (desc->get_weight() > 0.0)
        {
            const Int nDesc = desc->get_sizeDescriptor(0);
            numberDescriptors[key] = std::make_shared<Vec1Int>(nTypes, nDesc);
        }
        // currently all types have same number of descriptors.
        // this will always be true when setting up a new refit
    }

    return;
}

void Refit::prepare_kernelPolynomial()
{
    kernelPolynomial.set_numberDescriptors(numberDescriptors);
    kernelPolynomial.generate_featureSpaceMap();
    kernelPolynomial.make_ShmemArrayMap();
    setReferenceConfigurations();

    return;
}

void Refit::setReferenceConfigurations()
{
    const Vec2Int& lrcStructure = mlab->get<Vec2Int>("lrcStructure");
    const Vec2Int& lrcAtom = mlab->get<Vec2Int>("lrcAtom");
    strucToBatch.setNumberOfStructures(mlab->get<Int>("numStructures"));
    std::shared_ptr<DescriptorCollector> referenceConfigurations =
        batch_tools::createDescriptorCollectorReferenceConfigs(frames,
                                                               lrcStructure,
                                                               lrcAtom,
                                                               atomBatchMap,
                                                               strucToBatch);
    kernelPolynomial.set_descriptorsRefConfs(referenceConfigurations);

    return;
}

void Refit::computeFitMatrix()
{
    std::cout << "Starting to computeFitMatrix" << std::endl;
    Vec1ShCRec MLABstrucs = mlab->vcget<Vec1ShRec>("structures");
    switch (fitAlgorithm)
    {
    case FitAlgorithm::DesignMatrix:
        std::cout << "Using design matrix approach for fitting" << std::endl;
        for (Size i = 0; i < frames.size(); i++)
        {
            std::cout << "Update kernel " << i << std::endl;
            kernelPolynomial.updatePolynomialKernelRefit(frames[i]);
            std::cout << "Add structure to design Matrix " << i << std::endl;
            fitMatrix.addStructureDesignMatrix(kernelPolynomial, i);
            //throw;
        }
        break;
    case FitAlgorithm::NormalEquation:
        for (Size i = 0; i < frames.size(); i++)
        {
            std::cout << "Update kernel " << i << std::endl;
            kernelPolynomial.updatePolynomialKernelRefit(frames[i]);
            std::cout << "Add structure to design Matrix " << i << std::endl;
            fitMatrix.addStructureNormalEquation(kernelPolynomial,
                                                 MLABstrucs[i]->cget<Real>("energy"),
                                                 MLABstrucs[i]->cget<Vec1Real>("forces"),
                                                 MLABstrucs[i]->cget<Vec1Real>("stress"));
        }
        break;
    default:
        global::tutor.error(
            "Unknown fitting algorithm specified! Please choose between 'design' and 'normal'\n");
    }
    fitMatrix.write_fitMatrix2D("DesignMatrix.dat");

    return;
}

void Refit::computeFit()
{
#ifndef VASPML_USE_FIT_DEBUG
    switch (fitAlgorithm)
    {
    case FitAlgorithm::DesignMatrix:
        std::cout << "Using design matrix approach for fitting" << std::endl;
        SVDDesignMatrixFitting();
        break;
    case FitAlgorithm::NormalEquation:
        std::cout << "Using normal equation approach for fitting" << std::endl;
        SVDNormalequationFitting();
        break;
    default:
        global::tutor.error(
            "Unknown fitting algorithm specified! Please choose between 'design' and 'normal'\n");
    }
#else
    fillRegCoefficients(regressionCoefficients, (Int)fitMatrix.get_fitMatrixSize1());
#endif

    return;
}

void Refit::SVDDesignMatrixFitting()
{
    //Vec1Real target = dataset::prepareEnergyRefitIsolatedAtom( *mlab, incar->cget<Vec1Real>(
    //"ML_EATOM_REF" ) );
#ifdef VASPML_USE_FIT_DEBUG
    Vec1Real target = dataset::extractEnergies(*mlab);
#else
    Vec1Real target =
        dataset::prepareEnergyRefitIsolatedAtom(*mlab, incar->cget<Vec1Real>("ML_EATOM_REF"));
#endif
    Vec1Real forces = vector_tools::Vec2ToVec1(dataset::extractForces(*mlab));
    target.insert(target.end(), forces.cbegin(), forces.cend());
    Size nTarget = target.size();
    regressionCoefficients.resize(nRefConfs, 0.0);
    const Real regParam = incar->cget<Real>("ML_SIGW0") / incar->cget<Real>("ML_SIGV0");
    std::cout << "regression parameter " << regParam << "  " << nTarget << std::endl;

    auto file = file_io::openFileO("Target.dat");
    for (Size i = 0; i < nTarget; i++) { file << str("%24.16E ", target[i]) << std::endl; }
    file.close();

#ifdef VASPML_USE_FIT_DEBUG
    // copy matrix here to a Vec1Real
    const Real* matrix = fitMatrix.get_fitMatrixPtr();
    Vec1Real    matCopy(fitMatrix.get_fitMatrixSize0() * fitMatrix.get_fitMatrixSize1());
    std::copy(matrix,
              matrix + fitMatrix.get_fitMatrixSize0() * fitMatrix.get_fitMatrixSize1(),
              matCopy.begin());
    Real* matCopyPtr = matCopy.data();
    linalg::solveLeastSquaresSVD(matCopyPtr,
                                 nTarget,
                                 nRefConfs,
                                 target,
                                 regressionCoefficients,
                                 regParam);
#else
    Real* matrix = fitMatrix.get_fitMatrixPtr();
    linalg::solveLeastSquaresSVD(matrix,
                                 nTarget,
                                 nRefConfs,
                                 target,
                                 regressionCoefficients,
                                 regParam);
#endif
    writeRegressionCoefficients("RegressionCoefficients.dat");

    return;
}

void Refit::SVDNormalequationFitting()
{
    std::cout << "Starting routine SVDDesignMatrixFitting" << std::endl;
    Real*           matrix = fitMatrix.get_fitMatrixPtr();
    const Vec1Real& target = fitMatrix.get_normalEquationTarget();
    regressionCoefficients.resize(nRefConfs, 0.0);
    const Real regParam = incar->cget<Real>("ML_SIGW0") / incar->cget<Real>("ML_SIGV0");
    linalg::solveLeastSquaresSVD(matrix,
                                 nRefConfs,
                                 nRefConfs,
                                 target,
                                 regressionCoefficients,
                                 regParam);
    writeRegressionCoefficients("RegressionCoefficients.dat");

    return;
}

void Refit::writeRegressionCoefficients(const String& filename) const
{
    auto file = file_io::openFileO(filename);
    for (const auto& x : regressionCoefficients) { file << str("%24.16E ", x) << std::endl; }

    return;
}

void Refit::prepare_mlffRecord()
{
    dataset::createMlffFromSetup(*incar, *mlff);

    // add version number
    mlff->put<Int>("versionMajor", 0);
    mlff->put<Int>("versionMinor", 2);
    mlff->put<Int>("versionPatch", 1);

    Vec2Real regCoeff = vector_tools::Vec1ToVec2(regressionCoefficients, nLocRefConfs);
    mlff->put<Vec2Real>("regression-coeff", regCoeff);
    // TODO put real rmse values between fit and reference
    mlff->put<Real>("rmseEnergy", 0.0);
    mlff->put<Real>("rmseForces", 0.0);
    mlff->put<Real>("rmseStress", 0.0);

    mlff->put<bool>("ML_LSIC", false);

    std::map<String, Real> weights = kernelPolynomial.get_weights();

    // adding descriptor list
    if (weights["SHS3-3-body"] > 0)
    {
        const Vec2Int& descriptorList =
            static_cast<const DescriptorSHS3*>(
                frames[0].get_descriptors().at("SHS3-3-body").get())
                ->get_descriptorList();
        // NOTE ML_FF file starts counting at 1 for descriptorList
        mlff->put<Vec2Int>("SHS3-3-body-descriptor-list", descriptorList);
        Vec2Int& descListManip = mlff->get<Vec2Int>("SHS3-3-body-descriptor-list");
        for (auto& x : descListManip)
            for (auto& y : x) y += 1;
    }

    std::map<String, std::shared_ptr<Vec1Int>> descSizes = kernelPolynomial.get_numberDescriptors();

    if (weights["SHS3-3-body"] > 0)
        mlff->put<Vec1Int>("SHS3-3-body-number-descriptors-per-type", *descSizes.at("SHS3-3-body"));

    if (weights["SHS2-2-body"] > 0)
        mlff->put<Int>("SHS2-2-body-number-descriptors-per-type",
                       (*descSizes.at("SHS2-2-body"))[0]);

    const auto& refConfs = kernelPolynomial.get_descriptorsRefConfs();
    for (const auto& [key, desc] : refConfs)
    {
        rec::copyAddShmemArray(key + "-reference-configs", *desc, *mlff);
    }
    // set number of zeros of the bessel functions
    const Vec1Size& nRoots = basisFunctions.at("3-body")->get_nValidRoots();
    mlff->put<Vec1Int>("SHS3-3-body-number-radial-basis", Vec1Int(nRoots.cbegin(), nRoots.cend()));

    // add mlab things
    dataset::addMlabToMlff(*mlab, *mlff);

    dataset::mlabStatistics(*mlff, *mlab);
    mlff->put<Real>("sigv", incar->cget<Real>("ML_SIGV0"));
    mlff->put<Real>("sigw", incar->cget<Real>("ML_SIGW0"));
    mlff->put<Real>(
        "average-energy-per-atom",
        dataset::computeEnergyShiftIsolatedAtom(*mlab, incar->cget<Vec1Real>("ML_EATOM_REF")));
    mlff->put<Vec2Real>( "inverse-kernel-matrix", 
                         Vec2Real( mlab->cget<Int>("maxTypes"), Vec1Real( 1, 0.0 ) ));
    mlffPrepared = true;

    return;
}

void Refit::write_mlffFile(const String& fileName)
{
    if (mlffPrepared)
    {
        std::fstream file;
        io::open(file, fileName, "wb");
        io::processMlff(*mlff, file, true);
        io::close(file);
    }
    else
    {
        global::tutor.error(
            "MLFF not prepared, cannot write to file. Please call prepare_mlffRecord() first!\n");
    }

    return;
}

void Refit::trainDataInference()
{
    if (!mlffPrepared)
    {
        global::tutor.error("MLFF not prepared, cannot do inference on train set data. "
                            "Please call prepare_mlffRecord() first!\n");
        return;
    }
    inference.init(*mlff);
    inference.predict(structures);
    Vec1ShCRec data = inference.getPredictions()->vcget<Vec1ShRec>("predictions");
    Vec1ShCRec refConfs = mlab->vcget<Vec1ShRec>("structures");
    auto [energy, force, stress] = dataset::rmseDatsets(data, refConfs);

    auto file = file_io::openFileO("EnergyPred.dat");
    for (Size i = 0; i < data.size(); i++) { file << data[i]->cget<Real>("energy") << std::endl; }
    file.close();

    std::cout << "RMSE Energy " << energy << std::endl;
    std::cout << "RMSE Force " << force << std::endl;
    std::cout << "RMSE Stress " << stress << std::endl;
}

void Refit::designMatrixTest()
{
#ifdef VASPML_USE_FIT_DEBUG
    if (!mlffPrepared)
    {
        global::tutor.error("MLFF not prepared, cannot do inference on train set data. "
                            "Please call prepare_mlffRecord() first!\n");
    }

    Vec1Real target = dataset::extractEnergies(*mlab);
    Size     nEnergy = target.size();
    Vec1Real forces = vector_tools::Vec2ToVec1(dataset::extractForces(*mlab));
    target.insert(target.end(), forces.cbegin(), forces.cend());
    Size nTarget = target.size();

    inference.init(*mlff);
    inference.predict(structures);

    Vec1ShCRec data = inference.getPredictions()->vcget<Vec1ShRec>("predictions");
    Vec1ShCRec refConfs = mlab->vcget<Vec1ShRec>("structures");

    auto [energy, force, stress] = dataset::rmseDatsets(data, refConfs);

    Vec1Real    predicted(nTarget, (Real)0.0);
    const Real* designMatrix = fitMatrix.get_fitMatrixPtr();
    linalg::matVec(linalg::Transpose::TRANS,
                   predicted.size(),
                   fitMatrix.get_fitMatrixSize1(),
                   1.0,
                   designMatrix,
                   fitMatrix.get_fitMatrixSize1(),
                   regressionCoefficients,
                   1,
                   1.0,
                   predicted,
                   1);

    for (Size i = 0; i < nEnergy; i++)
    {
        std::cout << "Energy " << i << "   " << str(" Real value = %24.16E", target[i]) << "   "
                  << str(" refit = %24.16E", predicted[i]) << "   "
                  << str(" inference = %24.16E",
                         data[i]->cget<Real>("energy") / (Real)structures[i].get_numAtoms())
                  << "    " << std::endl;
    }
    auto file = file_io::openFileO("EnergyComp.dat");
    for (Size i = 0; i < nEnergy; i++)
    {
        file << " " << i << "   " << str("%24.16E", target[i]) << "   "
             << str("%24.16E", predicted[i]) << "   "
             << str("%24.16E", data[i]->cget<Real>("energy") / (Real)structures[i].get_numAtoms())
             << "    " << std::endl;
    }
    file.close();

    Vec1Real predForce(nTarget);
    Int      col = 0;
    for (Int i = 0; i < nEnergy; i++)
    {
        const Vec1Real& force = data[i]->cget<Vec1Real>("forces");
        for (Int j = 0; j < force.size(); j++)
        {
            predForce[col + nEnergy] = force[j];
            col++;
        }
    }

    file.open("ForceComp.dat");
    for (Size i = nEnergy; i < nTarget; i++)
    {
        file << " " << str("%24.16E", target[i]) << "   " << str("%24.16E", predicted[i]) << "   "
             << str("%24.16E", predForce[i]) << "    " << std::endl;
    }
    file.close();

#else
    std::cout << "Function " << flf(VASPML_FLF) << " can only be used in debug mode." << std::endl;
#endif

    return;
}

void Refit::testFittingDesignMatrix()
{
#ifdef VASPML_USE_FIT_DEBUG
    //Vec1Real target = dataset::prepareEnergyRefitIsolatedAtom( *mlab, incar->cget<Vec1Real>(
    //"ML_EATOM_REF" ) );
    Vec1Real target = dataset::extractEnergies(*mlab);
    Size     nEnergy = target.size();
    Vec1Real forces = vector_tools::Vec2ToVec1(dataset::extractForces(*mlab));
    target.insert(target.end(), forces.cbegin(), forces.cend());
    Size nTarget = target.size();

    Vec1Real predicted(nTarget, (Real)0.0);

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    const Real* designMatrix = fitMatrix.get_fitMatrixPtr();

    linalg::matVec(linalg::Transpose::NOTRANS,
                   predicted.size(),
                   fitMatrix.get_fitMatrixSize1(),
                   1.0,
                   designMatrix,
                   fitMatrix.get_fitMatrixSize1(),
                   regressionCoefficients,
                   1,
                   1.0,
                   predicted,
                   1);

    for (Size i = 0; i < nEnergy; i++)
    {
        std::cout << "Energy " << i << "  " << target[i] << "   " << predicted[i] << std::endl;
    }

    for (Size i = nEnergy; i < nTarget; i++)
    {
        std::cout << "Force " << i - nEnergy << "   " << target[i] << "   " << predicted[i] << "   "
                  << std::abs(target[i] - predicted[i]) << "   " << target[i] - predicted[i]
                  << std::endl;
    }

    std::cout << "Compute RMSE " << std::endl;
    std::cout << math::stats::calculateRMSE(target, predicted) << std::endl;
#else
    std::cout << "Function can only be used in debug mode. Design Matrix is destroyed."
              << std::endl;
#endif

    return;
}

void Refit::SVDFittingTest()

{
    std::vector<std::vector<double>> A = {
        {1, 2},
        {3, 4},
        {5, 6}
    };
    // Example vector b (3x1)
    std::vector<double> b = {7, 8, 9};
    int                 m = A.size();    // Number of rows in A
    int                 n = A[0].size(); // Number of columns in A
    std::vector<double> A_flat(m * n);
    for (int i = 0; i < m; ++i)
    {

        for (int j = 0; j < n; ++j)
        {

            A_flat[i * n + j] = A[i][j]; // Row-major order
        }
    }
    std::vector<double> x(n, 0.0);
    std::cout << "x size " << x.size() << std::endl;
    std::cout << "matrix size " << A_flat.size() << std::endl;
    std::cout << "b size " << b.size() << std::endl;
    std::cout << "m size " << m << std::endl;
    std::cout << "n size " << n << std::endl;
    linalg::solveLeastSquaresSVD(A_flat, m, n, b, x);

    for (double xi : x) { std::cout << xi << " "; }
    std::cout << std::endl;

    return;
}

void Refit::dgemmTests()
{
    Int n = 3;
    Int m = 4;
    std::cout << "Testing transpose computation and stuff like that " << std::endl;
    Vec1Real A = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    Vec1Real ATA(m * m);
    Vec1Real b = {1.0, 2.0, 3.0};
    Vec1Real y(m);

    for (Int i = 0; i < n; i++)
    {
        for (Int j = 0; j < m; j++) { std::cout << A[i * m + j] << " "; }
        std::cout << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~" << std::endl;
    for (Int j = 0; j < m; j++)
    {
        for (Int i = 0; i < n; i++) { std::cout << A[i * m + j] << " "; }
        std::cout << std::endl;
    }

    std::cout << "Matrix times vector test " << std::endl;
    linalg::matVec(linalg::Transpose::TRANS, n, m, 1.0, A, m, b, 1, 0.0, y, 1);

    std::cout << "###################" << std::endl;
    for (const auto& yi : y) { std::cout << yi << std::endl; }

    linalg::compute_ATransTimesVector(n, m, 1.0, A, b, 0.0, y);

    std::cout << "###################" << std::endl;
    for (const auto& yi : y) { std::cout << yi << std::endl; }

    std::cout << "My matrix computation " << std::endl;
    //linalg::matMul( linalg::Transpose:TRANS, linalg::Transpose:NOTRANS,
    //                m, m, n,
    //                1.0,
    //                A, m,
    //                A, m,
    //                0.0,
    //                ATA, m );
    //for ( Int i = 0; i < m; i++ )
    //{
    //    for ( Int j = 0; j < m; j++ )
    //    {
    //        std::cout << ATA[i*m+j] << " ";
    //    }
    //    std::cout << std::endl;
    //}
    std::cout << "My matrix computation " << std::endl;
    linalg::computeATransA(n, m, 1.0, A, 0.0, ATA);
    for (Int i = 0; i < m; i++)
    {
        for (Int j = 0; j < m; j++) { std::cout << ATA[i * m + j] << " "; }
        std::cout << std::endl;
    }

    throw;
}
