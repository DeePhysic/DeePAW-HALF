#include "KernelPolynomial.hpp"

#include "DescriptorCollector.hpp"
#include "Frame.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironment.hpp"
#include "TypeMap.hpp"
#include "Record.hpp"
#include "constants.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm> // for for_each
#include <cmath>     // for sqrt
#include <tuple>     // for tuple, get

using namespace vaspml;

MapString vaspml::dataMapKernelPolynomial()
{
    MapString m = {
        {"polynomialDerivative", "Vec2Real"},
        {"polynomialKernelNorm", "Vec2Real"}
    };

    for (const auto& dataMapEntry : dataMapKernel()) m.insert(dataMapEntry);

    return m;
}

KernelPolynomial::KernelPolynomial(ShRec record) :
    Kernel(assignOrMakeRecord(record, dataMapKernelPolynomial())),
    polynomialDerivative(data->get<Vec2Real>("polynomialDerivative")),
    polynomialKernelNorm(data->get<Vec2Real>("polynomialKernelNorm")),
    polynomialKernel(kernelMatrix)
{}

KernelPolynomial::KernelPolynomial(const Record&                 inputParameters,
                                   const std::shared_ptr<MlMPI>& mpiIn,
                                   ShRec                         record) :
    Kernel(inputParameters, mpiIn, assignOrMakeRecord(record, dataMapKernelPolynomial())),
    polynomialDerivative(data->get<Vec2Real>("polynomialDerivative")),
    polynomialKernelNorm(data->get<Vec2Real>("polynomialKernelNorm")),
    polynomialKernel(kernelMatrix)
{
    kernelPower = inputParameters.cget<Int>("ML_NHYP");
}

KernelPolynomial::KernelPolynomial(const Record&                 incar,
                                   const Record&                 mlab,
                                   const std::shared_ptr<MlMPI>& mpiIn,
                                   ShRec                         record) :
    Kernel(incar, mlab, mpiIn, assignOrMakeRecord(record, dataMapKernelPolynomial())),
    polynomialDerivative(data->get<Vec2Real>("polynomialDerivative")),
    polynomialKernelNorm(data->get<Vec2Real>("polynomialKernelNorm")),
    polynomialKernel(kernelMatrix)
{
    kernelPower = incar.cget<Int>("ML_NHYP");
}

void KernelPolynomial::init(const Record& inputParameters, const std::shared_ptr<MlMPI>& mpiIn)
{
    Kernel::init(inputParameters, mpiIn);
    kernelPower = inputParameters.cget<Int>("ML_NHYP");
}

void KernelPolynomial::allocateDerivativeArrays()
{
    if (global::parallel.off())
    {
        polynomialDerivative.resize(this->kernelMatrix.size());
        polynomialKernelNorm.resize(this->kernelMatrix.size());
        for (Size type = 0; type < polynomialDerivative.size(); type++)
        {
            polynomialDerivative[type].resize(kernelMatrix[type].size());
            polynomialKernelNorm[type].resize(kernelMatrix[type].size());
        }
        polynomialKernel = kernelMatrix;
    }
    else { allocateDerivativeArraysGPU(); }
}

void KernelPolynomial::allocateDerivativeArraysRefit()
{

    if (global::parallel.off())
    {
        polynomialDerivative.resize(this->kernelMatrix.size());
        for (Size type = 0; type < polynomialDerivative.size(); type++)
        {
            polynomialDerivative[type].resize(kernelMatrix[type].size());
        }
        polynomialKernel = kernelMatrix;
    }
    else { allocateDerivativeArraysRefitGPU(); }
}

void KernelPolynomial::updatePolynomialKernel(const Frame& frame)
{
    updateKernel(frame);
    allocateDerivativeArrays();
    //nvtxRangePush( "computePolynomialKernelTerms" );
    computePolynomialKernelTerms(*frame.get_typeMap());
    //nvtxRangePop();
}

void KernelPolynomial::updatePolynomialKernelRefit(const Frame& frame)
{
    updateKernel(frame);
    allocateDerivativeArraysRefit();
    //nvtxRangePush( "computePolynomialKernelTerms" );
    computePolynomialKernelTermsRefit(*frame.get_typeMap());
    //nvtxRangePop();
}

void KernelPolynomial::computePolynomialKernelTerms(const TypeMap& typeMap)
{
    // actual calculations
    if (global::parallel.off())
    {
        // transfer data in a std algorithm form
        std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data(polynomialDerivative.size());
        for (Size i = 0; i < polynomialKernel.size(); i++)
        {
            for (Size j = 0; j < polynomialKernel[i].size(); j++)
            {
                data[i].push_back({polynomialKernel[i][j],
                                   polynomialDerivative[i][j],
                                   polynomialKernelNorm[i][j]});
            }
        }
        computePolynomialKernelTermsCPU(typeMap, data);
    }
    else
    {
        makeLookUpTables(typeMap);
        computePolynomialKernelTerms_stdlib();
    }
}

void KernelPolynomial::computePolynomialKernelTermsRefit(const TypeMap& typeMap)
{
    // actual calculations
    if (global::parallel.off())
    {
        // transfer data in a std algorithm form
        std::vector<std::vector<std::tuple<Real&, Real&>>> data(polynomialDerivative.size());
        for (Size i = 0; i < polynomialKernel.size(); i++)
        {
            for (Size j = 0; j < polynomialKernel[i].size(); j++)
            {
                data[i].push_back({polynomialKernel[i][j], polynomialDerivative[i][j]});
            }
        }
        computePolynomialKernelTermsRefitCPU(typeMap, data);
    }
    else
    {
        makeLookUpTables(typeMap);
        computePolynomialKernelTermsRefit_stdlib();
    }
}

void KernelPolynomial::computePolynomialKernelTermsCPU(
    const TypeMap&                                            typeMap,
    std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data)
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();
    const Vec1Int&  numberLocalRefConfs = (*this->numberLocalRefConfs);
    const Int&      kernelPower = (this->kernelPower);
    Size            atomCounter = 0;
    for (Size type = 0; type < typeMap.countStructureTypes(); type++)
    {
        Size typeFF = typeMap.toType(type);
        Size numberLocalRefConfsType = numberLocalRefConfs[typeFF];
        Size increaseCounter = 0;
        for (Size centralAtom = 0; centralAtom < (Size)nAtomsPerType[type]; centralAtom++)
        {
            Real norm = 0;
            if (normAtoms[atomCounter] > constants::EPS_TOL)
            {
                norm = (Real)kernelPower / normAtoms[atomCounter];
            }
            std::for_each( // par,
                data[type].begin() + increaseCounter,
                data[type].begin() + increaseCounter + (Size)numberLocalRefConfsType,
                [&norm, &kernelPower](std::tuple<Real&, Real&, Real&>& kernelValue)
                {
                    Real factor = math::intPowMixAlgo(std::get<0>(kernelValue), kernelPower - 1);
                    std::get<1>(kernelValue) = factor * norm;
                    std::get<0>(kernelValue) = factor * std::get<0>(kernelValue);
                    std::get<2>(kernelValue) = std::get<0>(kernelValue) * norm;
                });
            atomCounter++;
            increaseCounter += numberLocalRefConfsType;
        }
    }
}

void KernelPolynomial::computePolynomialKernelTermsRefitCPU(
    const TypeMap&                                     typeMap,
    std::vector<std::vector<std::tuple<Real&, Real&>>> data)
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();
    const Vec1Int&  numberLocalRefConfs = (*this->numberLocalRefConfs);
    const Int&      kernelPower = (this->kernelPower);
    Size            atomCounter = 0;
    for (Size type = 0; type < typeMap.countStructureTypes(); type++)
    {
        Size typeFF = typeMap.toType(type);
        Size numberLocalRefConfsType = numberLocalRefConfs[typeFF];
        Size increaseCounter = 0;
        for (Size centralAtom = 0; centralAtom < (Size)nAtomsPerType[type]; centralAtom++)
        {
            Real factor1 = (Real)kernelPower;
            if (normAtoms[atomCounter] > constants::EPS_TOL)
            {
                factor1 = factor1 / normAtoms[atomCounter];
            }
            std::for_each( // par,
                data[type].begin() + increaseCounter,
                data[type].begin() + increaseCounter + (Size)numberLocalRefConfsType,
                [&factor1, &kernelPower](std::tuple<Real&, Real&>& kernelValue)
                {
                    Real factor = math::intPowMixAlgo(std::get<0>(kernelValue), kernelPower - 1);
                    std::get<1>(kernelValue) = factor * factor1;
                    std::get<0>(kernelValue) = factor * std::get<0>(kernelValue);
                });
            atomCounter++;
            increaseCounter += numberLocalRefConfsType;
        }
    }
}

void KernelPolynomial::computeDerivativeKernel(Vec1Real&     derivativeMatrix,
                                               Vec1Real&     forceVector,
                                               const Real*   descriptorsRefConfsWeighted,
                                               const Real*   regressionCoefficients,
                                               const Size&   typeStruc,
                                               const Size&   typeForceField,
                                               const String& key) const
{
    //cudaDeviceSynchronize();
    linalg::matMul(linalg::Transpose::NOTRANS,
                   linalg::Transpose::NOTRANS,
                   nAtomsPerType[typeStruc],
                   (*numberDescriptors.at(key))[typeForceField],
                   (*numberLocalRefConfs)[typeForceField],
                   std::sqrt(weights.at(key)),
                   polynomialDerivative[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   descriptorsRefConfsWeighted,
                   (*numberDescriptors.at(key))[typeForceField],
                   (Real)0.0,
                   derivativeMatrix,
                   (*numberDescriptors.at(key))[typeForceField]);

    linalg::matVec(linalg::Transpose::NOTRANS,
                   nAtomsPerType[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   -std::sqrt(weights.at(key)),
                   polynomialKernelNorm[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   regressionCoefficients,
                   1,
                   0.0,
                   forceVector,
                   1);

    const std::map<String, ShVec2Real>& descriptors =
        descriptorCollection->get_descriptorsNormalized();
    // this is needed. the reason is not sure but I think:
    // We have the host codes with the pointers. But there is no autosync as Matthias Noack
    // mentioned in his email:
    // "First, cuBLAS (computation) calls are non-blocking, that means they return immediately after
    // issuing the work for the GPU. cuBLASgetVector() is blocking though, which means without it,
    // you might access and print the results before the computations on the GPU are finished. In
    // general you can use cudaDeviceSynchronize() to wait for the GPU computations on the host, or
    // more fine grained per handle/stream via something like:"
#ifdef VASPML_PALGO_GPU_NVIDIA
    if (global::parallel.gpu()) cudaDeviceSynchronize();
#endif
    for (Size atom = 0; atom < (Size)nAtomsPerType[typeStruc]; atom++)
    {
        Real* derivativeMatrixPtr =
            &derivativeMatrix[atom * (*numberDescriptors.at(key))[typeForceField]];
        const Real* descriptorPtr =
            &(*descriptors.at(key))[typeStruc][atom * (*numberDescriptors.at(key))[typeForceField]];
        linalg::scaleVectorPlusVector(forceVector[atom],
                                      descriptorPtr,
                                      derivativeMatrixPtr,
                                      (*numberDescriptors.at(key))[typeForceField]);
    }
}

void KernelPolynomial::write_polynomialDerivative(const String& fname) const
{
    file_io::writeRealVec2D(polynomialDerivative, fname);
}
void KernelPolynomial::write_polynomialKernelNorm(const String& fname) const
{
    file_io::writeRealVec2D(polynomialKernelNorm, fname);
}
void KernelPolynomial::write_polynomialKernel(const String& fname) const
{
    file_io::writeRealVec2D(polynomialKernel, fname);
}

const Vec2Real& KernelPolynomial::get_polynomialDerivative() const
{
    return polynomialDerivative;
}

const Vec1Real& KernelPolynomial::get_polynomialDerivative(const Int& typeIndx) const
{
    return polynomialDerivative[typeIndx];
}
