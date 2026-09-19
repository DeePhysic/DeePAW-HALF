#include "FitMatrix.hpp"

#include "Linalg.hpp"
#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "KernelPolynomial.hpp"
#include "Record.hpp"
#include "ShmemArray.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm> // for copy, fill, max_element, for_each
#include <fstream>   // for basic_ostream, basic_ofstream
#include <iostream>  // for cout
#include <numeric>   // for partial_sum

using namespace vaspml;

Int vaspml::computeDesignMatrixSize(const Vec1Int& nAtoms)
{
    Int result = 0;
    for (const auto& n : nAtoms) result += 3 * n + 1;

    return result;
}

FitMatrix::FitMatrix(ShRec record) :
    data(record == nullptr ? std::make_shared<Record>() : record),
    nStrucs(0),
    nRefConfs(0),
    fitMatrix(nullptr),
    mpi(nullptr),
    nAtomsPerStructure(Vec1Int(0)),
    forceOffsets(0)
{}

FitMatrix::FitMatrix(const Vec1Int&         nAtoms,
                     const Int              nRefConfs,
                     std::shared_ptr<MlMPI> mpiIn,
                     ShRec                  record) :
    data(record == nullptr ? std::make_shared<Record>() : record),
    nStrucs(nAtoms.size()),
    nRefConfs(nRefConfs),
    fitMatrix(nullptr),
    mpi(mpiIn),
    nAtomsPerStructure(nAtoms),
    forceOffsets(nStrucs)
{}

void FitMatrix::initDesignMatrixMode()
{
    fitMatrix = std::make_shared<ShmemArray2D<Real>>(computeDesignMatrixSize(nAtomsPerStructure),
                                                     nRefConfs,
                                                     mpi);

    std::cout << "FitMatrix Size " << computeDesignMatrixSize(nAtomsPerStructure) << "   "
              << nRefConfs << std::endl;
    forceOffsets[0] = nStrucs;
    std::partial_sum(nAtomsPerStructure.cbegin(),
                     nAtomsPerStructure.cend() - 1,
                     forceOffsets.begin() + 1);
    std::transform(forceOffsets.begin() + 1,
                   forceOffsets.end(),
                   forceOffsets.begin() + 1,
                   [&](const Int& n) { return 3 * n + forceOffsets[0]; });
}

void FitMatrix::initNormalEqMode()
{
    forceOffsets[0] = nStrucs;
    fitMatrix = std::make_shared<ShmemArray2D<Real>>(nRefConfs, nRefConfs, mpi);
    Int maxAtoms = *std::max_element(nAtomsPerStructure.cbegin(), nAtomsPerStructure.cend());
    tempDesignMatrix.resize(computeDesignMatrixSize(Vec1Int(1, maxAtoms)) * nRefConfs, 0.0);
    fitMatrix = std::make_shared<ShmemArray2D<Real>>(nRefConfs, nRefConfs, mpi);
    data->put<Vec1Real>("normal-equation-target", Vec1Real(nRefConfs));
}

void FitMatrix::addStructure(const KernelPolynomial& kernelPoly, const Int strucIdx)
{
    addStructureEnergy(kernelPoly, strucIdx, fitMatrix->get_dataPointer());
    allocate_bracketTermDerivative(kernelPoly);
    addStructureForce(kernelPoly, fitMatrix->get_dataPointer(), forceOffsets[strucIdx]);
}

void FitMatrix::addStructureDesignMatrix(const KernelPolynomial& kernelPoly, const Int strucIdx)
{
    addStructureEnergy(kernelPoly, strucIdx, fitMatrix->get_dataPointer());
    allocate_bracketTermDerivative(kernelPoly);
    addStructureForce(kernelPoly, fitMatrix->get_dataPointer(), forceOffsets[strucIdx]);
}

void FitMatrix::addStructureNormalEquation(const KernelPolynomial& kernelPoly,
                                           const Real&             energyTarget,
                                           const Vec1Real&         forceTarget,
                                           const Vec1Real&         stressTarget)
{
    addStructureEnergy(kernelPoly, 0, tempDesignMatrix.data());
    addStructureForce(kernelPoly, tempDesignMatrix.data(), 1);
    const Vec1Int& nAtomsType = kernelPoly.get_nAtomsType();
    Int            nSize0 = vector_tools::sum(nAtomsType) * 3 + 7;
    linalg::computeATransA(nSize0, nRefConfs, 1.0, tempDesignMatrix, 1.0, *fitMatrix);
    // treat target values
    Vec1Real tempTarget(nSize0);
    // copy targets into contiguous array
    tempTarget[0] = energyTarget;
    std::copy(forceTarget.cbegin(), forceTarget.cend(), tempTarget.begin() + 1);
    Int offset = 1 + forceTarget.size();
    std::copy(stressTarget.cbegin(), stressTarget.cend(), tempTarget.begin() + offset);
    // now compute A^T * y and add to target
    linalg::compute_ATransTimesVector(nSize0,
                                      nRefConfs,
                                      1.0,
                                      tempDesignMatrix,
                                      tempTarget,
                                      1.0,
                                      data->get<Vec1Real>("normal-equation-target"));
}

void FitMatrix::addStructureEnergy(const KernelPolynomial& kernelPoly,
                                   const Int               strucIdx,
                                   Real*                   fitMatrix)
{
    const Vec1Int& nAtomsType = kernelPoly.get_nAtomsType();
    Real           nAtoms = (Real)vector_tools::sum(nAtomsType);

    const Vec2Real& kernelMatrix = kernelPoly.get_kernelMatrix();
    const Vec1Int&  nRefConfsType = *kernelPoly.get_numberLocalRefConfs();
    Size counter = 0; // running over all reference configurations (ref confs counted per type)
    for (Size type = 0; type < nAtomsType.size(); type++)
    {
        const Vec1Real& kernelMatrixType = kernelMatrix[type];
        for (Int refConf = 0; refConf < nRefConfsType[type]; refConf++)
        {
            for (Int atom = 0; atom < nAtomsType[type]; atom++)
            {
                fitMatrix[strucIdx * nRefConfs + counter] +=
                    kernelMatrixType[atom * nRefConfsType[type] + refConf] / nAtoms;
            }
            counter++;
        }
    }
}

void FitMatrix::addStructureForce(const KernelPolynomial& kernelPoly,
                                  Real* /* fitMatrix */,
                                  const Int /* forceOffset */)
{
    // extract kernel derivative
    //const Vec2Real& polyDerivative = kernelPoly.get_polynomialDerivative();
    // poly derivative is [nTypes][nAtoms*nRefConfs]
    // now I need the descriptor vectors normalized
    //const std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>>&
    //const auto& refConfs = kernelPoly.get_descriptorsRefConfs();
    const Vec1Int& nRefConfsType = *kernelPoly.get_numberLocalRefConfs();

    Vec1Int nRefConfOffset(nRefConfsType.size(), 0);
    std::partial_sum(nRefConfsType.cbegin(), nRefConfsType.cend() - 1, nRefConfOffset.begin() + 1);

    // [key][cntralAtom][basisFunctions]
    const std::map<String, std::shared_ptr<Descriptor>> descriptors =
        kernelPoly.get_descriptorCollection().get_descriptors();

    const std::shared_ptr<DescriptorCollector> descriptorCollector =
        kernelPoly.get_descriptorCollectionPtr();

    //const Vec1Real& norm = kernelPoly.get_descriptorCollection().get_normAtom();
    const NearestNeighborNSquare& maxNeighborList =
        *kernelPoly.get_descriptorCollection().getMaxNeighborList();

    const Vec1Int& nAtomsType = kernelPoly.get_nAtomsType();
    // loop over atom for which to compute force
    Size kAtom = 0;
    //Size refConfOffset = 0;
    for (Size kType = 0; kType < nAtomsType.size(); kType++)
    {
        for (Int kAtomInType = 0; kAtomInType < nAtomsType[kType]; kAtomInType++)
        {
            for (const String& descKey : constants::descriptorKeyList)
            {
                if (descriptors.at(descKey)->get_weight() <= 0) continue;
                // compute outer product
                //Int descSize =  descriptors.at(descKey)->get_sizeDescriptor( kAtom );
                descriptors.at(descKey)->compute_RefitDescDerivatives(kAtom, maxNeighborList);
                //descriptors.at( descKey )->computeTotalBracketTerm(
                //bracketTermDerivative[descKey], norm, kAtom );
            }
            descriptorCollector->computeBracketTerms(kAtom, descriptors, maxNeighborList);
            kAtom++;
            //std::cout << "###########################3" << std::endl;
        }
    }
    //throw;
}

void FitMatrix::computeOuterProductTerm(const String& key, const Real* XhatAtom, const Int nBasis)
{
    allocate_outerProduct(key, nBasis);

    std::fill(outerProduct[key].begin(), outerProduct[key].end(), 0.0);

    linalg::outerProduct(nBasis, nBasis, -1.0, XhatAtom, XhatAtom, outerProduct[key]);
    for (Int i = 0; i < nBasis; i++)
    {
        outerProduct[key][i * (nBasis + 1)] = 1.0 + outerProduct[key][i * (nBasis + 1)];
    }
}

void FitMatrix::write_fitMatrix(const String& filename) const
{
    fitMatrix->writeToFile(filename);
}

void FitMatrix::allocate_outerProduct(const String& key, const Size size)
{
    outerProduct[key].resize(size * size);
}

void FitMatrix::allocate_vectorTerm(const String& key, const Size size)
{
    vectorTerm[key].resize(size * 3);
}

void FitMatrix::allocate_bracketTermDerivative(const KernelPolynomial& kernelPoly)
{
    const std::map<String, std::shared_ptr<Descriptor>> descriptors =
        kernelPoly.get_descriptorCollection().get_descriptors();

    for (const String& key : constants::descriptorKeyList)
    {
        if (descriptors.at(key)->get_weight() <= 0) continue;
        const Vec1Int& numberNeighbors = descriptors.at(key)->get_neighborList().get_size();
        const Size     basisSetSize = descriptors.at(key)->get_sizeDescriptor(0);
        //const Int nAtoms = descriptors.at( key ) -> get_nAtoms();

        // TODO change maxNeighbor to maxNeighbor per type
        Int maxNeighbor = *std::max_element(numberNeighbors.cbegin(), numberNeighbors.cend());

        const Size nTypes = descriptors.at(key)->get_nAtomsType().size();
        //bracketTermDerivative[ key ].resize( nAtoms * maxNeighbor * basisSetSize * 3 );
        bracketTermDerivative[key].resize(nTypes);
        std::for_each(bracketTermDerivative[key].begin(),
                      bracketTermDerivative[key].end(),
                      [&](Vec1Real& slice)
                      {
                          slice.resize(maxNeighbor * basisSetSize * 3);
                          std::fill(slice.begin(), slice.end(), 0);
                      });
    }
}

void FitMatrix::write_fitMatrix2D(const String& filename) const
{
    auto& fitMatrix = *(this->fitMatrix);
    auto  file = file_io::openFileO(filename);
    for (Size i = 0; i < fitMatrix.get_size0(); i++)
    {
        for (Size j = 0; j < fitMatrix.get_size1(); j++)
        {
            file << str("%20.16E", fitMatrix.get_value(i, j)) << "  ";
        }
        file << std::endl;
    }
    file.close();
}

Real*& FitMatrix::get_fitMatrixPtr()
{
    return fitMatrix->get_dataPointer();
}

const Real* FitMatrix::get_fitMatrixPtr() const
{
    return fitMatrix->get_dataPointer();
}

const Vec1Real& FitMatrix::get_normalEquationTarget() const
{
    return data->cget<Vec1Real>("normal-equation-target");
}

Size FitMatrix::get_fitMatrixSize0() const
{
    return fitMatrix->get_size0();
}

Size FitMatrix::get_fitMatrixSize1() const
{
    return fitMatrix->get_size1();
}
