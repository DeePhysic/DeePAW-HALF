#include "DescriptorSHS2.hpp"

#include "BasisFunctions.hpp"
#include "BatchMap.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Timer.hpp"
#include "ThermodynamicIntegration.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm>                     // for transform, for_each, fill
#include <functional>                    // for multiplies, minus, bind, plus
#include <numeric>                       // for iota
#include <stdexcept>                     // for runtime_error

using namespace vaspml;

MapString vaspml::dataMapDescriptorSHS2()
{
    MapString m = MapString{
        {"clnmPair",                  "Vec2Real"},
        {"clnmPairDerivative",        "Vec2Real"},
        {"clnmVasp",                  "Vec2Real"},
        {"clnmDerivativeCentralVasp", "Vec2Real"},
        {"clnmPairDerivativeNO",      "Vec2Real"},
        {"clnmDerivativeRefitVasp",   "Vec2Real"}
    };

    for (const auto& dataMapEntry : dataMapDescriptor()) m.insert(dataMapEntry);

    return m;
}

DescriptorSHS2::DescriptorSHS2(const Real           weight,
                               const bool           isNormalized,
                               const DescriptorType descriptorType,
                               ShRec                record) :
    Descriptor(weight,
               isNormalized,
               descriptorType,
               assignOrMakeRecord(record, dataMapDescriptorSHS2())),
    clnmVasp(data->get<Vec2Real>("clnmVasp")),
    clnmDerivativeCentralVasp(data->get<Vec2Real>("clnmDerivativeCentralVasp")),
    clnmDerivativeRefitVasp(data->get<Vec2Real>("clnmDerivativeRefitVasp")),
    clnmPair(data->get<Vec2Real>("clnmPair")),
    clnmPairDerivative(data->get<Vec2Real>("clnmPairDerivative")),
    clnmPairDerivativeNO(data->get<Vec2Real>("clnmPairDerivativeNO"))
{
    if (record != nullptr)
    {
        if (clnmPair.size() != 0)
        {
            pairsComputed = true;
            nAtoms = clnmPair.size();
        }
        if (clnmVasp.size() != 0) { vaspFormatComputed = true; }
    }
    else
    {
        pairsComputed = false;
        vaspFormatComputed = false;
        basisFunctionsSet = false;
    }

    clnmPairOrder = ClnmPairOrder::paramMajor;
    //clnmPairOrder = ClnmPairOrder::neighborMajor;
}

Size DescriptorSHS2::compute_BasisSetSize(const Size lmax, const Vec1Size& nRoots)
{
    Size basisSetSize = 0;
    for (Size l = 0; l < lmax; l++)
    {
        for (Size m = 0; m < 2 * l + 1; m++)
        {
            for (Size root = 0; root < nRoots[l]; root++) { basisSetSize++; }
        }
    }
    return basisSetSize;
}

Size DescriptorSHS2::get_basisSetSize() const
{
    return basisSetSize;
}

void DescriptorSHS2::updatePairCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list)
{
    if (!nn_list->is_typeSorted())
    {
        global::tutor.bug("Neighbor list is not type sorted. Error in ComputeCLMs::updateCLMs");
    }
    set_neighborList(nn_list);
    nAtoms = neighborList->get_numberAtoms();
    allocatePairCoefficientArrays(basisFunctions->get_ldim());
    if (global::parallel.off())
    {
        VASPML_PROFILING_START("DescriptorSHS2::updatePairCoefficients");
        updatePairCoefficientsCPU(basisFunctions);
        VASPML_PROFILING_STOP("DescriptorSHS2::updatePairCoefficients");
    }
    else
    {
        resize_radialBasisFunction();
        resize_ylm(basisFunctions->get_ldim());
        VASPML_PROFILING_START("DescriptorSHS2::updatePairCoefficients");
        VASPML_PROFILING_START("lookUpRadialSpline.update");
        lookUpRadialSplineNN.update(nAtoms, maxOrder, nRoots, lookUpRadialSplineNNSize);
        VASPML_PROFILING_STOP("lookUpRadialSpline.update");

        VASPML_PROFILING_START("lookUpYlm.update");
        lookUpYlmNN.update(nAtoms, maxOrder, lookUpYlmNNSize);
        VASPML_PROFILING_STOP("lookUpYlm.update");

        VASPML_PROFILING_START("lookUpclnmPair.update");
        lookUpclnmPairNN.update(nAtoms, maxOrder, nRoots, lookUpclnmPairNNSize);
        VASPML_PROFILING_STOP("lookUpclnmPair.update");

        if (clnmPairOrder == ClnmPairOrder::paramMajor)
        {
            VASPML_PROFILING_START("updatePairCoefficients_stdpar");
            updatePairCoefficients_stdparNN();
            VASPML_PROFILING_STOP("updatePairCoefficients_stdpar");
        }
        else if (clnmPairOrder == ClnmPairOrder::neighborMajor)
        {
            VASPML_PROFILING_START("updatePairCoefficients_stdpar");
            updatePairCoefficients_stdparNO();
            VASPML_PROFILING_STOP("updatePairCoefficients_stdpar");
            VASPML_PROFILING_START("reorder_clnmPairDerivative_stdpar");
            reorder_clnmPairDerivative_stdpar();
            VASPML_PROFILING_STOP("reorder_clnmPairDerivative_stdpar");
        }
        VASPML_PROFILING_STOP("DescriptorSHS2::updatePairCoefficients");
    }
    pairsComputed = true;
}

void DescriptorSHS2::set_basisFunctions(const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    this->basisFunctions = basisFunctions;
    // this should be removed
    nRoots = basisFunctions->get_nValidRoots();
    // plus 1 to be able to write loops with < smbol
    maxOrder = basisFunctions->get_maxOrder() + 1;
    basisSetSize = compute_BasisSetSize(maxOrder, nRoots);
    // allocate pair coefficient arrays
    nRadialBasis = basisFunctions->get_totalNumberBesselFunctions();
    computeOffsets();
    basisFunctionsSet = true;
}

void DescriptorSHS2::updatePairCoefficientsCPU(
    const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        computePairCoefficientSingleAtom(atom, basisFunctions, *neighborList);
    }
}

void DescriptorSHS2::computeVaspCoefficientSingleAtom(const Size atom,
                                                      Vec1Real&  radialBasisFunction,
                                                      Vec1Real&  radialBasisFunctionDerivative)
{
    auto& clnmVasp = this->clnmVasp[atom];
    auto& clnmDerivativeCentralVasp = this->clnmDerivativeCentralVasp[atom];
    auto& clnmPairDerivative = this->clnmPairDerivative[atom];

    // compute angular descriptor
    const Vec1Real& distances = neighborList->get_distances(atom);
    const Vec1Real& connectionVectorNormalized = neighborList->get_connectionVectorNormalized(atom);
    const Vec1Int&  neighborTypeIndex = neighborList->get_typeIndex(atom);

    basisFunctions->computeAngularBasis(connectionVectorNormalized,
                                        distances,
                                        ylm[atom],
                                        ylmDerivative[atom]);

    Size nbasis = 0;
    Size ylm_index = 0;
    Size ylm_index_xyz = 0;
    Size ylm_index_offset = 0;
    Size ylm_index_xyz_offset = 0;

    Size neighbor = 0;
    // compute radial basis functions
    neighbor = 0;
    std::for_each(
        VASPML_POLICY_SEQ
        distances.cbegin(),
        distances.cend(),
        [&](const Real& r) mutable
        {
            // compute radial basis for atom pair
            basisFunctions->interpolate(r, radialBasisFunction, radialBasisFunctionDerivative);

            Real normed_x = connectionVectorNormalized[3 * neighbor];
            Real normed_y = connectionVectorNormalized[3 * neighbor + 1];
            Real normed_z = connectionVectorNormalized[3 * neighbor + 2];

            Size besselOffset = 0;
            Size nderivativeBasis = 0;
            Size typeOffset = neighborTypeIndex[neighbor] * basisSetSize;
            Size typeOffset3x = 3 * neighborTypeIndex[neighbor] * basisSetSize;
            for (Size ll = 0; ll < maxOrder; ll++)
            {
                for (Size nRoot = 0; nRoot < nRoots[ll]; nRoot++)
                {
                    ylm_index = ylm_index_offset;
                    ylm_index_xyz = ylm_index_xyz_offset;
                    for (Size mm = 0; mm < 2 * ll + 1; mm++)
                    {
                        // compute expansion coefficient
                        clnmVasp[typeOffset + nderivativeBasis] +=
                            radialBasisFunction[besselOffset + nRoot] * ylm[atom][ylm_index];
                        // compute derivatives
                        Real x_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_x
                                         + ylmDerivative[atom][ylm_index_xyz]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real y_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_y
                                         + ylmDerivative[atom][ylm_index_xyz + 1]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real z_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_z
                                         + ylmDerivative[atom][ylm_index_xyz + 2]
                                               * radialBasisFunction[besselOffset + nRoot];
                        clnmDerivativeCentralVasp[typeOffset3x + nderivativeBasis] -= x_component;
                        clnmDerivativeCentralVasp[typeOffset3x + nderivativeBasis] -= x_component;
                        clnmDerivativeCentralVasp[typeOffset3x + basisSetSize + nderivativeBasis] -=
                            y_component;
                        clnmDerivativeCentralVasp[typeOffset3x + 2 * basisSetSize
                                                  + nderivativeBasis] -= z_component;

                        clnmPairDerivative[3 * neighbor * basisSetSize + nderivativeBasis] =
                            x_component;
                        clnmPairDerivative[3 * neighbor * basisSetSize + basisSetSize
                                           + nderivativeBasis] = y_component;
                        clnmPairDerivative[3 * neighbor * basisSetSize + 2 * basisSetSize
                                           + nderivativeBasis] = z_component;

                        nbasis++;
                        nderivativeBasis++;
                        ylm_index++;
                        ylm_index_xyz += 3;
                    }
                }
                ylm_index_offset += 2 * ll + 1;
                ylm_index_xyz_offset += 6 * ll + 3;
                besselOffset += nRoots[ll];
            }
            neighbor++;
        });
}

void DescriptorSHS2::computePairCoefficientSingleAtom(
    const Size                             atom,
    const std::shared_ptr<BasisFunctions>& basisFunctions,
    const NearestNeighborNSquare&          NNL)
{
    // compute angular descriptor
    const Vec1Real& distances = NNL.get_distances(atom);
    const Vec1Real& connectionVectorNormalized = NNL.get_connectionVectorNormalized(atom);
    basisFunctions->computeAngularBasis(connectionVectorNormalized,
                                        distances,
                                        ylm[atom],
                                        ylmDerivative[atom]);
    Vec1Real& radialBasisFunction = (this->radialBasisFunction[atom]);
    Vec1Real& radialBasisFunctionDerivative = (this->radialBasisFunctionDerivative[atom]);

    Size nbasis = 0;
    Size ylm_index = 0;
    Size ylm_index_xyz = 0;
    Size ylm_index_offset = 0;
    Size ylm_index_xyz_offset = 0;

    Size neighbor = 0;
    // compute radial basis functions
    neighbor = 0;
    std::for_each(
        VASPML_POLICY_SEQ
        distances.cbegin(),
        distances.cend(),
        [&](const Real& r) mutable
        {
            // compute radial basis for atom pair
            basisFunctions->interpolate(r, radialBasisFunction, radialBasisFunctionDerivative);
            Real normed_x = connectionVectorNormalized[3 * neighbor];
            Real normed_y = connectionVectorNormalized[3 * neighbor + 1];
            Real normed_z = connectionVectorNormalized[3 * neighbor + 2];
            Size besselOffset = 0;
            Size nderivativeBasis = 0;
            for (Size ll = 0; ll < maxOrder; ll++)
            {
                for (Size nRoot = 0; nRoot < nRoots[ll]; nRoot++)
                {
                    ylm_index = ylm_index_offset;
                    ylm_index_xyz = ylm_index_xyz_offset;
                    for (Size mm = 0; mm < 2 * ll + 1; mm++)
                    {
                        clnmPair[atom][nbasis] =
                            radialBasisFunction[besselOffset + nRoot] * ylm[atom][ylm_index];
                        Real x_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_x
                                         + ylmDerivative[atom][ylm_index_xyz]
                                               * radialBasisFunction[besselOffset + nRoot];
                        Real y_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_y
                                         + ylmDerivative[atom][ylm_index_xyz + 1]
                                               * radialBasisFunction[besselOffset + nRoot];
                        Real z_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_z
                                         + ylmDerivative[atom][ylm_index_xyz + 2]
                                               * radialBasisFunction[besselOffset + nRoot];
                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + nderivativeBasis] =
                            x_component;
                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + basisSetSize
                                                 + nderivativeBasis] = y_component;
                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + 2 * basisSetSize
                                                 + nderivativeBasis] = z_component;
                        nbasis++;
                        nderivativeBasis++;
                        ylm_index++;
                        ylm_index_xyz += 3;
                    }
                }
                ylm_index_offset += 2 * ll + 1;
                ylm_index_xyz_offset += 6 * ll + 3;
                besselOffset += nRoots[ll];
            }
            neighbor++;
        });
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients()
{
    if (!pairsComputed)
    {
        throw std::runtime_error(
            "ERROR in routine DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients \n"
            "Pairs have to be computed first. Adapt your code such that \n"
            "void DescriptorSHS2::updatePairCoefficients( const "
            "std::shared_ptr<NearestNeighborNSquare>& nn_list,const "
            "std::unique_ptr<DescriptorBase>& descriptor )\n"
            "is called before \n");
    }
    allocateArraysVaspFormat();
    if (global::parallel.off()) { computeVaspCoefficientsFromPairCoefficientsCPU(); }
    else
    {
        lookUpClnmVasp.update(nAtoms, maxOrder, nRoots, lookUpClnmVaspSize);
        if (clnmPairOrder == ClnmPairOrder::paramMajor)
        {
            computeVaspCoefficientsFromPairCoefficientsNN_stdpar();
        }
        else if (clnmPairOrder == ClnmPairOrder::neighborMajor)
        {
            computeVaspCoefficientsFromPairCoefficientsNO_stdpar();
        }
    }
    vaspFormatComputed = true;
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficientsCPU()
{
    for (Size index_atom = 0; index_atom < neighborList->get_nAtoms(); index_atom++)
    {
        for (Size neighbor = 0; neighbor < neighborList->get_size(index_atom); neighbor++)
        {
            const auto& clnmPair_begin = clnmPair[index_atom].begin() + neighbor * basisSetSize;
            const auto& clnm_derivative_pair_begin =
                clnmPairDerivative[index_atom].begin() + 3 * neighbor * basisSetSize;
            const Size& neighborType = neighborList->get_typeIndex(index_atom, neighbor);
            Size        typeOffset = neighborType * basisSetSize;
            Size        typeOffset3x = 3 * neighborType * basisSetSize;
            std::transform( // par_unseq,
                clnmPair_begin,
                clnmPair_begin + basisSetSize,
                clnmVasp[index_atom].begin() + typeOffset,
                clnmVasp[index_atom].begin() + typeOffset,
                std::plus<Real>());
            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                clnm_derivative_pair_begin,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x,
                std::minus<Real>());

            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                clnm_derivative_pair_begin + basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                std::minus<Real>());

            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 3 * basisSetSize,
                clnm_derivative_pair_begin + 2 * basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                std::minus<Real>());
        }
    }
    vaspFormatComputed = true;
}

void DescriptorSHS2::computeVaspRefitCoefficientsFromVaspCoefficients()
{
    if (!vaspFormatComputed)
    {
        global::tutor.bug(
            "ERROR:" + flf(VASPML_FLF)
            + "clnmDerivativeCentralVasp have to be computed first. Adapt your code such that \n"
            + "void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients() is called before "
              "\n");
    }
    allocateArraysVaspRefitFormat();
    if (global::parallel.off()) computeVaspRefitCoefficientsFromVaspCoefficientsCPU();
    else
    {
        //VASPML_PARALLEL(
        //    lookUpClnmVasp.update(nAtoms, maxOrder, nRoots, lookUpClnmVaspSize);
        //    if (clnmPairOrder == ClnmPairOrder::paramMajor)
        //        computeVaspCoefficientsFromPairCoefficientsNN_stdpar();
        //    else if (clnmPairOrder == ClnmPairOrder::neighborMajor)
        //        computeVaspCoefficientsFromPairCoefficientsNO_stdpar();
        //);
        global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " GPU version not implemented yet.");
    }
    //write_clnmPair( "clnmPair.dat" );
    //write_clnmPairDerivative( "clnmPairDerivative.dat" );
    //write_clnmVasp( "clnmVasp.dat" );
    //write_clnmDerivativeCentralVasp( "clnmDerivativeCentralVasp.dat" );
}

void DescriptorSHS2::computeVaspRefitCoefficientsFromVaspCoefficientsCPU()
{
    for (Size atom = 0; atom < neighborList->get_nAtoms(); atom++)
    {
        for (Size type = 0; type < neighborList->get_numberTypes(); type++)
        {
            const Size& neighborType = type;
            Size        typeOffset3x = 3 * neighborType * basisSetSize;

            // x
            std::transform( // par_unseq,
                clnmDerivativeRefitVasp[atom].begin(),
                clnmDerivativeRefitVasp[atom].begin() + basisSetSize,
                clnmDerivativeCentralVasp[atom].begin() + typeOffset3x,
                clnmDerivativeRefitVasp[atom].begin(),
                std::minus<Real>());
            // y
            std::transform( // par_unseq,
                clnmDerivativeRefitVasp[atom].begin() + basisSetSize,
                clnmDerivativeRefitVasp[atom].begin() + 2 * basisSetSize,
                clnmDerivativeCentralVasp[atom].begin() + typeOffset3x + basisSetSize,
                clnmDerivativeRefitVasp[atom].begin() + basisSetSize,
                std::minus<Real>());
            // z
            std::transform( // par_unseq,
                clnmDerivativeRefitVasp[atom].begin() + 2 * basisSetSize,
                clnmDerivativeRefitVasp[atom].begin() + 3 * basisSetSize,
                clnmDerivativeCentralVasp[atom].begin() + typeOffset3x + 2 * basisSetSize,
                clnmDerivativeRefitVasp[atom].begin() + 2 * basisSetSize,
                std::minus<Real>());
        }
    }

    //auto file = file_io::openFileO( "DescNew.dat" );
    //for (Size atom = 0; atom < neighborList->get_nAtoms(); atom++)
    //{
    //    for ( Int i = 0; i < clnmDerivativeRefitVasp[atom].size(); i++ )
    //    {
    //        //file << str( " %24.16E ", clnmDerivativeRefitVasp[0][i] ) << std::endl;
    //        std::cout << "Desci " << str( " %24.16E ", clnmDerivativeRefitVasp[atom][i] ) <<
    //        std::endl;
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //file.close();
}

void DescriptorSHS2::computeThermodynamicIntegration(ThermodynIntegration& thermodynVars)
{
    const auto& atomList = thermodynVars.get_atomList();
    const auto& couplingConstants = thermodynVars.get_couplingConstants();

    for (Size i = 0; i < atomList.size(); i++)
    {
        Size centralAtom = atomList[i];
        std::fill(clnmVasp[centralAtom].begin(), clnmVasp[centralAtom].end(), (Real)0);
        for (Size neighbor = 0; neighbor < neighborList->get_size(centralAtom); neighbor++)
        {
            // compute basis set offset
            Size pair_offset = neighbor * basisSetSize;
            // get neighbor index and type
            const Int&  neighborIndex = neighborList->get_globalIndex(centralAtom, neighbor);
            const Size& neighborType = neighborList->get_typeIndex(centralAtom, neighbor);
            // type offset in clnmVasp
            Size typeOffset = neighborType * basisSetSize;
            for (Size nbasis = 0; nbasis < basisSetSize; nbasis++)
            {
                thermodynVars.set_clnmCoupling(centralAtom,
                                               typeOffset + nbasis,
                                               clnmVasp[centralAtom][typeOffset + nbasis]);
                clnmVasp[centralAtom][typeOffset + nbasis] +=
                    couplingConstants[neighborIndex] * clnmPair[i][pair_offset + nbasis];
            }
        }
        // rescale parts which only depend on central atom coupling constant
        std::transform(clnmPairDerivative[centralAtom].begin(),
                       clnmPairDerivative[centralAtom].end(),
                       clnmPairDerivative[centralAtom].begin(),
                       std::bind(std::multiplies<Real>(),
                                 std::placeholders::_1,
                                 couplingConstants[centralAtom]));

        std::transform(clnmDerivativeCentralVasp[centralAtom].begin(),
                       clnmDerivativeCentralVasp[centralAtom].end(),
                       clnmDerivativeCentralVasp[centralAtom].begin(),
                       std::bind(std::multiplies<Real>(),
                                 std::placeholders::_1,
                                 couplingConstants[centralAtom]));
    }
}

void DescriptorSHS2::resize_centralAtomIndex()
{
    const Size& nAtoms = neighborList->get_numberAtoms();
    bool        resize = centralAtomIndexSize.checkResize(nAtoms);
    if (resize)
    {
        centralAtomIndex.resize(nAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }
}

void DescriptorSHS2::resize_ylm(const Size& ldim)
{
    const Size&    nAtoms = neighborList->get_numberAtoms();
    const Vec1Int& numberNeighbors = neighborList->get_size();
    bool           resize = ylmSize.checkResize1Dim(nAtoms);
    if (resize)
    {
        ylmSize.act1Dim = nAtoms;
        ylmSize.max1Dim = nAtoms;
        ylmSize.resizeArray1Dim(ylm, numberNeighbors, ldim);
        ylmDerivativeSize.act1Dim = nAtoms;
        ylmDerivativeSize.max1Dim = nAtoms;
        ylmDerivativeSize.resizeArray1Dim(ylmDerivative, numberNeighbors, 3 * ldim);
    }
    else
    {
        resize = ylmSize.checkResize2Dim(numberNeighbors, ldim);
        if (resize)
        {
            ylmSize.resizeArray2Dim(ylm, numberNeighbors, ldim);
            ylmDerivativeSize.resizeArray2Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        }
    }
}

void DescriptorSHS2::allocatePairCoefficientArrays(const Size& ldim)
{
    if (global::parallel.off())
    {
        const Size& nAtoms = neighborList->get_numberAtoms();
        clnmPair.resize(nAtoms);
        clnmPairDerivative.resize(nAtoms);
        ylm.resize(nAtoms);
        ylmDerivative.resize(nAtoms);
        radialBasisFunction.resize(nAtoms);
        radialBasisFunctionDerivative.resize(nAtoms);
        for (Size atom = 0; atom < nAtoms; atom++)
        {
            const Size& numberNeighbors = neighborList->get_size(atom);
            clnmPair[atom].resize(basisSetSize * numberNeighbors);
            clnmPairDerivative[atom].resize(3 * numberNeighbors * basisSetSize);
            ylm[atom].resize(numberNeighbors * ldim);
            ylmDerivative[atom].resize(3 * numberNeighbors * ldim);
            radialBasisFunction[atom].resize(nRadialBasis);
            radialBasisFunctionDerivative[atom].resize(nRadialBasis);
        }
    }
    else resizePairArraysGPU(ldim);
}

void DescriptorSHS2::allocateArraysVaspFormat()
{
    if (global::parallel.off())
    {
        const Size& numberTypes = neighborList->get_numberTypes();
        clnmVasp.resize(neighborList->get_nAtoms());
        clnmDerivativeCentralVasp.resize(neighborList->get_nAtoms());
        for (Size i = 0; i < neighborList->get_nAtoms(); i++)
        {
            clnmVasp[i].resize(numberTypes * basisSetSize);
            std::fill(clnmVasp[i].begin(), clnmVasp[i].end(), (Real)0.0);
            clnmDerivativeCentralVasp[i].resize(3 * basisSetSize * numberTypes);
            std::fill(clnmDerivativeCentralVasp[i].begin(),
                      clnmDerivativeCentralVasp[i].end(),
                      (Real)0.0);
        }
    }
    else resizeArraysVaspFormatGPU();
}

void DescriptorSHS2::allocateArraysVaspRefitFormat()
{
    if (global::parallel.off())
    {
        clnmDerivativeRefitVasp.resize(neighborList->get_nAtoms());
        for (Size i = 0; i < neighborList->get_nAtoms(); i++)
        {
            clnmDerivativeRefitVasp[i].resize(3 * basisSetSize);
            std::fill(clnmDerivativeRefitVasp[i].begin(),
                      clnmDerivativeRefitVasp[i].end(),
                      (Real)0.0);
        }
    }
    else { global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " not implemented.\n"); }
}

void DescriptorSHS2::allocate_typeIndexArray(const Size& numberTypes)
{
    if (global::parallel.off()) resize_typeIndexArrayCPU(numberTypes);
    else resize_typeIndexArrayGPU(numberTypes);
}

void DescriptorSHS2::resize_typeIndexArrayCPU(const Size& numberTypes)
{
    typeIndexArray.resize(numberTypes);
    std::iota(typeIndexArray.begin(), typeIndexArray.end(), 0);
}

void DescriptorSHS2::computeOffsets()
{
    lOffset.resize(maxOrder);
    for (Size ll = 1; ll < maxOrder; ll++)
    {
        lOffset[ll] = lOffset[ll - 1] + (2 * (ll - 1) + 1) * nRoots[ll - 1];
    }
}

Size DescriptorSHS2::get_lOffset(const Size l) const
{
    VASPML_DEBUG_L1(
        if (l >= maxOrder)
        {
            throw std::runtime_error("ERROR in const Size DescriptorSHS2::get_lOffset"
                                     "(const Size l) const l index out of bounds");
        }
    );

    return lOffset[l];
}

Size DescriptorSHS2::get_nRadialBasis() const
{
    return nRadialBasis;
}
bool DescriptorSHS2::get_pairsComputed() const
{
    return pairsComputed;
}
bool DescriptorSHS2::get_vaspFormatComputed() const
{
    return vaspFormatComputed;
}
bool DescriptorSHS2::get_basisFunctionsSet() const
{
    return basisFunctionsSet;
}

Size DescriptorSHS2::get_maxOrder() const
{
    return maxOrder - 1;
}

const Size& DescriptorSHS2::get_nRootsOrder(const Size l) const
{
    return nRoots[l];
}

const Vec1Size& DescriptorSHS2::get_nRoots() const
{
    return nRoots;
}

const Real& DescriptorSHS2::get_clnmVasp(const Size central_atom,
                                         const Size type_number,
                                         const Size l,
                                         const Size n,
                                         const Size m) const
{
    VASPML_DEBUG_L1(
        if (!vaspFormatComputed)
        {
            throw std::runtime_error(
                "DescriptorSHS2::get_clnmVasp vasp format coefficients not computed \n"
                "please call computeVaspCoefficientsFromPairCoefficients() beforehand");
        }
        if (central_atom >= nAtoms)
        {
            throw std::runtime_error(
                "DescriptorSHS2::get_clnmVasp index central_atom out of bounds");
        }
        if (l >= maxOrder)
        {
            throw std::runtime_error("DescriptorSHS2::get_clnmVasp index l out of bounds");
        }
        if (n >= nRoots[l])
        {
            throw std::runtime_error("DescriptorSHS2::get_clnmVasp index n out of bounds");
        }
        if (m >= 2 * l + 1)
        {
            throw std::runtime_error("DescriptorSHS2::get_clnmVasp index m out of bounds");
        }
        if (!pairsComputed)
        {
            throw std::runtime_error("DescriptorSHS2::get_clnmVasp vasp pairs are not computed");
        }
    );

    return clnmVasp[central_atom][get_Index(type_number, l, n, m)];
}

VASPML_EXEC_SPACE_SPECIFIER
Size DescriptorSHS2::get_Index(const Size atom, const Size l, const Size n, const Size m) const
{
    return atom * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

VASPML_EXEC_SPACE_SPECIFIER
Size DescriptorSHS2::get_IndexX(const Size atom, const Size l, const Size n, const Size m) const
{
    return 3 * atom * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

VASPML_EXEC_SPACE_SPECIFIER
Size DescriptorSHS2::get_IndexY(const Size atom, const Size l, const Size n, const Size m) const
{
    return 3 * atom * basisSetSize + basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

VASPML_EXEC_SPACE_SPECIFIER
Size DescriptorSHS2::get_IndexZ(const Size atom, const Size l, const Size n, const Size m) const
{
    return 3 * atom * basisSetSize + 2 * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

Size DescriptorSHS2::get_IndexDir(const Size atom,
                                  const Size xyz,
                                  const Size l,
                                  const Size n,
                                  const Size m) const
{
    return 3 * atom * basisSetSize + xyz * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

//Size DescriptorSHS2::get_IndexRefit(const Size xyz, const Size atomOrType, const Size l, const
//Size n, const Size m) const
//{
//    return 3*xyz * basisSetSize + atomOrType * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
//}

const Vec1Real& DescriptorSHS2::get_clnmPair(const Size centralAtom) const
{
    return clnmPair[centralAtom];
}

const Vec2Real& DescriptorSHS2::get_clnmPair() const
{
    return clnmPair;
}

const Vec1Real& DescriptorSHS2::get_clnmPairDerivative(const Size centralAtom) const
{
    return clnmPairDerivative[centralAtom];
}

const Vec2Real& DescriptorSHS2::get_clnmPairDerivative() const
{
    return clnmPairDerivative;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& DescriptorSHS2::get_clnmVasp(const Size centralAtom) const
{
    return clnmVasp[centralAtom];
}

const Vec2Real& DescriptorSHS2::get_clnmVasp() const
{
    return clnmVasp;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& DescriptorSHS2::get_descriptor(const Size centralAtom) const
{
    return clnmVasp[centralAtom];
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec2Real& DescriptorSHS2::get_descriptor() const
{
    return clnmVasp;
}

VASPML_EXEC_SPACE_SPECIFIER
Int DescriptorSHS2::get_sizeDescriptor(const Size centralAtom) const
{
    return clnmVasp[centralAtom].size();
}

void DescriptorSHS2::rescale_descriptor(const Size centralAtom, const Real scaleFactor)
{
    std::transform( // par_unseq,
        clnmVasp[centralAtom].begin(),
        clnmVasp[centralAtom].end(),
        clnmVasp[centralAtom].begin(),
        std::bind(std::multiplies<Real>(), std::placeholders::_1, scaleFactor));
}

const Vec1Real& DescriptorSHS2::get_clnmDerivativeCentralVasp(const Size centralAtom) const
{
    return clnmDerivativeCentralVasp[centralAtom];
}

const Vec2Real& DescriptorSHS2::get_clnmDerivativeCentralVasp() const
{
    return clnmDerivativeCentralVasp;
}

const Vec1Real& DescriptorSHS2::get_clnmDerivativeRefitVasp(const Size centralAtom) const
{
    return clnmDerivativeRefitVasp[centralAtom];
}

const Vec2Real& DescriptorSHS2::get_clnmDerivativeRefitVasp() const
{
    return clnmDerivativeRefitVasp;
}

Size DescriptorSHS2::computeIndexDerivativeMatrix(Size centralAtom,
                                                  Size typeIndexForceField,
                                                  Size radialIndex,
                                                  Size totalShift)
{
    return centralAtom * totalShift + typeIndexForceField * totalShift + radialIndex;
}

void DescriptorSHS2::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                              Vec2Real&       forcePreContract,
                                              const TypeMap&  typeMap) const
{
    VASPML_PROFILING_START("DescriptorSHS2::compute_forcePreContract");
    if (global::parallel.off())
    {
        compute_forcePreContractCPU(derivativeMatrix, forcePreContract, typeMap);
    }
    else { compute_forcePreContract_stdpar(derivativeMatrix, forcePreContract); }
    VASPML_PROFILING_STOP("DescriptorSHS2::compute_forcePreContract");
}

void DescriptorSHS2::compute_forcePreContractCPU(const Vec2Real& derivativeMatrix,
                                                 Vec2Real&       forcePreContract,
                                                 const TypeMap&  typeMap) const
{
    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& get_centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();

    Size atomCounter = 0;
    Size atomShift = basisSetSize * typeMap.countForceFieldTypes();
    std::for_each(forcePreContract.begin(),
                  forcePreContract.begin() + nAtoms,
                  [&](Vec1Real& slice)
                  {
                      Size centralType = typeIndexCentral[atomCounter];
                      Size atomInType = get_centralAtomIndexPerType[atomCounter];

                      std::for_each(typeIndexArray.cbegin(),
                                    typeIndexArray.cend(),
                                    [&](Size typeInStructure)
                                    {
                                        Int neighborTypeFF = typeMap.toType(typeInStructure);
                                        for (Size nDesc = 0; nDesc < basisSetSize; nDesc++)
                                        {
                                            Size indexMatrix = atomInType * atomShift
                                                             + neighborTypeFF * basisSetSize
                                                             + nDesc;
                                            slice[get_Index(typeInStructure, 0, nDesc, 0)] =
                                                derivativeMatrix[centralType][indexMatrix];
                                        }
                                    });

                      atomCounter++;
                  });
}

Size DescriptorSHS2::get_forcePreContractSize(const Size atomIndex) const
{
    return clnmVasp[atomIndex].size();
}

const Vec1Size& DescriptorSHS2::get_lOffset() const
{
    return lOffset;
}

void DescriptorSHS2::computeForceTerms(const Vec2Real& forcePreContract,
                                       Vec2Real&       pairForces,
                                       Vec1Real&       centralForces) const
{
    VASPML_PROFILING_START("DescriptorSHS2::computeForceTerms");

    if (global::parallel.off())
    {
        computeForceTermsCPU(forcePreContract, pairForces, centralForces);
    }
    else { computeForceTerms_stdpar(forcePreContract, pairForces, centralForces); }

    VASPML_PROFILING_STOP("DescriptorSHS2::computeForceTerms");
}

void DescriptorSHS2::computeForceTermsCPU(const Vec2Real& forcePreContract,
                                          Vec2Real&       pairForces,
                                          Vec1Real&       centralForces) const
{
    const Vec2Int& neighborTypes = neighborList->get_typeIndex();
    for (Size atom = 0; atom < neighborList->get_nAtoms(); atom++)
    {
        std::for_each(typeIndexArray.cbegin(),
                      typeIndexArray.cend(),
                      [&](const Size& neighborType)
                      {
                          // x-component
                          const Real* forcePreContractPtr =
                              &forcePreContract[atom][basisSetSize * neighborType];
                          const Real* shecPtr =
                              &clnmDerivativeCentralVasp[atom][basisSetSize * neighborType * 3];
                          centralForces[3 * atom] -=
                              linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                          // y-component
                          forcePreContractPtr =
                              &forcePreContract[atom][basisSetSize * neighborType];
                          shecPtr = &clnmDerivativeCentralVasp[atom][basisSetSize * neighborType * 3
                                                                     + basisSetSize];
                          centralForces[3 * atom + 1] -=
                              linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                          // z-component
                          forcePreContractPtr =
                              &forcePreContract[atom][basisSetSize * neighborType];
                          shecPtr = &clnmDerivativeCentralVasp[atom][basisSetSize * neighborType * 3
                                                                     + 2 * basisSetSize];
                          centralForces[3 * atom + 2] -=
                              linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                      });
        Size counter = 0;
        Size neighborAtom = 0;
        std::for_each(
            neighborTypes[atom].cbegin(),
            neighborTypes[atom].cend(),
            [&](const Int& neighborType)
            {
                const Real* forcePreContractPtr =
                    &forcePreContract[atom][basisSetSize * neighborType];
                // x component
                const Real* shecPtr = &clnmPairDerivative[atom][basisSetSize * neighborAtom * 3];
                Real dotProduct = linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                // y component
                shecPtr = &clnmPairDerivative[atom][basisSetSize * neighborAtom * 3 + basisSetSize];
                dotProduct = linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                // z component
                shecPtr =
                    &clnmPairDerivative[atom][basisSetSize * neighborAtom * 3 + 2 * basisSetSize];
                dotProduct = linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                neighborAtom++;
            });
    }
}

void DescriptorSHS2::write_clnmPair(const String& fname) const
{
    file_io::writeRealVec2D(clnmPair, fname);
}

void DescriptorSHS2::write_clnmVasp(const String& fname) const
{
    file_io::writeRealVec2D(clnmVasp, fname);
}

void DescriptorSHS2::write_clnmDerivativeCentralVasp(const String& fname) const
{
    file_io::writeRealVec2D(clnmDerivativeCentralVasp, fname);
}

void DescriptorSHS2::write_clnmPairDerivative(const String& fname) const
{
    file_io::writeRealVec2D(clnmPairDerivative, fname);
}

void DescriptorSHS2::fillDescriptorRefConfSHS2(const Size&     atom,
                                               const Vec1Real& clnmPair,
                                               const Vec1Real& clnmPairDerivative,
                                               const Vec1Real& clnmVasp,
                                               const Vec1Real& clnmDerivativeCentralVasp)
{
    this->clnmPair[atom] = clnmPair;
    this->clnmPairDerivative[atom] = clnmPairDerivative;
    this->clnmVasp[atom] = clnmVasp;
    this->clnmDerivativeCentralVasp[atom] = clnmDerivativeCentralVasp;
}

void DescriptorSHS2::setParameters(const Real                             weight,
                                   const bool                             isNormalized,
                                   const DescriptorType                   descriptorType,
                                   const Size                             basisSetSize,
                                   const Vec1Size&                        nRoots,
                                   const Size                             nRadialBasis,
                                   const Size                             maxOrder,
                                   const Vec1Size&                        lOffset,
                                   const bool                             pairsComputed,
                                   const bool                             vaspFormatComputed,
                                   const bool                             basisFunctionsSet,
                                   const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    this->weight = weight;
    this->isNormalized = isNormalized;
    this->descriptorType = descriptorType;
    this->basisSetSize = basisSetSize;
    this->nRoots = nRoots;
    this->nRadialBasis = nRadialBasis;
    this->maxOrder = maxOrder;
    this->lOffset = lOffset;
    this->pairsComputed = pairsComputed;
    this->vaspFormatComputed = vaspFormatComputed;
    this->basisFunctionsSet = basisFunctionsSet;
    this->basisFunctions = basisFunctions;
}

std::shared_ptr<BasisFunctions> DescriptorSHS2::get_basisFunctions() const
{
    return basisFunctions;
}

void DescriptorSHS2::setData(const Vec2Real& clnmPair,
                             const Vec2Real& clnmPairDerivative,
                             const Vec2Real& clnmVasp,
                             const Vec2Real& clnmDerivativeCentralVasp)
{
    this->clnmPair = clnmPair;
    this->clnmPairDerivative = clnmPairDerivative;
    this->clnmVasp = clnmVasp;
    this->clnmDerivativeCentralVasp = clnmDerivativeCentralVasp;
    nAtoms = clnmPair.size();
}

void DescriptorSHS2::extendData(const Vec2Real& clnmPair,
                                const Vec2Real& clnmPairDerivative,
                                const Vec2Real& clnmVasp,
                                const Vec2Real& clnmDerivativeCentralVasp)
{
    this->clnmPair.insert(this->clnmPair.end(), clnmPair.cbegin(), clnmPair.cend());
    this->clnmPairDerivative.insert(this->clnmPairDerivative.end(),
                                    clnmPairDerivative.cbegin(),
                                    clnmPairDerivative.cend());
    this->clnmVasp.insert(this->clnmVasp.end(), clnmVasp.cbegin(), clnmVasp.cend());
    this->clnmDerivativeCentralVasp.insert(this->clnmDerivativeCentralVasp.end(),
                                           clnmDerivativeCentralVasp.cbegin(),
                                           clnmDerivativeCentralVasp.cend());
    nAtoms += clnmPair.size();
}

void DescriptorSHS2::addElement(const Vec1Real& clnmPair,
                                const Vec1Real& clnmPairDerivative,
                                const Vec1Real& clnmVasp,
                                const Vec1Real& clnmDerivativeCentralVasp)
{

    this->clnmPair.push_back(clnmPair);
    this->clnmPairDerivative.push_back(clnmPairDerivative);
    this->clnmVasp.push_back(clnmVasp);
    this->clnmDerivativeCentralVasp.push_back(clnmDerivativeCentralVasp);
}

void DescriptorSHS2::reorderSHS2DataWithBatchMap(const AtomBatchMap& batchMap)
{
    Vec2Real clnmPairCopy(nAtoms);
    Vec2Real clnmPairDerivativeCopy(nAtoms);
    Vec2Real clnmVaspCopy(nAtoms);
    Vec2Real clnmDerivativeCentralVaspCopy(nAtoms);
    for (Size i = 0; i < clnmPairCopy.size(); i++)
    {
        const Int& posOld = batchMap.get_mapTO_OrigOrder(i);
        clnmPairCopy[i] = clnmPair[posOld];
        clnmPairDerivativeCopy[i] = clnmPairDerivative[posOld];
        clnmVaspCopy[i] = clnmVasp[posOld];
        clnmDerivativeCentralVaspCopy[i] = clnmDerivativeCentralVasp[posOld];
    }
    clnmPair = std::move(clnmPairCopy);
    clnmPairDerivative = std::move(clnmPairDerivativeCopy);
    clnmVasp = std::move(clnmVaspCopy);
    clnmDerivativeCentralVasp = std::move(clnmDerivativeCentralVaspCopy);
}

Int DescriptorSHS2::get_totalNumberBasisFunctions() const
{
    VASPML_DEBUG_L1(
        if (!basisFunctionsSet)
        {
            global::tutor.bug(
                "ERROR in DescriptorSHS2::get_totalNumberBasisFunctions() basis functions not set");
        }
    );
    return basisFunctions->get_totalNumberBasisFunctions();
}

void DescriptorSHS2::compute_RefitDescDerivatives(const Int                     kAtom,
                                                  const NearestNeighborNSquare& maxNeighborList)
{
    compute_dC00(kAtom, maxNeighborList);
    compute_dC00Head(kAtom);
}

void DescriptorSHS2::compute_dC00(const Int kAtom, const NearestNeighborNSquare& maxNeighborList)
{
    allocate_dC00(maxNeighborList);

    const Vec1Int& globalIndex = maxNeighborList.get_globalIndex(kAtom);
    const Int      typeShift = basisSetSize * neighborList->get_numberTypes();

    //std::ofstream outputF( "CLNM.dat", std::ios::app );
    //for ( const auto& x : clnmVasp ) for ( const auto& y : x ) outputF << y << "  " << std::endl;

    //throw;
    // IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0))
    for (Size centralAtom = 0; centralAtom < maxNeighborList.get_size(kAtom); centralAtom++)
    {
        const Int centralGlobal = globalIndex[centralAtom];
        //const Int centralType = typeIndex[centralAtom];
        Int col = 0;
        for (Size jType = 0; jType < neighborList->get_numberTypes(); jType++)
        {
            // if higher orders of l are needed. please insert l and m loop here
            for (Size nRoot = 0; nRoot < nRoots[0]; nRoot++)
            {
                const Int indx = get_Index(jType, (Size)0, nRoot, (Size)0);
                dC00[centralAtom][col] = clnmVasp[centralGlobal][indx];
                col++;
            }
        }

        if (centralAtom < neighborList->get_size(kAtom))
        {
            col = 0;
            for (Size nRoot = 0; nRoot < nRoots[0]; nRoot++)
            {
                const Int indxX = get_IndexX(centralAtom, (Size)0, nRoot, (Size)0);
                const Int indxY = get_IndexY(centralAtom, (Size)0, nRoot, (Size)0);
                const Int indxZ = get_IndexZ(centralAtom, (Size)0, nRoot, (Size)0);

                dC00[centralAtom][typeShift + col] = clnmPairDerivative[kAtom][indxX];
                dC00[centralAtom][2 * typeShift + col] = clnmPairDerivative[kAtom][indxY];
                dC00[centralAtom][3 * typeShift + col] = clnmPairDerivative[kAtom][indxZ];
                col++;
            }
        }
    }
    //std::ofstream outputFile( "dC00.dat", std::ios::app );
    //for ( Int centralAtom = 0; centralAtom < maxNeighborList.get_size( kAtom ); centralAtom++ )
    //{
    //    //for ( Int j = 0; j <  neighborList->get_numberTypes()*basisSetSize; j++ )
    //    for ( Int j = 0; j <  dC00[centralAtom].size(); j++ )
    //    {
    //        outputFile << dC00[centralAtom][j] << std::endl;
    //    }
    //}
}

void DescriptorSHS2::compute_dC00Head(const Int kAtom)
{
    allocate_dC00Head();
    Int       kType = neighborList->get_typeIndexCentral(kAtom);
    const Int typeShift = basisSetSize * neighborList->get_numberTypes();
    //std::cout << "TYPE SHIFT " << typeShift << std::endl;
    Int col = 0;
    for (Size jType = 0; jType < neighborList->get_numberTypes(); jType++)
    {
        // if higher orders of l are needed. please insert l and m loop here
        for (Size nRoot = 0; nRoot < nRoots[0]; nRoot++)
        {
            const Int indx = get_Index(jType, (Size)0, nRoot, (Size)0);
            dC00Head[kType][col] = clnmVasp[kAtom][indx];
            col++;
        }
    }
    col = 0;
    for (Size jType = 0; jType < neighborList->get_numberTypes(); jType++)
    {
        // if higher orders of l are needed. please insert l and m loop here
        for (Size nRoot = 0; nRoot < nRoots[0]; nRoot++)
        {
            const Int indxX = get_IndexX(jType, (Size)0, nRoot, (Size)0);
            const Int indxY = get_IndexY(jType, (Size)0, nRoot, (Size)0);
            const Int indxZ = get_IndexZ(jType, (Size)0, nRoot, (Size)0);
            //const Int indx
            dC00Head[kType][typeShift + col] = clnmDerivativeCentralVasp[kAtom][indxX];
            dC00Head[kType][2 * typeShift + col] = clnmDerivativeCentralVasp[kAtom][indxY];
            dC00Head[kType][3 * typeShift + col] = clnmDerivativeCentralVasp[kAtom][indxZ];
            col++;
        }
    }

    //std::ofstream outputFile("dC00Head.dat", std::ios::app);
    //for ( Int j = 0; j <  dC00Head[kType].size(); j++ )
    ////for ( Int j = 0; j < typeShift; j++ )
    //{
    //    //std::cout << << dC00Head[kType][j]<< std::endl;
    //    outputFile << dC00Head[kType][j]<< std::endl;
    //}
    //outputFile.close();
}

void DescriptorSHS2::allocate_dC00(const NearestNeighborNSquare& maxNeighborList)
{
    if (global::parallel.off())
    {
        const Vec1Int& numberNeighbors = maxNeighborList.get_size();
        Int maxNeighbors = *std::max_element(numberNeighbors.cbegin(), numberNeighbors.cend());
        dC00.resize(maxNeighbors);
        std::for_each(dC00.begin(),
                      dC00.end(),
                      [&](Vec1Real& slice)
                      {
                          slice.resize(neighborList->get_numberTypes() * basisSetSize * 4);
                          std::fill(slice.begin(), slice.end(), 0);
                      });
    }
    else global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " GPU not implemented yet.");
}

void DescriptorSHS2::allocate_dC00Head()
{
    if (global::parallel.off())
    {
        dC00Head.resize(neighborList->get_numberTypes());
        std::for_each(dC00Head.begin(),
                      dC00Head.end(),
                      [&](Vec1Real& slice)
                      {
                          slice.resize(neighborList->get_numberTypes() * basisSetSize * 4);
                          std::fill(slice.begin(), slice.end(), 0);
                      });
    }
    else global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " GPU not implemented yet.");
}

const Vec2Real& DescriptorSHS2::get_refitHeadTerm() const
{
    return dC00Head;
}

void DescriptorSHS2::reorderElements(const Vec1Int& newOrder)
{
    if (pairsComputed)
    {
        vector_tools::inplaceReorderVector(clnmPair, newOrder);
        vector_tools::inplaceReorderVector(clnmPairDerivative, newOrder);
    }
    if (vaspFormatComputed)
    {
        vector_tools::inplaceReorderVector(clnmVasp, newOrder);
        vector_tools::inplaceReorderVector(clnmDerivativeCentralVasp, newOrder);
    }
}

namespace vaspml
{

void copyParameters(const DescriptorSHS2& descA, DescriptorSHS2& descB)
{
    descB.setParameters(descA.get_weight(),
                        descA.get_isNormalized(),
                        descA.get_descriptorType(),
                        descA.get_basisSetSize(),
                        descA.get_nRoots(),
                        descA.get_nRadialBasis(),
                        descA.get_maxOrder(),
                        descA.get_lOffset(),
                        descA.get_pairsComputed(),
                        descA.get_vaspFormatComputed(),
                        descA.get_basisFunctionsSet(),
                        descA.get_basisFunctions());
}

void copyNeighborList(const DescriptorSHS2& descA, DescriptorSHS2& descB)
{
    descB.set_neighborList(descA.get_neighborList_ptr());
}

void copyBasisFunctions(const DescriptorSHS2& descA, DescriptorSHS2& descB)
{
    descB.set_basisFunctions(descA.get_basisFunctions());
}

void copyData(const DescriptorSHS2& descA, DescriptorSHS2& descB)
{
    descB.setData(descA.get_clnmPair(),
                  descA.get_clnmPairDerivative(),
                  descA.get_clnmVasp(),
                  descA.get_clnmDerivativeCentralVasp());
}

void extendData(const DescriptorSHS2& descA, DescriptorSHS2& descB)
{
    descB.extendData(descA.get_clnmPair(),
                     descA.get_clnmPairDerivative(),
                     descA.get_clnmVasp(),
                     descA.get_clnmDerivativeCentralVasp());
}

void addElementData(const DescriptorSHS2& descA, DescriptorSHS2& descB, const Size atomIdx)
{
    descB.addElement(descA.get_clnmPair(atomIdx),
                     descA.get_clnmPairDerivative(atomIdx),
                     descA.get_clnmVasp(atomIdx),
                     descA.get_clnmDerivativeCentralVasp(atomIdx));
}

} //namespace vaspml
