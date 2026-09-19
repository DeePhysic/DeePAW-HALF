#include "Descriptor.hpp"

#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "DescriptorSHS3ReducedLinElem.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>

using namespace vaspml;

MapString vaspml::dataMapDescriptor()
{
    return MapString{};
}

Descriptor::Descriptor(ShRec record) :
    data(assignOrMakeRecord(record, dataMapDescriptor())),
    descriptorType(DescriptorType::none)
{
    this->weight = (Real)0;
    this->isNormalized = false;
    this->derivativeMapComputed = false;
    this->neighborList = nullptr;
}

Descriptor::Descriptor(const Real           weight,
                       const bool           isNormalized,
                       const DescriptorType descriptorType,
                       ShRec                record) :
    data(assignOrMakeRecord(record, dataMapDescriptor())),
    descriptorType(descriptorType)
{
    this->weight = weight;
    this->isNormalized = isNormalized;
    this->derivativeMapComputed = false;
    this->neighborList = nullptr;
}

VASPML_EXEC_SPACE_SPECIFIER
Real Descriptor::get_weight() const
{
    return weight;
}

VASPML_EXEC_SPACE_SPECIFIER
bool Descriptor::get_isNormalized() const
{
    return isNormalized;
}

DescriptorType Descriptor::get_descriptorType() const
{
    return descriptorType;
}

void Descriptor::set_neighborList(const std::shared_ptr<const NearestNeighborNSquare>& neighborList)
{
    this->neighborList = neighborList;
}

VASPML_EXEC_SPACE_SPECIFIER
const Int& Descriptor::get_typeIndexCentral(const Size atomIndex) const
{
    return neighborList->get_typeIndexCentral(atomIndex);
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Int& Descriptor::get_typeIndexCentral() const
{
    return neighborList->get_typeIndexCentral();
}

const Int& Descriptor::get_typeIndex(const Size atomIndex, const Size indxNeigh) const
{
    return neighborList->get_typeIndex(atomIndex, indxNeigh);
}

Int Descriptor::get_nAtoms() const
{
    return neighborList->get_nAtoms();
}

const Vec1Int& Descriptor::get_nAtomsType() const
{
    return neighborList->get_nAtomsType();
}

VASPML_EXEC_SPACE_SPECIFIER
const NearestNeighborNSquare& Descriptor::get_neighborList() const
{
    return *neighborList;
}

VASPML_EXEC_SPACE_SPECIFIER
const std::shared_ptr<const NearestNeighborNSquare>& Descriptor::get_neighborList_ptr() const
{
    return neighborList;
}

const ShRec& Descriptor::get_dataRecord() const
{
    return data;
}

ShRec& Descriptor::get_dataRecord()
{
    return data;
}

//
// implementation of member functions which were virtual
//
VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& Descriptor::get_descriptor(const Size centralAtom) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_descriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_descriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_descriptor(centralAtom);
    }
    else
    {
#ifndef VASPML_PALGO_GPU
        throw std::runtime_error("ERROR: Invalid descriptor type.");
#endif
    }
}

VASPML_EXEC_SPACE_SPECIFIER
Int Descriptor::get_sizeDescriptor(const Size centralAtom) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_sizeDescriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_sizeDescriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_sizeDescriptor(
            centralAtom);
    }
    else return -1;
}

void Descriptor::rescale_descriptor(const Size centralAtom, const Real scaleFactor)
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<DescriptorSHS2*>(this)->rescale_descriptor(centralAtom, scaleFactor);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<DescriptorSHS3*>(this)->rescale_descriptor(centralAtom, scaleFactor);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<DescriptorSHS3ReducedLinElem*>(this)->rescale_descriptor(centralAtom,
                                                                                    scaleFactor);
    }
    else return;
}

void Descriptor::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                          Vec2Real&       forcePreContract,
                                          const TypeMap&  typeMap) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->compute_forcePreContract(derivativeMatrix,
                                                                                  forcePreContract,
                                                                                  typeMap);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->compute_forcePreContract(derivativeMatrix,
                                                                                  forcePreContract,
                                                                                  typeMap);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->compute_forcePreContract(
            derivativeMatrix,
            forcePreContract,
            typeMap);
    }
    else return;
}

Size Descriptor::get_forcePreContractSize(const Size atomIndex) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_forcePreContractSize(atomIndex);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_forcePreContractSize(atomIndex);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_forcePreContractSize(
            atomIndex);
    }
    else return std::numeric_limits<Size>::max();
}

void Descriptor::computeForceTerms(const Vec2Real& forcePreContract,
                                   Vec2Real&       pairForces,
                                   Vec1Real&       centralForces) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->computeForceTerms(forcePreContract,
                                                                           pairForces,
                                                                           centralForces);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->computeForceTerms(forcePreContract,
                                                                           pairForces,
                                                                           centralForces);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->computeForceTerms(
            forcePreContract,
            pairForces,
            centralForces);
    }
    else return;
}

void Descriptor::resizeArraysRefConf(const Vec1Size& sizes)
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<DescriptorSHS2*>(this)->resizeArraysRefConf(sizes);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<DescriptorSHS3*>(this)->resizeArraysRefConf(sizes);
    }
}

void Descriptor::fillDescriptorRefConfSHS3(const Size& atom, const Vec1Real& spectrumVasp)
{
    if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<DescriptorSHS3*>(this)->fillDescriptorRefConfSHS3(atom, spectrumVasp);
    }
    else
    {
        throw std::runtime_error(
            "ERROR: Descriptor::fillDescriptorRefConfSHS3 only implemented for DescriptorSHS3!!\n");
    }
}

void Descriptor::fillDescriptorRefConfSHS2(const Size&     atom,
                                           const Vec1Real& clnmPair,
                                           const Vec1Real& clnmPairDerivative,
                                           const Vec1Real& clnmVasp,
                                           const Vec1Real& clnmDerivativeCentralVasp)
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<DescriptorSHS2*>(this)->fillDescriptorRefConfSHS2(
            atom,
            clnmPair,
            clnmPairDerivative,
            clnmVasp,
            clnmDerivativeCentralVasp);
    }
    else
    {
        throw std::runtime_error(
            "ERROR: Descriptor::fillDescriptorRefConfSHS2 only implemented for DescriptorSHS2!!\n");
    }
}

Descriptor Descriptor::copy(const ShRec& data) const
{
    if (descriptorType == DescriptorType::none)
    {
        Descriptor descNew(this->get_weight(),
                           this->get_isNormalized(),
                           DescriptorType::none,
                           data);
        return descNew;
    }
    else
    {
        global::tutor.bug(
            "ERROR: " + flf(VASPML_FLF)
            + " doing copy with descriptorType!=DescriptorType::none is impossible!.");
        throw;
    }
}

Int Descriptor::get_totalNumberBasisFunctions() const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_totalNumberBasisFunctions();
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_totalNumberBasisFunctions();
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)
            ->get_totalNumberBasisFunctions();
    }
    else
    {
        throw std::runtime_error("ERROR: Descriptor::get_totalNumberBasisFunctions() only "
                                 "implemented for DescriptorSHS2 and DescriptorSHS3!!\n");
    }
}

void Descriptor::compute_RefitDescDerivatives(const Int                     kAtom,
                                              const NearestNeighborNSquare& maxNeighborList)
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<DescriptorSHS2*>(this)->compute_RefitDescDerivatives(kAtom,
                                                                                maxNeighborList);
    }
    if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<DescriptorSHS3*>(this)->compute_RefitDescDerivatives(kAtom,
                                                                                maxNeighborList);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        std::cout << "NOT implemented yet not implemented yet " << std::endl;
    }
    else
    {
        global::tutor.bug("Error: " + flf(VASPML_FLF) + ". No valid descriptor type selected \n");
    }
}

const Vec2Real& Descriptor::get_refitHeadTerm() const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_refitHeadTerm();
    }
    if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_refitHeadTerm();
    }
    else
    {
        global::tutor.bug("Error: " + flf(VASPML_FLF) + ". No valid descriptor type selected \n");
        throw;
    }
}
