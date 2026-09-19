#include "TypeMap.hpp"

#include "ParallelEnvironment.hpp"
#include "debug.hpp"
#include "utils.hpp"

#include <iostream>
#include <stdexcept>

using namespace vaspml;

TypeMap::TypeMap(const Vec1String& types, const Vec1String& subTypes) :
    types(types),
    subTypes(subTypes)
{
    createMapping();
}

void TypeMap::update(const Vec1String& types, const Vec1String& subTypes)
{
    this->types = types;
    this->subTypes = subTypes;
    createMapping();

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
Int TypeMap::toType(const Int subType) const
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (subType >= (Int)mapToType.size() || subType < 0)
        {
            throw std::runtime_error("ERROR: TypeMap::toType: supplied index out of range");
        }
    );
#endif

    return mapToType[subType];
}

const Vec1Int& TypeMap::get_toType() const
{
    return mapToType;
}

Int TypeMap::toType(const String name) const
{
    auto const& entry = mapNameToType.find(name);
    if (entry == mapNameToType.end()) return -1;
    else return entry->second;
}

VASPML_EXEC_SPACE_SPECIFIER
Int TypeMap::toSubType(const Int type) const
{
#ifndef VASPML_PALGO_GPU
    VASPML_DEBUG_L1(
        if (type >= (Int)mapToSubType.size() or type < 0)
        {
            throw std::runtime_error("ERROR: TypeMap::toSubType: supplied index out of range");
        }
    );
#endif

    return mapToSubType[type];
}

Int TypeMap::toSubType(const String name) const
{
    auto const& entry = mapNameToSubType.find(name);
    if (entry == mapNameToSubType.end()) return -1;
    else return entry->second;
}

String TypeMap::nameType(const Int indx) const
{
    VASPML_DEBUG_L1(
        if (indx >= (Int)types.size() or indx < 0)
        {
            throw std::runtime_error("ERROR: TypeMap::get_typeName: supplied index out of range");
        }
    );

    return types[indx];
}

String TypeMap::nameSubType(const Int indx) const
{
    VASPML_DEBUG_L1(
        if (indx >= (Int)subTypes.size() or indx < 0)
        {
            throw std::runtime_error("ERROR: TypeMap::get_typeName: supplied index out of range");
        }
    );

    return subTypes[indx];
}

Size TypeMap::countForceFieldTypes() const
{
    return types.size();
}

Size TypeMap::countStructureTypes() const
{
    return subTypes.size();
}

void TypeMap::createMapping()
{
    mapToSubType.clear();
    mapToType.clear();
    mapToSubType.resize(types.size(), -1);
    mapToType.resize(subTypes.size(), -1);
    for (Size i = 0; i < subTypes.size(); i++)
    {
        bool found = false;
        for (Size j = 0; j < types.size(); j++)
        {
            if (string_tools::trim(subTypes[i]) == string_tools::trim(types[j]))
            {
                mapToType[i] = j;
                mapToSubType[j] = i;
                found = true;
                break;
            }
        }
        if (!found)
        {
            throw std::runtime_error("ERROR: TypeMap::make_typeMap: Could not find element "
                                     + subTypes[i] + ".");
        }
    }
    mapNameToType.clear();
    for (Size i = 0; i < types.size(); i++) { mapNameToType[types[i]] = i; }
    mapNameToSubType.clear();
    for (Size i = 0; i < subTypes.size(); i++) { mapNameToSubType[subTypes[i]] = i; }

    return;
}

void TypeMap::printTypeMap() const
{
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Writing TypeMap." << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "List of subTypes: " << std::endl;
    for (const auto& x : subTypes) std::cout << x << std::endl;
    std::cout << "List of types: " << std::endl;
    for (const auto& x : types) std::cout << x << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl;

    return;
}

const Vec1String& TypeMap::get_types() const
{
    return types;
}

const Vec1String& TypeMap::get_subTypes() const
{
    return subTypes;
}
