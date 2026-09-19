#ifndef TYPEMAP_HPP
#define TYPEMAP_HPP

#include "types.hpp"

namespace vaspml
{

/***************************************************************************************************
 * Provides an integer mapping between atom types of structure and force field.
 *
 * This class is used when the number of atom types in the already generated machine
 * learning force field and the actual structure do not match. The structure types
 * have to be a subset of the machine learning force field types.
 *
 * Example:
 *
 * * Types in force field file:  `H O Na Cl S`
 * * Types in structure:  `Cl O S`
 *
 * Then, the mapping looks like this:
 *
 * * `mapToType = { 3 1 4 }`
 * * `mapToSubType = { -1  1 -1  0  2 }`
 *
 * @note If there is no mapping between a type and subtype (because the force field defines more
 * types than the structure) then the toSubType() function returns -1.
 **************************************************************************************************/
class TypeMap
{
  public:
    /***********************************************************************************************
     * Construct empty TypeMap.
     *
     * @warning The member function update() must be called before use.
     **********************************************************************************************/
    TypeMap() = default;
    /***********************************************************************************************
     * Construct TypeMap from input types.
     *
     * @param types The types in the machine learning force field.
     * @param subTypes The types in the actual structure.
     *
     * The condition `subTypes` @f$\subseteq@f$ `types` must be statisfied.
     **********************************************************************************************/
    TypeMap(const Vec1String& types, const Vec1String& subTypes);
    /***********************************************************************************************
     * Update TypeMap, equivalent to creating new instance with updated types.
     *
     * @param types The types in the machine learning force field.
     * @param subTypes The types in the actual structure.
     *
     * The condition `subTypes` @f$\subseteq@f$ `types` must be statisfied.
     **********************************************************************************************/
    void update(const Vec1String& types, const Vec1String& subTypes);
    /***********************************************************************************************
     * Mapping from structure type `subType` to type in force field.
     *
     * @param subType Type index in the structure.
     * @return Type index of `subType` in the machine learning force field.
     **********************************************************************************************/
    Int toType(const Int subType) const;
    /*******************************************************************************************
     * Mapping from structure type `subType` to type in force field.
     *
     * #return a Vec1Int which can be used to mimic same behviour as Int toType(const Int subType)
     *           const
     *
     * This function can be helpfull when using the TypeMap in combination with gpus
     * because only the vector and not the whole class have to be copied
     *******************************************************************************************/
    const Vec1Int& get_toType() const;
    /***********************************************************************************************
     * Convert type name (string representation) to force field type.
     *
     * @param name Input type name.
     * @return Force field type index corresponding to name, or, -1 if not present.
     **********************************************************************************************/
    Int toType(const String name) const;
    /***********************************************************************************************
     * Mapping from machine learned force field type to type in current structure.
     *
     * @param type Type index in the machine learning force field.
     * @return Subtype index of the input type index in the structure, or, -1 if the type is not in
     * the structure.
     **********************************************************************************************/
    Int toSubType(const Int type) const;
    /***********************************************************************************************
     * Convert type name (string representation) to structure type.
     *
     * @param name Input type name.
     * @return Structure type index corresponding to name, or, -1 if not present.
     **********************************************************************************************/
    Int toSubType(const String name) const;
    /***********************************************************************************************
     * Get type name for given input type.
     *
     * @param type Type index in the machine learning force field.
     * @return String representation of given type.
     **********************************************************************************************/
    String nameType(const Int type) const;
    /***********************************************************************************************
     * Get subtype name for given input subtype.
     *
     * @param subType Type index in the structure.
     * @return String representation of given subtype.
     **********************************************************************************************/
    String nameSubType(const Int subType) const;
    /***********************************************************************************************
     * Number of different atom types present in the current force field.
     *
     * @return Number of types (force field types).
     **********************************************************************************************/
    Size countForceFieldTypes() const;
    /***********************************************************************************************
     * Number of different atom types present in the current structure.
     *
     * @return Number of subtypes (structure types).
     **********************************************************************************************/
    Size countStructureTypes() const;
    /*******************************************************************************************
     * write types in type maps to the screen
     *******************************************************************************************/
    void printTypeMap() const;
    /*******************************************************************************************
     * @brief Retrieves the list of types.
     *
     * @return A constant reference to the vector of types.
     *******************************************************************************************/
    const Vec1String& get_types() const;
    /*******************************************************************************************
     * @brief Retrieves the list of subtypes.
     *
     * @return A constant reference to the vector of subtypes.
     *******************************************************************************************/
    const Vec1String& get_subTypes() const;

  private:
    /***********************************************************************************************
     * Create mapping between the types in the force field and the structure.
     **********************************************************************************************/
    void createMapping();
    /***********************************************************************************************
     * String representation (names) of force field types.
     **********************************************************************************************/
    Vec1String types;
    /***********************************************************************************************
     * String representation (names) of structure types.
     **********************************************************************************************/
    Vec1String subTypes;
    /***********************************************************************************************
     * Map from type index to subType index.
     **********************************************************************************************/
    Vec1Int mapToSubType;
    /***********************************************************************************************
     * Map from subtype index to type index.
     **********************************************************************************************/
    Vec1Int mapToType;
    /***********************************************************************************************
     * Map from name to type index.
     **********************************************************************************************/
    std::map<String, Int> mapNameToType;
    /***********************************************************************************************
     * Map from name to subtype index.
     **********************************************************************************************/
    std::map<String, Int> mapNameToSubType;
};

} // namespace vaspml

#endif
