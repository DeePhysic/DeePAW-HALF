#ifndef STRINGALTERNATIVE_HPP
#define STRINGALTERNATIVE_HPP

#include "types.hpp"

namespace vaspml
{

/**************************************************************************************************
 * Helper class to deal with strings and their alternatives (e.g. current and deprecated tags).
 *
 * Here, we assume that each string has exactly one "main" variant. All other variants are
 * "alternatives". The main and alternative strings need to be set up via the add() functions.
 * Once set up, the operator overloads for () and [] may be used to get the main variant for a given
 * input string.
 **************************************************************************************************/
class StringAlternative
{
  public:
    explicit StringAlternative(ShRec record = nullptr);
    /**********************************************************************************************
     * Return main variant for alternative string (soft version).
     *
     * @param arg Given input string (may also be an unknown string or the main variant).
     * @return The main variant for a given input string or the input string for unknown strings.
     **********************************************************************************************/
    String operator()(const String& arg) const;
    /**********************************************************************************************
     * Return main variant for alternative string (strict version).
     *
     * @param arg Given input string (may also be the main variant).
     * @return The main variant for a given input string.
     *
     * @note This operator overload will throw an error if the input string is not in the list of
     * variants.
     **********************************************************************************************/
    String operator[](const String& arg) const;
    /**********************************************************************************************
     * Set up multiple entries for main and alternative strings.
     *
     * @param alternatives A 2-dimensional vector containing the string variants.
     *
     * See 1-dimensional version of add().
     **********************************************************************************************/
    void add(const Vec2String& alternatives);
    /**********************************************************************************************
     * Set up a single entry for one main and its alternative strings.
     *
     * @param alternatives A vector containing the string variants.
     *
     * @note The main string must be the first item in the input vector and should not be empty
     * (otherwise the list is ignored). Empty alternatives are ignored.
     **********************************************************************************************/
    void add(const Vec1String& alternatives);

  private:
    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /// Number of unique strings currently hosted.
    Int& numUnique;
    /// List of all strings including alternatives.
    Vec1String& astring;
    /// Indices assigned to all string alternatives (main string gets negative index).
    Vec1Int& index;

    String searchMain(Size i, String arg) const;
};

/**************************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **************************************************************************************************/
MapString dataMapStringAlternative();

} // namespace vaspml

#endif
