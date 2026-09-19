#ifndef RECORD_HPP
#define RECORD_HPP

#include "Item.hpp"
#include "SmartEnum.hpp"
#include "debug.hpp"
#include "meta.hpp"
#include "types.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace vaspml
{

using RecordMap = std::map<String, Item>;

/**************************************************************************************************
 * Generic data structure used to store most data of VASPml.
 *
 * A Record contains a dictionary from string keys to vaspml::Item which is a specialization of
 * `std::variant`. Item contains only the common types we use in vaspml (e.g. Vec1Real, Vec2Int).
 * Two special entries are ShRec and Vec1ShRec, i.e., a single or array of smart pointers to another
 * instance of a Record. In this way records can be nested. Basically this data structure can be
 * thought of an in-memory representation of a JSON file (with limitations regarding array nesting).
 *
 * The Record class adds to this dictionary functionality to safely read and write access the data.
 * To add a new Record entry use the bracket notation equal to `std::map`.
 *
 * ~~~{.cpp}
 * Record myRecord;
 * myRecord["description"] = String{"This is my first Record entry"};
 * myRecord["age"] = 99;
 * myRecord["size"] = 2.00;
 * ~~~
 *
 * @warning Always wrap literal strings with `String{...}`, otherwise the compiler will try to
 * convert them to integer or bool type! It is best practice to be as explicit as possible to avoid
 * automatic type determination issues. For example, use `2.0` or `Real(2)`, otherwise `2` will be
 * stored in the Item as integer.
 *
 * Using the bracket access operator "[]" does not provide any form of checking whether the given
 * key-value pair already exists and would eventually be overwritten. For safe alternatives use:
 *
 * * add(key, type) ... Adds an entry with given type (ItemIndex or string representation).
 * * put(key, value) ... Directly puts a value for given key in Record.
 *
 * These two functions will throw a warning if an existing key would be touched by the operation.
 * An additional third argument allows to silence the warning in case the entry already exists and
 * has the desired type. A warning is always given if the existing and desired type do not match.
 *
 * Depending on the underlying type use the following functions to retrieve data:
 *
 * Type        |  get()  |  cget() |  dget() | dcget() | vcget() |
 * ------------|:-------:|:-------:|:-------:|:-------:|:-------:|
 * Real        |    x    |    x    |         |         |         |
 * Int         |    x    |    x    |         |         |         |
 * String      |    x    |    x    |         |         |         |
 * bool        |    x    |    x    |         |         |         |
 * ShRec       |    x    |    x    |    x    |    x    |         |
 * Vec1Real    |    x    |    x    |         |         |         |
 * Vec1Int     |    x    |    x    |         |         |         |
 * Vec1String  |    x    |    x    |         |         |         |
 * Vec1ShRec   |    x    |         |         |         |    x    |
 * Vec2Real    |    x    |    x    |         |         |         |
 * Vec2Int     |    x    |    x    |         |         |         |
 * Vec2String  |    x    |    x    |         |         |         |
 *
 * Here, the variations of the get() functions have the following meaning:
 *
 * * get() ... read/write access to data,
 * * cget() ... read-only access to data,
 * * dget() ... read/write access to dereferenced data (for `shared_ptr` content),
 * * dcget() ... read-only access to dereferenced data (for `shared_ptr` content).
 * * vcget() ... special operation for Vec1ShRec, returns vector with pointers to read-only records.
 *
 * @attention Avoid using the `std::get<>()` function directly on the dictionary values because this
 * allows to accidentially modify data which is meant to be constant.
 *
 * **Examples:**
 *
 * ~~~{.cpp}
 * // Create a record intended to hold a data set of structures.
 * Record dataset;
 *
 * // Add basic data set information.
 * dataset["_type"] = String("dataset");
 * dataset["numStructures"] = 2;
 * dataset["types"] = Vec1String{"Fe", "O"};
 *
 * // Create a vector of "sub"-records for the structures.
 * dataset["structures"] = Vec1ShRec(2, std::make_shared<Record>());
 *
 * // Fill the first structure with some data.
 * Record& s1 = *dataset.get<Vec1ShRec>("structures")[0];
 * s1["lattice"] = Vec2Real{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
 * s1["positions"] = Vec1Real{0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5};
 *
 * // Fill the second structure with some data.
 * Record& s2 = *dataset.get<Vec1ShRec>("structures")[1];
 * s2["lattice"] = Vec2Real{{1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
 * s2["positions"] = Vec1Real{0.5, 0.5, 0.3, 0.3, 0.0, 0.3, 0.0, 0.7, 0.1};
 *
 * // Access data and print to screen.
 * const Int& n = dataset.cget<Int>("numStructures");
 * //    Int& m = dataset.cget<Int>("numStructures"); // Does not compile, cget() returns const!
 * std::cout << "Record contains " << n << " structures with elements \"";
 * for (auto e : dataset.cget<Vec1String>("types")) std::cout << e << " ";
 * std::cout << "\": \n";
 * Vec1ShCRec structures = dataset.vcget<Vec1ShRec>("structures"); // Read-only access.
 * for (Int i = 0; i < n; ++i)
 * {
 *     std::cout << "Structure " << i << ": positions = ";
 *     for (auto p : structures[i]->cget<Vec1Real>("positions")) std::cout << p << " ";
 *     std::cout << "\n";
 * }
 * ~~~
 *
 * @note The getter function get() and dget() by design cannot be called for instances of `const
 * Record`, only "const"-Versions (cget(), dcget(), vcget()) can be used in this case. This helps
 * enforcing const-correctness.
 *
 * Records can be transformed into JSON and serialized binary buffer like this:
 * ~~~
 *
 * +----------+                                        +----------+
 * |          |     String toJson(...)                 |          |
 * |          |     void toJson(..., ostream&)         |          |
 * |          | -------------------------------------> |   JSON   |
 * |          | <------------------------------------- |          |
 * |          |     void fromJson(String&, ...)        |          |
 * |          |                                        +----------+
 * |  Record  |
 * |          |                                        +----------+
 * |          |     void toBuffer(..., Buffer&)        |          |
 * |          | -------------------------------------> |  Binary  |
 * |          | <------------------------------------- |   blob   |
 * |          |     void fromBuffer(Buffer&, ...)      |          |
 * |          |                                        +----------+
 * +----------+
 *
 * ~~~
 **************************************************************************************************/
struct Record
{
    /// Dictionary holding all data in this Record.
    RecordMap dict;

    Record() = default;
    /// Copy constructor (uses DeepCopy functor to perform recursive deep copy).
    Record(const Record& other);
    /// Move constructor.
    Record(Record&& other) noexcept = default;
    /// Copy assignment operator (uses DeepCopy functor to perform recursive deep copy).
    Record& operator=(const Record& other);
    /// Move assignment operator.
    Record& operator=(Record&& other) noexcept = default;
    /// Forwards bracket operator of dictionary to Record.
    Item& operator[](String key);
    /// Get number of entries in this Record.
    Size size() const;
    /// Forwards at() member function of dictionary to Record.
    const Item& at(String key) const;
    /// Check if key is present in Record.
    bool contains(String key) const;
    /// Check if key is present in Record.
    bool empty() const;
    /// Return list of Record keys.
    Vec1String keys() const;
    /// Return ItemIndex representation of type of entry with given key.
    ItemIndex itemIndexOf(String key) const;
    /// Return string representation of type of entry with given key.
    String typeOf(String key) const;
    /// Add empty entry via ItemIndex type.
    Record& add(String key, ItemIndex type, bool allowKeyPresent = false);
    /// Add empty entry via type string.
    Record& add(String key, String type, bool allowKeyPresent = false);
    /// Copy in entry of given Record with given key.
    Record& copyFrom(const Record& input, String key, bool allowKeyPresent = false);
    /// Merge in all entries of given Record.
    Record& merge(const Record& input, bool allowKeyPresent = false);
    /// Add empty entry via template function (not allowed for shared pointers and vectors thereof).
    template<typename T,
             typename std::enable_if_t<std::conjunction_v<std::negation<isSharedPtr<T>>,
                                                          std::negation<isVectorSharedPtr<T>>>>* =
                 nullptr>
    Record& put(String key, T value = T(), bool allowKeyPresent = false)
    {
        add(key, itemIndex<T>(), allowKeyPresent);
        get<T>(key) = value;

        return *this;
    }
    /// Remove entry from Record.
    void erase(String key);
    /// Wrapper around std::get<>().
    template<typename T>
    T& get(String key)
    {
        VASPML_DEBUG_L1(checkKey<T>(key););
        return std::get<T>(dict[key]);
    }
    /// For `shared_ptr`-type content only, return reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto& dget(String key)
    {
        VASPML_DEBUG_L2(checkKey<T>(key););
        return *std::get<T>(dict[key]);
    }
    /// For non-`shared_ptr` and non-`vector<shared_ptr>` content only, return const reference to
    /// underlying data.
    template<typename T,
             typename std::enable_if_t<std::conjunction_v<std::negation<isSharedPtr<T>>,
                                                          std::negation<isVectorSharedPtr<T>>>>* =
                 nullptr>
    T const& cget(String key) const
    {
        VASPML_DEBUG_L2(checkKey<T>(key););
        return std::get<T>(dict.at(key));
    }
    /**********************************************************************************************
     * For `shared_ptr`-type content only, return `shared_ptr` to const version of underlying data.
     *
     * @note Creates a new `shared_ptr`, hence increases `use_count()`. Prefer dcget<>() if
     * possible.
     **********************************************************************************************/
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    auto cget(String key) const
    {
        VASPML_DEBUG_L2(checkKey<T>(key););
        using constT = const typename std::pointer_traits<T>::element_type;
        return std::shared_ptr<constT>(std::get<T>(dict.at(key)));
    }
    /// For `shared_ptr`-type content only, return const reference to underlying data.
    template<typename T, typename std::enable_if_t<isSharedPtr<T>::value>* = nullptr>
    const typename std::pointer_traits<T>::element_type& dcget(String key) const
    {
        VASPML_DEBUG_L2(checkKey<T>(key););
        return *std::get<T>(dict.at(key));
    }
    /**********************************************************************************************
     * Special getter function for Vec1ShRec Record entries.
     *
     * This getter is used to obtain a "const"-version of an Vec1ShRec Record entry. Background: a
     * naive version of such a getter could return a "const Vec1ShRec". However, this would still
     * allow to change the Record behind the shared pointers which we do NOT want! Hence, this
     * function instead creates a new vector of shared pointers to "const"-casted Record entries.
     *
     * @return Vector of shared pointers to const Records.
     **********************************************************************************************/
    template<typename T, typename std::enable_if_t<isVectorSharedPtr<T>::value>* = nullptr>
    auto vcget(String key) const
    {
        if (!std::holds_alternative<Vec1ShRec>(dict.at(key)))
        {
            throw std::runtime_error("ERROR: Item mapped to dictionary entry \"" + key
                                     + "\" does not hold a vector of shared pointers,"
                                     + "cannot execute vcget() function.");
        }
        VASPML_DEBUG_L2(checkKey<T>(key););
        // Unpack types, e.g. from Vec1ShRec extract Record.
        using sharedPtrT = typename T::value_type;
        using baseT = typename std::pointer_traits<sharedPtrT>::element_type;
        // Version of template type which is const at innermost level, e.g. ShCRec.
        using constT = std::shared_ptr<const baseT>;
        T const&            in = std::get<T>(dict.at(key));
        std::vector<constT> out(in.size());
        std::copy(in.begin(), in.end(), out.begin());
        return out;
    }
    /**********************************************************************************************
     * Check if key already exists in Record.
     **********************************************************************************************/
    template<typename T>
    void checkKey(const String& key) const
    {
        try
        {
            std::get<T>(dict.at(key));
        }
        catch (const std::out_of_range& e)
        {
            const String message =
                "Key '" + key + "' not found in Record. Please check supplied key.\n";
            throw std::runtime_error("\033[31m" + message + "\033[0m" + e.what());
        }
        catch (const std::bad_variant_access& e)
        {
            const String    type = typeOf(key);
            const ItemIndex indx = itemIndex<T>();
            const String    supplied = SmartEnum<ItemIndex>::toString(indx);
            const String    message =
                "For key '" + key
                + "' wrong template argument was supplied. The supplied template argument is '"
                + supplied + "' and the actual template argument is '" + type
                + "'.\nPlease supply proper template argument.\n";
            throw std::runtime_error("\033[31m" + message + "\033[0m" + e.what());
        }
        return;
    }
};

/**************************************************************************************************
 * Create Record with entries according to input type map.
 *
 * @param keyTypeMap A map from String to String containing the desired (key, type)-pairs.
 *
 * The type strings correspond to the ones defined in SmartEnum<ItemIndex>::names (see Record.cpp).
 * The key-type map is typically provided as a literal, e.g.
 *
 * @code
 * Record r = makeRecord({{"natoms", "Int"}, {"positions", "Vec2Real"});
 * @endcode
 **************************************************************************************************/
Record makeRecord(const MapString keyTypeMap);
/***************************************************************************************************
 * Helper function which either returns input or creates new shared pointer to Record.
 *
 * @param in Input shared pointer or `nullptr`.
 *
 * If `in` is the null pointer a new object and a `shared_ptr` to it is created with
 * makeRecord(). Otherwise, the input pointer is just returned. Use in constructor for classes
 * which accept their memory to be managed from outside, e.g.:
 *
 * @code
 * MyClass::MyClass(inputRecord) :
 *     myRecord(assignOrMakeRecord(inputRecord,
 *                                 {
 *                                     {"myInt", "Int"     },
 *                                     {"myVec", "Vec1Real"}
 * })) {};
 * @endcode
 **************************************************************************************************/
ShRec assignOrMakeRecord(ShRec in, const MapString keyTypeMap);

} // namespace vaspml

#endif
