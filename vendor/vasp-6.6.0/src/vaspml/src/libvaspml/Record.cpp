#include "Record.hpp"

#include "DeepCopy.hpp"
#include "RecordFunctor.hpp"
#include "Tutor.hpp"

using namespace vaspml;

Record::Record(const Record& other)
{
    rec::detail::DeepCopy deepCopy(*this);
    traverse(other, deepCopy);
}

Record& Record::operator=(const Record& other)
{
    if (this != &other)
    {
        rec::detail::DeepCopy deepCopy(*this);
        traverse(other, deepCopy);
    }

    return *this;
}

Item& Record::operator[](String key)
{
    return dict[key];
}

Size Record::size() const
{
    return dict.size();
}

const Item& Record::at(String key) const
{
    return dict.at(key);
}

bool Record::contains(String key) const
{
    return dict.find(key) != dict.end();
}

bool Record::empty() const
{
    return dict.empty();
}

Vec1String Record::keys() const
{
    Vec1String vkeys(dict.size());
    std::transform(dict.begin(),
                   dict.end(),
                   vkeys.begin(),
                   [](auto const& entry) { return entry.first; });
    return vkeys;
}

ItemIndex Record::itemIndexOf(String key) const
{
    ItemIndex i = ItemIndex::UNDEFINED;
    if (contains(key)) i = static_cast<ItemIndex>(dict.at(key).index() + 1);

    return i;
}

String Record::typeOf(String key) const
{
    return SmartEnum<ItemIndex>::toString(itemIndexOf(key));
}

Record& Record::add(String key, ItemIndex type, bool allowKeyPresent)
{
    VASPML_DEBUG_L2(
        if (contains(key)
            && (!allowKeyPresent || typeOf(key) != SmartEnum<ItemIndex>::toString(type)))
        {
            global::tutor.warning("Overwriting existing Record entry \"" + key
                                  + "\", old type: " + typeOf(key)
                                  + ", new type: " + SmartEnum<ItemIndex>::toString(type) + ".");
        }
    );

    // Early exit in case an entry to that key already exists.
    if (contains(key) && allowKeyPresent) return *this;

    // clang-format off
    if      (type == ItemIndex::REAL)       dict[key] = Real();
    else if (type == ItemIndex::INT)        dict[key] = Int();
    else if (type == ItemIndex::STRING)     dict[key] = String();
    else if (type == ItemIndex::BOOL)       dict[key] = bool();
    else if (type == ItemIndex::SHREC)      dict[key] = std::make_shared<Record>();
    else if (type == ItemIndex::VEC1REAL)   dict[key] = Vec1Real();
    else if (type == ItemIndex::VEC1INT)    dict[key] = Vec1Int();
    else if (type == ItemIndex::VEC1STRING) dict[key] = Vec1String();
    else if (type == ItemIndex::VEC1SHREC)  dict[key] = Vec1ShRec();
    else if (type == ItemIndex::VEC2REAL)   dict[key] = Vec2Real();
    else if (type == ItemIndex::VEC2INT)    dict[key] = Vec2Int();
    else if (type == ItemIndex::VEC2STRING) dict[key] = Vec2String();
    // clang-format on
    else
    {
        global::tutor.warning("ERROR: Cannot add entry to Record, supplied type \""
                              + SmartEnum<ItemIndex>::toString(type) + "\" is not supported.");
    }

    return *this;
}

Record& Record::copyFrom(const Record& input, String key, bool allowKeyPresent)
{
    const ItemIndex type = input.itemIndexOf(key);
    add(key, type, allowKeyPresent);

    // clang-format off
    if      (type == ItemIndex::REAL)       get<Real>(key)       = input.cget<Real>(key);
    else if (type == ItemIndex::INT)        get<Int>(key)        = input.cget<Int>(key);
    else if (type == ItemIndex::STRING)     get<String>(key)     = input.cget<String>(key);
    else if (type == ItemIndex::BOOL)       get<bool>(key)       = input.cget<bool>(key);
    else if (type == ItemIndex::VEC1REAL)   get<Vec1Real>(key)   = input.cget<Vec1Real>(key);
    else if (type == ItemIndex::VEC1INT)    get<Vec1Int>(key)    = input.cget<Vec1Int>(key);
    else if (type == ItemIndex::VEC1STRING) get<Vec1String>(key) = input.cget<Vec1String>(key);
    else if (type == ItemIndex::VEC2REAL)   get<Vec2Real>(key)   = input.cget<Vec2Real>(key);
    else if (type == ItemIndex::VEC2INT)    get<Vec2Int>(key)    = input.cget<Vec2Int>(key);
    else if (type == ItemIndex::VEC2STRING) get<Vec2String>(key) = input.cget<Vec2String>(key);
    // clang-format on
    else if (type == ItemIndex::SHREC)
    {
        get<ShRec>(key) = std::make_shared<Record>(input.dcget<ShRec>(key));
    }
    else if (type == ItemIndex::VEC1SHREC)
    {
        Vec1ShCRec other = input.vcget<Vec1ShRec>(key);
        Vec1ShRec& my = get<Vec1ShRec>(key);
        my.resize(other.size());
        std::transform(other.cbegin(),
                       other.cend(),
                       my.begin(),
                       [](const ShCRec& rec) { return std::make_shared<Record>(*rec); });
    }
    else
    {
        global::tutor.warning("ERROR: Cannot copy entry to Record, supplied type \""
                              + SmartEnum<ItemIndex>::toString(type) + "\" is not supported.");
    }

    return *this;
}

Record& Record::merge(const Record& input, bool allowKeyPresent)
{
    for (const String& key : input.keys()) { copyFrom(input, key, allowKeyPresent); }

    return *this;
}

Record& Record::add(String key, String type, bool allowKeyPresent)
{
    add(key, SmartEnum<ItemIndex>::toEnum(type), allowKeyPresent);

    return *this;
}

void Record::erase(String key)
{
    VASPML_DEBUG_L2(
        if (!contains(key))
        {
            global::tutor.warning("Attempt to erase Record entry \"" + key
                                  + "\" which is not present.");
        }
    );

    dict.erase(key);

    return;
}

Record vaspml::makeRecord(const MapString keyTypeMap)
{
    Record d;
    for (auto& [key, typeString] : keyTypeMap) d.add(key, typeString);

    return d;
}

ShRec vaspml::assignOrMakeRecord(ShRec in, const MapString keyTypeMap)
{
    if (in == nullptr) return std::make_shared<Record>(makeRecord(keyTypeMap));
    else
    {
        for (auto& [key, typeString] : keyTypeMap)
        {
            if (in->contains(key))
            {
                if (in->typeOf(key) != typeString)
                {
                    global::tutor.bug("Conflicting Record entries for key \"" + key
                                      + "\", exiting entry type: " + in->typeOf(key)
                                      + ", desired type: " + typeString + ".");
                }
            }
            else in->add(key, typeString);
        }
        return in;
    }
}
