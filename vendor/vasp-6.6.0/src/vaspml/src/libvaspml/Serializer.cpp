#include "Serializer.hpp"

#include "Record.hpp"
#include "debug.hpp"

#include <stdexcept> // for runtime_error

using namespace vaspml;
using namespace vaspml::rec::detail;
using namespace vaspml::io::detail;

Serializer::Serializer(Buffer& buffer) : buffer(buffer)
{}

void Serializer::operator()(const Record& /* record */, bool entry)
{
    ItemIndex marker;
    if (entry) marker = ItemIndex::RECORD_BEGIN;
    else marker = ItemIndex::RECORD_END;
    write(marker, buffer);

    return;
}

// clang-format off
void Serializer::operator()(const String& key, const Int&        item) { add(key, item); }
void Serializer::operator()(const String& key, const Real&       item) { add(key, item); }
void Serializer::operator()(const String& key, const String&     item) { add(key, item); }
void Serializer::operator()(const String& key, const bool&       item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec1Real&   item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec1Int&    item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec1String& item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec2Real&   item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec2Int&    item) { add(key, item); }
void Serializer::operator()(const String& key, const Vec2String& item) { add(key, item); }
// clang-format on

void Serializer::down(const String& key, Int arrayIndex, Size /* arraySize */)
{
    level++;
    this->key = key;

    // Single Record (ShRec):
    // ----------------------
    // Like any other item the following will be written before the actual data:
    // * RECORD_ITEM
    // * key
    // * SHREC (type marker)
    // Single Record is signaled by arrayIndex == -1.
    //
    // Vector of Records (Vec1ShRec) is treated in shouldLoop()!

    if (arrayIndex < 0)
    {
        // Add RECORD_ITEM marker to buffer to signal start of a new item.
        write(ItemIndex::RECORD_ITEM, buffer);

        // Add key to buffer.
        write(key, buffer);

        // Write item type marker to buffer.
        write(itemIndex<ShRec>(), buffer);
    }

    return;
}
bool Serializer::shouldLoop(const String& key, Size arraySize, bool /* skipRecordArrays */)
{
    // Vector of Records (Vec1ShRec):
    // ------------------------------
    // Only before the first Record the following will be written:
    // * RECORD_ITEM
    // * key
    // * VEC1SHREC (type marker)
    // * number of Record entries in vector

    // Add RECORD_ITEM marker to buffer to signal start of a new item.
    write(ItemIndex::RECORD_ITEM, buffer);

    // Add key to buffer.
    write(key, buffer);

    // Write item type marker to buffer.
    write(itemIndex<Vec1ShRec>(), buffer);

    // Add number of Records in vector to buffer.
    write(arraySize, buffer);

    return true;
}

Size rec::detail::deserialize(Buffer const& buffer, Record& record, Size position)
{
    ItemIndex mark = readItemIndex(buffer, position);
    VASPML_DEBUG_L1(
        if (mark != ItemIndex::RECORD_BEGIN)
        {
            throw std::runtime_error("ERROR: Expected begin of new record.");
        }
    );

    bool readNextItem = true;
    while (readNextItem)
    {
        mark = readItemIndex(buffer, position);
        VASPML_DEBUG_L1(
            if (!(mark == ItemIndex::RECORD_ITEM || mark == ItemIndex::RECORD_END))
            {
                throw std::runtime_error("ERROR: Expected item or end of record.");
            }
        );

        // If there are no more items, break here.
        if (mark == ItemIndex::RECORD_END) break;

        String key;
        read(buffer, key, position);

        // Now follows a marker determining the type of the item.
        mark = readItemIndex(buffer, position);

        // Create new key-value pair in Record.
        record.add(key, mark);
        // Read in data for new Record entry.
        // clang-format off
        if      (mark == ItemIndex::REAL)       read(buffer, record.get<Real>(key), position);
        else if (mark == ItemIndex::INT)        read(buffer, record.get<Int>(key), position);
        else if (mark == ItemIndex::STRING)     read(buffer, record.get<String>(key), position);
        else if (mark == ItemIndex::BOOL)       read(buffer, record.get<bool>(key), position);
        else if (mark == ItemIndex::VEC1REAL)   read(buffer, record.get<Vec1Real>(key), position);
        else if (mark == ItemIndex::VEC1INT)    read(buffer, record.get<Vec1Int>(key), position);
        else if (mark == ItemIndex::VEC1STRING) read(buffer, record.get<Vec1String>(key), position);
        else if (mark == ItemIndex::VEC2REAL)   read(buffer, record.get<Vec2Real>(key), position);
        else if (mark == ItemIndex::VEC2INT)    read(buffer, record.get<Vec2Int>(key), position);
        else if (mark == ItemIndex::VEC2STRING) read(buffer, record.get<Vec2String>(key), position);
        // clang-format on
        // If new entry is another Record, descend into recursive deserialize().
        else if (mark == ItemIndex::SHREC)
        {
            position = deserialize(buffer, record.dget<ShRec>(key), position);
        }
        else if (mark == ItemIndex::VEC1SHREC)
        {
            Size n;
            read(buffer, n, position);
            for (Size i = 0; i < n; ++i)
            {
                record.get<Vec1ShRec>(key).push_back(std::make_shared<Record>());
                position = deserialize(buffer, *(record.get<Vec1ShRec>(key).back()), position);
            }
        }
        else { throw std::runtime_error("ERROR: Unknown item type encountered."); }
    }

    return position;
}
