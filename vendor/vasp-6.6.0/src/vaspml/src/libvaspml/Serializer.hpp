#ifndef SERIALIZER_HPP
#define SERIALIZER_HPP

#include "Item.hpp"
#include "RecordFunctor.hpp"
#include "buffer.hpp"
#include "types.hpp"

namespace vaspml
{
struct Record;
}

namespace vaspml::rec::detail
{

struct Serializer : public RecordFunctor
{
    Serializer(Buffer& buffer);

    template<typename T>
    void add(const String& key, const T& item)
    {
        // Add RECORD_ITEM marker to buffer to signal start of a new item.
        io::detail::write(ItemIndex::RECORD_ITEM, buffer);

        // Add key to buffer.
        io::detail::write(key, buffer);

        // Write item type marker to buffer.
        io::detail::write(itemIndex<T>(), buffer);

        // Add item to buffer.
        io::detail::write(item, buffer);

        return;
    }

    // clang-format off
    void operator()(const Record& record, bool entry) override;
    void operator()(const String& key, const Int&        item) override;
    void operator()(const String& key, const Real&       item) override;
    void operator()(const String& key, const String&     item) override;
    void operator()(const String& key, const bool&       item) override;
    void operator()(const String& key, const Vec1Real&   item) override;
    void operator()(const String& key, const Vec1Int&    item) override;
    void operator()(const String& key, const Vec1String& item) override;
    void operator()(const String& key, const Vec2Real&   item) override;
    void operator()(const String& key, const Vec2Int&    item) override;
    void operator()(const String& key, const Vec2String& item) override;
    // clang-format on

    void down(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    bool shouldLoop(const String& key, Size arraySize, bool skipRecordArrays) override;

    Buffer& buffer;
};

Size deserialize(Buffer const& buffer, Record& record, Size position = 0);

} // namespace vaspml::rec::detail

#endif
