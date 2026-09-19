#ifndef DEEPCOPY_HPP
#define DEEPCOPY_HPP

#include "Record.hpp" // IWYU pragma: keep
#include "RecordFunctor.hpp"
#include "types.hpp"

namespace vaspml::rec::detail
{

struct DeepCopy : public RecordFunctor
{
    DeepCopy(Record& record);

    template<typename T>
    void copy(const String& key, const T& item)
    {
        recordNew->put<T>(key, item);
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
    void up(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    bool shouldLoop(const String& key, Size arraySize, bool skipRecordArrays) override;

    /**********************************************************************************************
     * Pointer to new Record which should be filled with data from traversed one.
     *
     * If traverse() descends down into sub-Records this pointer will be replaced with newly
     * generated sub-Records. The "higher" Record is then temporarliy stored in `recordHistory` in
     * down() and restored in up().
     *
     * @note A raw pointer must be used here because it is initialized from `this` pointer.
     **********************************************************************************************/
    Record* recordNew;
    /// History of Record pointers required for descending in sub-Records.
    std::vector<Record*> recordHistory;
};

} // namespace vaspml::rec::detail

#endif
