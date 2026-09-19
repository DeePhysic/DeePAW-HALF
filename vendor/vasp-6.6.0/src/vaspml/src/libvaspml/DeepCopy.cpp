#include "DeepCopy.hpp"

#include "Item.hpp"
#include "Record.hpp"

using namespace vaspml::rec::detail;

DeepCopy::DeepCopy(Record& record) : recordNew(&record)
{}

// clang-format off
void DeepCopy::operator()(const Record& /*record*/, bool /*entry*/) {}
void DeepCopy::operator()(const String& key, const Int&        item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Real&       item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const String&     item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const bool&       item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec1Real&   item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec1Int&    item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec1String& item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec2Real&   item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec2Int&    item) { copy(key, item); }
void DeepCopy::operator()(const String& key, const Vec2String& item) { copy(key, item); }
// clang-format on

void DeepCopy::down(const String& key, Int arrayIndex, Size /*arraySize*/)
{
    if (arrayIndex < 0)
    {
        recordNew->add(key, ItemIndex::SHREC);
        recordHistory.push_back(&recordNew->dget<ShRec>(key));
    }
    else
    {
        recordNew->get<Vec1ShRec>(key).push_back(std::make_shared<Record>());
        recordHistory.push_back(&*recordNew->get<Vec1ShRec>(key).back());
    }
    std::swap(recordNew, recordHistory.back());

    return;
}

void DeepCopy::up(const String& /*key*/, Int /*arrayIndex*/, Size /*arraySize*/)
{
    std::swap(recordNew, recordHistory.back());
    recordHistory.pop_back();

    return;
}

bool DeepCopy::shouldLoop(const String& key, Size /* arraySize */, bool /* skipRecordArrays */)
{
    recordNew->add(key, ItemIndex::VEC1SHREC);

    return true;
}
