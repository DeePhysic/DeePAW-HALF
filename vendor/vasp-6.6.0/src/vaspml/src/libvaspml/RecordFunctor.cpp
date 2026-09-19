#include "RecordFunctor.hpp"

#include "Record.hpp"

#include <algorithm>   // for find, sort, find_if_not
#include <iterator>    // for distance, __iterator_category
#include <type_traits> // for decay_t, is_same_v
#include <variant>     // for visit

using namespace vaspml::rec::detail;

// clang-format off
void RecordFunctor::operator()(const Record& /* record */, bool /* entry */) {}
void RecordFunctor::operator()(const String& /* key */, const Real&       /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Int&        /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const String&     /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const bool&       /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1Real&   /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1Int&    /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec1String& /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec2Real&   /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec2Int&    /* item */) {}
void RecordFunctor::operator()(const String& /* key */, const Vec2String& /* item */) {}
// clang-format on

void RecordFunctor::down(const String& key, Int /* arrayIndex */, Size /* arraySize */)
{
    level++;
    this->key = key;
    //std::cout << "Record \"" << key << "\"";
    //if (arrayIndex >= 0) std::cout << "[" << arrayIndex << "]";
    //std::cout << " (level " << level << ")\n";

    return;
}

void RecordFunctor::up(const String& /*key*/, Int /* arrayIndex */, Size /* arraySize */)
{
    level--;
    //std::cout << "End Record \"" << key << "\"\n";

    return;
}

void RecordFunctor::next(const String& /*key*/, bool /*last*/)
{
    return;
}

bool RecordFunctor::shouldLoop(const String& /* key */,
                               Size /* arraySize */,
                               bool /* skipRecordArrays */)
{
    return true;
}

bool RecordFunctor::shouldDescend(const Record& /* record */,
                                  const String& /* key */,
                                  Int /* arrayIndex */,
                                  Size /* arraySize */)
{
    return true;
}

void RecordFunctor::postLoop()
{
    return;
}

void vaspml::rec::detail::traverse(const Record& record, RecordFunctor& f, bool skipRecordArrays)
{
    f(record, true);
    Vec1String keys = record.keys();
    // Sort keys such that metadata comes first, otherwise alphabetically.
    std::sort(keys.begin(),
              keys.end(),
              //[](const String& lhs, const String& rhs)
              [](const String& lhs, const String& rhs)
              {
                  if (lhs[0] == '_' && rhs[0] != '_') return true;
                  else if (lhs[0] != '_' && rhs[0] == '_') return false;
                  else return lhs < rhs;
              });
    // If there is a specific sort order given in the metadata use that to sort further.
    if (record.contains("_sort"))
    {
        Vec1String order = record.cget<Vec1String>("_sort");
        // Find start of non-metadata keys.
        Vec1String::iterator start =
            std::find_if_not(keys.begin(), keys.end(), [](String& s) { return s[0] == '_'; });
        // Sort all non-metadata keys.
        std::sort(start,
                  keys.end(),
                  [&order](const String& lhs, const String& rhs)
                  {
                      // Determine position in sort order list (gives end() if not found).
                      Vec1String::const_iterator plhs = std::find(order.begin(), order.end(), lhs);
                      Vec1String::const_iterator prhs = std::find(order.begin(), order.end(), rhs);
                      // Since two keys cannot be the same, identical iterator means both were not
                      // found. Hence, sort again alphabetically.
                      if (plhs == prhs) return lhs < rhs;
                      // Otherwise, if one or both were found compare with their position in sort
                      // order list.
                      else return plhs < prhs;
                  });
    }
    const auto last = keys.cend() - 1;
    for (auto k = keys.cbegin(); k != keys.cend(); ++k)
    {
        std::visit(
            [&](auto&& arg)
            {
                using T = std::decay_t<decltype(arg)>;
                //******************
                // Item types: ShRec
                //******************
                if constexpr (std::is_same_v<T, ShRec>)
                {
                    if (f.shouldDescend(*arg, *k))
                    {
                        f.down(*k);
                        traverse(*arg, f);
                        f.up(*k);
                    }
                }
                //**********************
                // Item types: Vec1ShRec
                //**********************
                else if constexpr (std::is_same_v<T, Vec1ShRec>)
                {
                    Vec1ShRec::const_iterator begin = arg.cbegin();
                    Vec1ShRec::const_iterator end = arg.cend();
                    Size                      arraySize = arg.size();
                    if (!arg.empty() && skipRecordArrays)
                    {
                        end = arg.cbegin() + 1;
                        arraySize = 1;
                    }
                    if (f.shouldLoop(*k, arraySize, skipRecordArrays))
                    {
                        for (auto a1 = begin; a1 != end; ++a1)
                        {
                            Int arrayIndex = std::distance(begin, a1);
                            if (f.shouldDescend(**a1, *k, arrayIndex, arraySize))
                            {
                                f.down(*k, arrayIndex, arraySize);
                                traverse(**a1, f);
                                f.up(*k, arrayIndex, arraySize);
                            }
                        }
                        f.postLoop();
                    }
                }
                //******************
                // Actual data types
                //******************
                else f(*k, arg);
            },
            record.at(*k));
        f.next(*k, k == last);
    }
    f(record, false);

    return;
}
