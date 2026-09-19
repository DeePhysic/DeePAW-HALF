#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstddef> // IWYU pragma: export
#include <cstdint> // IWYU pragma: export
#include <map>     // IWYU pragma: export
#include <memory>  // IWYU pragma: export
#include <string>  // IWYU pragma: export
#include <vector>  // IWYU pragma: export

namespace vaspml
{

// Basic types.
using Real = double;
using Int = int;
using Int64 = int64_t;
using UInt = unsigned int;
using UInt64 = uint64_t;
using String = std::string;
using Size = std::size_t;

// Forward declare Record.
struct Record;
using ShRec = std::shared_ptr<Record>;
using ShCRec = std::shared_ptr<const Record>;
using Vec1ShRec = std::vector<ShRec>;
using Vec1ShCRec = std::vector<ShCRec>;

// Vectors of basic types.
using Vec1Real = std::vector<Real>;
using Vec2Real = std::vector<std::vector<Real>>;
using Vec3Real = std::vector<std::vector<std::vector<Real>>>;

using Vec1Int = std::vector<Int>;
using Vec2Int = std::vector<std::vector<Int>>;
using Vec3Int = std::vector<std::vector<std::vector<Int>>>;

using Vec1Int64 = std::vector<Int64>;
using Vec2Int64 = std::vector<std::vector<Int64>>;
using Vec3Int64 = std::vector<std::vector<std::vector<Int64>>>;

using Vec1UInt = std::vector<UInt>;
using Vec2UInt = std::vector<std::vector<UInt>>;
using Vec3UInt = std::vector<std::vector<std::vector<UInt>>>;

using Vec1UInt64 = std::vector<UInt64>;
using Vec2UInt64 = std::vector<std::vector<UInt64>>;
using Vec3UInt64 = std::vector<std::vector<std::vector<UInt64>>>;

using Vec1Size = std::vector<Size>;
using Vec2Size = std::vector<std::vector<Size>>;
using Vec3Size = std::vector<std::vector<std::vector<Size>>>;

using Vec1String = std::vector<String>;
using Vec2String = std::vector<std::vector<String>>;
using Vec3String = std::vector<std::vector<std::vector<String>>>;

// Shared pointers to vectors of basic types.
using ShVec1Real = std::shared_ptr<Vec1Real>;
using ShVec2Real = std::shared_ptr<Vec2Real>;
using ShVec3Real = std::shared_ptr<Vec3Real>;

using ShVec1Int = std::shared_ptr<Vec1Int>;
using ShVec2Int = std::shared_ptr<Vec2Int>;
using ShVec3Int = std::shared_ptr<Vec3Int>;

using ShVec1Int64 = std::shared_ptr<Vec1Int64>;
using ShVec2Int64 = std::shared_ptr<Vec2Int64>;
using ShVec3Int64 = std::shared_ptr<Vec3Int64>;

using ShVec1UInt = std::shared_ptr<Vec1UInt>;
using ShVec2UInt = std::shared_ptr<Vec2UInt>;
using ShVec3UInt = std::shared_ptr<Vec3UInt>;

using ShVec1UInt64 = std::shared_ptr<Vec1UInt64>;
using ShVec2UInt64 = std::shared_ptr<Vec2UInt64>;
using ShVec3UInt64 = std::shared_ptr<Vec3UInt64>;

using ShVec1String = std::shared_ptr<Vec1String>;
using ShVec2String = std::shared_ptr<Vec2String>;
using ShVec3String = std::shared_ptr<Vec3String>;

// Const content versions of shared pointers used in MultiType.
using ShCVec1Real = std::shared_ptr<const Vec1Real>;
using ShCVec1Int = std::shared_ptr<const Vec1Int>;
using ShCVec1String = std::shared_ptr<const Vec1String>;
using ShCVec2Real = std::shared_ptr<const Vec2Real>;
using ShCVec2Int = std::shared_ptr<const Vec2Int>;

// Raw bytes and corresponding buffer type.
using Byte = std::byte;
using Buffer = std::vector<Byte>;

// Shortcuts for various predefined maps.
using MapString = std::map<String, String>;

} //namespace vaspml

#endif
