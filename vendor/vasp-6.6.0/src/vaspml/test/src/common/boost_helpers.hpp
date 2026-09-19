#ifndef BOOST_HELPERS_HPP
#define BOOST_HELPERS_HPP

#include "Item.hpp"
#include "Record.hpp"
#include "SmartEnum.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <boost/test/unit_test.hpp>

#include <iterator>    // for distance, begin, end
#include <limits>      // for numeric_limits
#include <type_traits> // for is_same_v

/** Compare two iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<double>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS(aa, bb, tolerance, containerName)                          \
    {                                                                                        \
        using std::distance;                                                                 \
        using std::begin;                                                                    \
        using std::end;                                                                      \
        auto a = begin(aa), a0 = begin(aa), ae = end(aa);                                    \
        auto b = begin(bb);                                                                  \
        BOOST_TEST_INFO("Vector \"" << containerName << "\" length comparison");             \
        BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb)));                          \
        for (; a != ae; ++a, ++b)                                                            \
        {                                                                                    \
            BOOST_TEST_INFO("Vector \"" << containerName << "\" element " << distance(a0, a) \
                                        << " comparison: " << *a << " vs. " << *b);          \
            BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                         \
        }                                                                                    \
    }

/** Compare two two-dimensional iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<double>>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS_2D(aaa, bbb, tolerance, containerName)                          \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aa = begin(aaa), aa0 = begin(aaa), aae = end(aaa);                                   \
        auto bb = begin(bbb);                                                                     \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(bbb)));                           \
        for (; aa != aae; ++aa, ++bb)                                                             \
        {                                                                                         \
            auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                                  \
            auto b = begin(*bb);                                                                  \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa)     \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                          \
            for (; a != ae; ++a, ++b)                                                             \
            {                                                                                     \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa) \
                                           << "," << distance(a0, a) << " comparison: " << *a     \
                                           << " vs. " << *b);                                     \
                BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                          \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare two three-dimensional iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<std::vector<double>>>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS_3D(aaaa, bbbb, tolerance, containerName)                        \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aaa = begin(aaaa), aaa0 = begin(aaaa), aaae = end(aaaa);                             \
        auto bbb = begin(bbbb);                                                                   \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aaa, aaae), distance(bbb, end(bbbb)));                       \
        for (; aaa != aaae; ++aaa, ++bbb)                                                         \
        {                                                                                         \
            auto aa = begin(*aaa), aa0 = begin(*aaa), aae = end(*aaa);                            \
            auto bb = begin(*bbb);                                                                \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aaa0, aaa)   \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(*bbb)));                      \
            for (; aa != aae; ++aa, ++bb)                                                         \
            {                                                                                     \
                auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                              \
                auto b = begin(*bb);                                                              \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element "                      \
                                           << distance(aaa0, aaa) << "," << distance(aa0, aa)     \
                                           << " 3rd dimension comparison");                       \
                BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                      \
                for (; a != ae; ++a, ++b)                                                         \
                {                                                                                 \
                    BOOST_TEST_INFO("Array \"" << containerName << "\" element "                  \
                                               << distance(aaa0, aaa) << "," << distance(aa0, aa) \
                                               << "," << distance(a0, a) << " comparison: " << *a \
                                               << " vs. " << *b);                                 \
                    BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                      \
                }                                                                                 \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare equality of two iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<int>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS(aa, bb, containerName)                                     \
    {                                                                                        \
        using std::distance;                                                                 \
        using std::begin;                                                                    \
        using std::end;                                                                      \
        auto a = begin(aa), a0 = begin(aa), ae = end(aa);                                    \
        auto b = begin(bb);                                                                  \
        BOOST_TEST_INFO("Vector \"" << containerName << "\" length comparison");             \
        BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb)));                          \
        for (; a != ae; ++a, ++b)                                                            \
        {                                                                                    \
            BOOST_TEST_INFO("Vector \"" << containerName << "\" element " << distance(a0, a) \
                                        << " comparison: " << *a << " vs. " << *b);          \
            BOOST_REQUIRE_EQUAL(*a, *b);                                                     \
        }                                                                                    \
    }

/** Compare equality of two two-dimensional iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<int>>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS_2D(aaa, bbb, containerName)                                     \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aa = begin(aaa), aa0 = begin(aaa), aae = end(aaa);                                   \
        auto bb = begin(bbb);                                                                     \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(bbb)));                           \
        for (; aa != aae; ++aa, ++bb)                                                             \
        {                                                                                         \
            auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                                  \
            auto b = begin(*bb);                                                                  \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa)     \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                          \
            for (; a != ae; ++a, ++b)                                                             \
            {                                                                                     \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa) \
                                           << "," << distance(a0, a) << " comparison: " << *a     \
                                           << " vs. " << *b);                                     \
                BOOST_REQUIRE_EQUAL(*a, *b);                                                      \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare equality of two three-dimensional iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<std::vector<int>>>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS_3D(aaaa, bbbb, containerName)                                   \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aaa = begin(aaaa), aaa0 = begin(aaaa), aaae = end(aaaa);                             \
        auto bbb = begin(bbbb);                                                                   \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aaa, aaae), distance(bbb, end(bbbb)));                       \
        for (; aaa != aaae; ++aaa, ++bbb)                                                         \
        {                                                                                         \
            auto aa = begin(*aaa), aa0 = begin(*aaa), aae = end(*aaa);                            \
            auto bb = begin(*bbb);                                                                \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aaa0, aaa)   \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(*bbb)));                      \
            for (; aa != aae; ++aa, ++bb)                                                         \
            {                                                                                     \
                auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                              \
                auto b = begin(*bb);                                                              \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element "                      \
                                           << distance(aaa0, aaa) << "," << distance(aa0, aa)     \
                                           << " 3rd dimension comparison");                       \
                BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                      \
                for (; a != ae; ++a, ++b)                                                         \
                {                                                                                 \
                    BOOST_TEST_INFO("Array \"" << containerName << "\" element "                  \
                                               << distance(aaa0, aaa) << "," << distance(aa0, aa) \
                                               << "," << distance(a0, a) << " comparison: " << *a \
                                               << " vs. " << *b);                                 \
                    BOOST_REQUIRE_EQUAL(*a, *b);                                                  \
                }                                                                                 \
            }                                                                                     \
        }                                                                                         \
    }

/// Check equality with context message.
#define CHECK_EQUAL_MESSAGE(L, R, M) \
    {                                \
        BOOST_TEST_INFO(M);          \
        BOOST_CHECK_EQUAL(L, R);     \
    }
/// Check equality with context message (warning).
#define WARN_EQUAL_MESSAGE(L, R, M) \
    {                               \
        BOOST_TEST_INFO(M);         \
        BOOST_WARN_EQUAL(L, R);     \
    }
/// Check equality with context message (required).
#define REQUIRE_EQUAL_MESSAGE(L, R, M) \
    {                                  \
        BOOST_TEST_INFO(M);            \
        BOOST_REQUIRE_EQUAL(L, R);     \
    }

/// Check if expression throws and what() message contains expected substring.
#define REQUIRE_EXCEPTION_WHAT(expression, exceptionType, partOfWhat)              \
    {                                                                              \
        BOOST_TEST_INFO(std::string("Exception what() message must contain: \"")   \
                        + partOfWhat + "\".");                                     \
        auto exceptionWhatContains = [](const std::exception& e)                   \
        {                                                                          \
            std::string what = e.what();                                           \
            bool contains = what.find(partOfWhat) != what.npos;                    \
            BOOST_TEST_INFO(std::string("Exception what() message is: \"")         \
                            + what + "\".");                                       \
            return contains;                                                       \
        };                                                                         \
        BOOST_REQUIRE_EXCEPTION(expression, exceptionType, exceptionWhatContains); \
    }

template<typename T>
void REQUIRE_RECORD_ENTRY(const vaspml::Record& record, std::string key, const T& expected)
{
    using namespace vaspml;
    BOOST_TEST_INFO("Testing if record contains key: \"" + key + "\"");
    BOOST_REQUIRE(record.contains(key));
    vaspml::ItemIndex type = itemIndex<T>();
    BOOST_TEST_INFO("Testing if record item for key \"" + key + "\" has type: \"" + toString(type)
                    + "\"");
    BOOST_REQUIRE_EQUAL(record.itemIndexOf(key), type);

    Real tol = 10 * std::numeric_limits<Real>::epsilon();
    if constexpr(std::is_same_v<T, Real>)
    {
        BOOST_TEST_INFO("Testing if record item for key \"" + key + "\" is close to: "
                        + str("%24.16E", expected) + " (tolerance: " + str("%24.16E", tol) + ")");
        BOOST_REQUIRE_SMALL(record.cget<Real>(key) - expected, tol);
    }
    else if constexpr(std::is_same_v<T, Int>)
    {
        BOOST_TEST_INFO("Testing if record item for key \"" + key
                        + "\" has integer value: " + std::to_string(expected));
        BOOST_REQUIRE_EQUAL(record.cget<Int>(key), expected);
    }
    else if constexpr(std::is_same_v<T, String>)
    {
        BOOST_TEST_INFO("Testing if record item for key \"" + key
                        + "\" is string: \"" + expected + "\"");
        BOOST_REQUIRE_EQUAL(record.cget<String>(key), expected);
    }
    else if constexpr(std::is_same_v<T, bool>)
    {
        BOOST_TEST_INFO("Testing if record item for key \"" + key
                        + "\" has boolean value: " + (expected ? "true" : "false"));
        BOOST_REQUIRE_EQUAL(record.cget<bool>(key), expected);
    }
    else if constexpr(std::is_same_v<T, Vec1Real>)
    {
        REQUIRE_CLOSE_COLLECTIONS(record.cget<Vec1Real>(key), expected, tol, key);
    }
    else if constexpr(std::is_same_v<T, Vec1Int>)
    {
        REQUIRE_EQUAL_COLLECTIONS(record.cget<Vec1Int>(key), expected, key);
    }
    else if constexpr(std::is_same_v<T, Vec1String>)
    {
        REQUIRE_EQUAL_COLLECTIONS(record.cget<Vec1String>(key), expected, key);
    }
    else if constexpr(std::is_same_v<T, Vec2Real>)
    {
        REQUIRE_CLOSE_COLLECTIONS_2D(record.cget<Vec2Real>(key), expected, tol, key);
    }
    else if constexpr(std::is_same_v<T, Vec2Int>)
    {
        REQUIRE_EQUAL_COLLECTIONS_2D(record.cget<Vec2Int>(key), expected, key);
    }
    else if constexpr(std::is_same_v<T, Vec2String>)
    {
        REQUIRE_EQUAL_COLLECTIONS_2D(record.cget<Vec2String>(key), expected, key);
    }

    return;
}

inline std::string testName()
{
    using namespace boost::unit_test;

    return std::string(framework::current_test_case().p_name);
}

inline std::string testNameFull()
{
    using namespace boost::unit_test;

    std::string result = framework::current_test_case().p_name;

    auto parent_id = framework::current_test_case().p_parent_id;
    test_suite* parent = nullptr;
    while (parent_id != INV_TEST_UNIT_ID)
    {
        parent = &framework::get<test_suite>(parent_id);
        result = std::string(parent->p_name) + "_" + result;
        parent_id = parent->p_parent_id;
    }

    return result;
}

#endif
