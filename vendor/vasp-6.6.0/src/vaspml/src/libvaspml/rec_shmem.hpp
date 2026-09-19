#ifndef REC_SHMEM_HPP
#define REC_SHMEM_HPP

#include "Record.hpp"     // IWYU pragma: keep
#include "ShmemArray.hpp" // IWYU pragma: keep
#include "types.hpp"      // IWYU pragma: keep

/// Functions dealing with general operations on Records.
namespace vaspml::rec
{

/*******************************************************************************************
 * Copy contents of ShmemArray into Record.
 *
 * @param key Key under which the data will be stored in the Record.
 * @param array ShmemArray to copy data from.
 * @param record Record to copy data into.
 *******************************************************************************************/
template<typename T>
void copyAddShmemArray(const String& key, const ShmemArray<T>& array, Record& record)
{
    record.put<std::vector<T>>(key, std::vector<T>(array.get_size()));
    std::vector<T>& data = record.get<std::vector<T>>(key);
    for (Size i = 0; i < array.get_size(); i++) data[i] = array.get_value(i);
}
/*******************************************************************************************
 * Copy contents of ShmemArray2D into Record.
 *
 * @param key Key under which the data will be stored in the Record.
 * @param array ShmemArray2D to copy data from.
 * @param record Record to copy data into.
 *******************************************************************************************/
template<typename T>
void copyAddShmemArray(const String& key, const ShmemArray2D<T>& array, Record& record)
{
    record.put<std::vector<std::vector<T>>>(key, std::vector<std::vector<T>>());
    std::vector<std::vector<T>>& data = record.get<std::vector<std::vector<T>>>(key);
    data.resize(array.get_size0());
    for (Size i = 0; i < array.get_size0(); i++)
    {
        data[i].resize(array.get_size1());
        for (Size j = 0; j < array.get_size1(); j++) data[i][j] = array.get_value(i, j);
    }
}
/*******************************************************************************************
 * Copy contents of ShmemArray2DVariableLen into Record.
 *
 * @param key Key under which the data will be stored in the Record.
 * @param array ShmemArray2DVariableLen to copy data from.
 * @param record Record to copy data into.
 *******************************************************************************************/
template<typename T>
void copyAddShmemArray(const String& key, const ShmemArray2DVariableLen<T>& array, Record& record)
{
    record.put<std::vector<std::vector<T>>>(key, std::vector<std::vector<T>>());
    std::vector<std::vector<T>>& data = record.get<std::vector<std::vector<T>>>(key);
    data.resize(array.get_size0());
    for (Size i = 0; i < array.get_size0(); i++)
    {
        data[i].resize(array.get_size1(i));
        for (Size j = 0; j < array.get_size1(i); j++) data[i][j] = array.get_value(i, j);
    }
}

} //namespace vaspml::rec

#endif
