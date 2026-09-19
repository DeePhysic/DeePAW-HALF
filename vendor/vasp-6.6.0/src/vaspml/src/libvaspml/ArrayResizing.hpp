#ifndef ARRAYRESIZING_HPP
#define ARRAYRESIZING_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <type_traits>

namespace vaspml
{

struct ArrayResizing1D
{
    ArrayResizing1D();
    ArrayResizing1D(const Size& n);
    void init(const Size& n);
    bool inline checkResize(const Size& newDim);

    Size actDim;
    Size maxDim;
};

bool inline ArrayResizing1D::checkResize(const Size& newDim)
{
    bool result = false;
    if (newDim > maxDim)
    {
        maxDim = newDim;
        result = true;
    }
    actDim = newDim;
    return result;
}

/*******************************************************************************************
 * ArrayResizing2D is a helping structure to check if certain arrays should be resized.
 *
 * This function is meant to be used in combination with the GPU port of vaspml.
 * with this function the user is capable of only resizing the manipulated data arrays
 * if really needed. This minimizes copying back and forth large data arrays from GPU to
 * CPU.
 *
 * Example:
 * @code
 * Vec2Real data( n );
 * ArrayResizing2D resizer( n );
 * newSizes = neighborList.get_numberNeighbors();
 * if ( resizer.resizeNecessary( newSizes ) ){
 *    for ( Size i = 0; i < n; i++ )
 *    {
 *        data[i].resize( resizer.maxSize[i] );
 *    }
 * }
 * @endcode
 *******************************************************************************************/
struct ArrayResizing2D
{
    /*******************************************************************************************
     * default constructor, zero initialization
     *******************************************************************************************/
    ArrayResizing2D();
    /*******************************************************************************************
     * constructor which initializes size arrays with given length and values zero
     *
     * @param[ in ] length of first dimension of array that has to be treated
     *******************************************************************************************/
    ArrayResizing2D(const Size& n);

    template<typename T>
    void init(const Size& dim1, const std::vector<T>& sizes, const Size& multi = 1);
    /*******************************************************************************************
     * Function to decide if  an array has to be resized. Which is important for GPU
     *
     * The variables maxSize and actSize will be set to new values in this function.
     * Function returns a bool telling the user if a resize has to be done or not.
     * True is returned if resize is necessary, and false if there is no resize needed
     *
     * @param[in] newSize actual size needed for the next computation step
     *******************************************************************************************/
    template<typename T>
    bool inline checkResize2Dim(const std::vector<T>& newSize, const Size& multi = 1);
    bool inline checkResize2Dim(const Size& newSize, const Size& multi = 1);
    bool inline checkResize1Dim(const Size& newSize);
    template<typename T>
    bool inline refill(const std::vector<T>& newSize, const Size& multi = 1);
    // void inline checkResize( const Vec1Size& newSize );
    /*******************************************************************************************
     * resize the data array which is connected to the ArrayResizing2D instance
     *
     * @param[inout] data 2dimensional std vector which will be resized to maxSize
     *******************************************************************************************/
    template<typename T>
    void resizeArray(std::vector<std::vector<T>>& data);
    template<typename T, typename U>
    void resizeArray1Dim(std::vector<std::vector<T>>& data,
                         const std::vector<U>&        newSize,
                         const Size&                  multi = 1);
    template<typename T>
    void resizeArray1Dim(std::vector<std::vector<T>>& data,
                         const Size&                  newSize,
                         const Size&                  multi = 1);
    template<typename T, typename U>
    void resizeArray2Dim(std::vector<std::vector<T>>& data,
                         const std::vector<U>&        newSize,
                         const Size&                  multi = 1);
    template<typename T>
    void resizeArray2Dim(std::vector<std::vector<T>>& data,
                         const Size&                  newSize,
                         const Size&                  multi = 1);

    /*******************************************************************************************
     * leadning dimension of the 2d dimensional array which size will be controlled
     *******************************************************************************************/
    Size act1Dim;
    Size max1Dim;
    /*******************************************************************************************
     * actual size of the second array dimension ( not maximal capacity )
     *
     * this array will be set with the function resizeNecessary
     *******************************************************************************************/
    Vec1Size actSize;
    /*******************************************************************************************
     * maximal size of the second array dimension, denotes the maximal capacity (can be increased)
     *
     * this array will be set with the function resizeNecessary
     *******************************************************************************************/
    Vec1Size maxSize;
    Size     actSizeRec;
    Size     maxSizeRec;
};

bool inline ArrayResizing2D::checkResize1Dim(const Size& newSize)
{
    bool result = false;
    if (newSize > max1Dim)
    {
        max1Dim = newSize;
        result = true;
    }
    act1Dim = newSize;
    return result;
}

bool inline ArrayResizing2D::checkResize2Dim(const Size& newSize, const Size& multi)
{
    bool result = false;
    if (newSize * multi > maxSizeRec)
    {
        maxSizeRec = newSize * multi;
        result = true;
    }
    actSizeRec = newSize * multi;
    return result;
}

template<typename T>
bool inline ArrayResizing2D::checkResize2Dim(const std::vector<T>& newSize, const Size& multi)
{
    VASPML_DEBUG_L1(
        String funcName = __func__;
        if (!std::is_integral<T>::value)
        {
            global::tutor.bug("ERROR in" + funcName + " ! \nSupplied data type is not integral");
        }
    );
    bool result = false;
    for (Size i = 0; i < actSize.size(); i++)
    {
        if (newSize[i] * multi > maxSize[i])
        {
            maxSize[i] = newSize[i] * multi;
            result = true;
        }
        actSize[i] = newSize[i];
    }
    return result;
}

template<typename T>
bool inline ArrayResizing2D::refill(const std::vector<T>& newSize, const Size& multi)
{
    VASPML_DEBUG_L1(
        String funcName = __func__;
        if (!std::is_integral<T>::value)
        {
            global::tutor.bug("ERROR in" + funcName + " ! \nSupplied data type is not integral");
        }
    );
    bool result = false;
    for (Size i = 0; i < actSize.size(); i++)
    {
        if (newSize[i] * multi != maxSize[i])
        {
            maxSize[i] = newSize[i] * multi;
            result = true;
        }
        actSize[i] = newSize[i];
    }
    return result;
}

template<typename T>
void ArrayResizing2D::resizeArray(std::vector<std::vector<T>>& data)
{
    for (Size i = 0; i < act1Dim; i++) { data[i].resize(maxSize[i]); }
}

template<typename T>
void ArrayResizing2D::init(const Size& dim1, const std::vector<T>& sizes, const Size& multi)
{
    // For integral types only:
    VASPML_DEBUG_L1(
        String funcName = __func__;
        if (!std::is_integral<T>::value)
        {
            global::tutor.bug("ERROR in" + funcName + " ! \nSupplied data type is not integral");
        }
    );
    act1Dim = dim1;
    max1Dim = dim1;
    actSize.resize(dim1);
    maxSize.resize(dim1);
    for (Size i = 0; i < dim1; i++)
    {
        actSize[i] = sizes[i] * multi;
        maxSize[i] = sizes[i] * multi;
    }
}

template<typename T, typename U>
void ArrayResizing2D::resizeArray1Dim(std::vector<std::vector<T>>& data,
                                      const std::vector<U>&        newSize,
                                      const Size&                  multi)
{
    VASPML_DEBUG_L1(
        String funcName = __func__;
        if (!std::is_integral<U>::value)
        {
            global::tutor.bug("ERROR in" + funcName
                              + " ! \nSupplied data type of newSize is not integral");
        }
    );

    actSize.resize(max1Dim);
    maxSize.resize(max1Dim);
    data.resize(max1Dim);
    resizeArray2Dim(data, newSize, multi);
}

template<typename T>
void ArrayResizing2D::resizeArray1Dim(std::vector<std::vector<T>>& data,
                                      const Size&                  newSize,
                                      const Size&                  multi)
{
    actSize.resize(max1Dim);
    maxSize.resize(max1Dim);
    data.resize(max1Dim);
    resizeArray2Dim(data, newSize, multi);
}

template<typename T, typename U>
void ArrayResizing2D::resizeArray2Dim(std::vector<std::vector<T>>& data,
                                      const std::vector<U>&        newSize,
                                      const Size&                  multi)
{
    VASPML_DEBUG_L1(
        String funcName = __func__;
        if (!std::is_integral<U>::value)
        {
            global::tutor.bug("ERROR in" + funcName
                              + " ! \nSupplied data type of newSize is not integral");
        }
    );

    for (Size i = 0; i < max1Dim; i++)
    {
        if (newSize[i] * multi > maxSize[i]) maxSize[i] = newSize[i] * multi;
        actSize[i] = newSize[i] * multi;
        data[i].resize(maxSize[i]);
    }
}

template<typename T>
void ArrayResizing2D::resizeArray2Dim(std::vector<std::vector<T>>& data,
                                      const Size&                  newSize,
                                      const Size&                  multi)
{
    if (newSize * multi > maxSizeRec)
    {
        maxSizeRec = newSize * multi;
        for (Size i = 0; i < max1Dim; i++) { data[i].resize(maxSizeRec); }
    }
    actSizeRec = newSize * multi;
}

} // namespace vaspml

#endif // ARRAYRESIZING_HPP
