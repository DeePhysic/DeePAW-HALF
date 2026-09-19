#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "types.hpp"

#include <iostream>
#include <stdexcept>

namespace vaspml
{

template<class T>
class Matrix
{
  public:
    Matrix(const Int n0, const Int n1, const std::shared_ptr<std::vector<T>>& data = nullptr);

    std::shared_ptr<std::vector<T>> get_data();

    Int  get_firstDimension() const;
    Int  get_secondDimension() const;
    void print();

  private:
    std::shared_ptr<std::vector<T>> data;
    Int                             dim0;
    Int                             dim1;
};

template<class T>
Matrix<T>::Matrix(const Int n0, const Int n1, const std::shared_ptr<std::vector<T>>& data)
{
    dim0 = n0;
    dim1 = n1;

    if (data == nullptr) { this->data = std::make_shared<std::vector<T>>(); }
    else
    {
        if ((Size)dim0 * dim1 != data->size())
        {
            throw std::runtime_error(
                "ERROR: Matrix<T>::Matrix( const Int n0, const Int n1, const ShVec1Real& data ) \n"
                "supplied dimensions n0 and n1 don't match length of supplied vector");
        }
        this->data = data;
    }

    std::cout << "Class created " << std::endl;
}

template<class T>
std::shared_ptr<std::vector<T>> Matrix<T>::get_data()
{
    return data;
}

template<class T>
Int Matrix<T>::get_firstDimension() const
{
    return dim0;
}

template<class T>
Int Matrix<T>::get_secondDimension() const
{
    return dim1;
}

template<class T>
void Matrix<T>::print()
{
    Size col = 0;
    for (Size i = 0; i < (Size)dim0; i++)
    {
        for (Size j = 0; j < (Size)dim1; j++)
        {
            std::cout << (*data)[col] << "  ";
            col++;
        }
        std::cout << std::endl;
    }
}

} //namespace vaspml

#endif
