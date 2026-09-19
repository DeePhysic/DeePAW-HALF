#include "ArrayResizing.hpp"

using namespace vaspml;

ArrayResizing1D::ArrayResizing1D() : actDim(0), maxDim(0)
{}

ArrayResizing1D::ArrayResizing1D(const Size& n) : actDim(n), maxDim(n)
{}

void ArrayResizing1D::init(const Size& n)
{
    actDim = n;
    maxDim = n;
}

ArrayResizing2D::ArrayResizing2D() : act1Dim(0), max1Dim(0), actSizeRec(0), maxSizeRec(0)
{}

ArrayResizing2D::ArrayResizing2D(const Size& n) : act1Dim(n), max1Dim(n), actSize(n), maxSize(n)
{}
