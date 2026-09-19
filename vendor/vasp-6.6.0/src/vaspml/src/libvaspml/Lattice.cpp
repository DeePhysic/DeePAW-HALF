#include "Lattice.hpp"

#include "Record.hpp"
#include "utils.hpp"

#include <cmath>
#include <iostream>

using namespace vaspml;

MapString vaspml::dataMapLattice()
{
    return MapString{
        {"components", "Vec1Real"}
    };
}

Lattice::Lattice(ShRec record) :
    data(assignOrMakeRecord(record, dataMapLattice())),
    components(data->get<Vec1Real>("components"))
{
    volume = 0.0;
}

Lattice::Lattice(Real  ax,
                 Real  ay,
                 Real  az,
                 Real  bx,
                 Real  by,
                 Real  bz,
                 Real  cx,
                 Real  cy,
                 Real  cz,
                 ShRec record) :
    data(assignOrMakeRecord(record, dataMapLattice())),
    components(data->get<Vec1Real>("components"))
{
    set_components(ax, ay, az, bx, by, bz, cx, cy, cz);
}

void Lattice::set_components(Real ax,
                             Real ay,
                             Real az,
                             Real bx,
                             Real by,
                             Real bz,
                             Real cx,
                             Real cy,
                             Real cz)
{
    auto& A = (this->components);
    A.resize(9);
    A[0] = ax;
    A[1] = ay;
    A[2] = az;
    A[3] = bx;
    A[4] = by;
    A[5] = bz;
    A[6] = cx;
    A[7] = cy;
    A[8] = cz;
    computeVolume();
}

void Lattice::set_components(const Vec1Real& componentsVector)
{
    components = componentsVector;
    computeVolume();
}

void Lattice::setComponent(Size i, Size j, Real value)
{
    Size col = i * 3 + j;
    components[col] = value;
}

Real Lattice::component(Size i, Size j) const
{
    Size col = i * 3 + j;
    return components[col];
}

Vec1Real Lattice::timesVector(const Vec1Real& vec) const
{
    const auto& A = this->components;

    const Real x = A[0] * vec[0] + A[3] * vec[1] + A[6] * vec[2];
    const Real y = A[1] * vec[0] + A[4] * vec[1] + A[7] * vec[2];
    const Real z = A[2] * vec[0] + A[5] * vec[1] + A[8] * vec[2];

    return Vec1Real{x, y, z};
}

void Lattice::timesVectorInPlace(Vec1Real& vec) const
{
    const auto& A = this->components;

    const Real x = A[0] * vec[0] + A[3] * vec[1] + A[6] * vec[2];
    const Real y = A[1] * vec[0] + A[4] * vec[1] + A[7] * vec[2];
    const Real z = A[2] * vec[0] + A[5] * vec[1] + A[8] * vec[2];

    vec[0] = x;
    vec[1] = y;
    vec[2] = z;

    return;
}

void Lattice::timesVectorInPlace(Real& v1, Real& v2, Real& v3) const
{
    const auto& A = this->components;

    const Real x1 = A[0] * v1 + A[3] * v2 + A[6] * v3;
    const Real x2 = A[1] * v1 + A[4] * v2 + A[7] * v3;
    const Real x3 = A[2] * v1 + A[5] * v2 + A[8] * v3;

    v1 = x1;
    v2 = x2;
    v3 = x3;

    return;
}

void Lattice::invert()
{
    const auto& A = this->components;

    // clang-format off
    const Real omega = A[0] * (A[4] * A[8] - A[5] * A[7])
                     + A[1] * (A[5] * A[6] - A[3] * A[8])
                     + A[2] * (A[3] * A[7] - A[4] * A[6]);

    const Real xx = (A[4] * A[8] - A[5] * A[7]) / omega;
    const Real xy = (A[7] * A[2] - A[8] * A[1]) / omega;
    const Real xz = (A[1] * A[5] - A[2] * A[4]) / omega;

    const Real yx = (A[5] * A[6] - A[3] * A[8]) / omega;
    const Real yy = (A[8] * A[0] - A[6] * A[2]) / omega;
    const Real yz = (A[2] * A[3] - A[0] * A[5]) / omega;

    const Real zx = (A[3] * A[7] - A[4] * A[6]) / omega;
    const Real zy = (A[6] * A[1] - A[7] * A[0]) / omega;
    const Real zz = (A[0] * A[4] - A[1] * A[3]) / omega;
    // clang-format on
    set_components(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

void Lattice::computeVolume()
{
    const auto& A = this->components;

    // clang-format off
    volume = std::abs(A[0] * (A[4] * A[8] - A[5] * A[7])
                    + A[1] * (A[5] * A[6] - A[3] * A[8])
                    + A[2] * (A[3] * A[7] - A[4] * A[6]));
    // clang-format on

    return;
}

Real Lattice::get_volume() const
{
    return volume;
}

void Lattice::writeToScreen() const
{
    const auto& A = this->components;
    std::cout << A.size() << std::endl;
    std::cout << str("[[%24.16E, %24.16E, %24.16E],\n", A[0], A[1], A[2]);
    std::cout << str(" [%24.16E, %24.16E, %24.16E],\n", A[3], A[4], A[5]);
    std::cout << str(" [%24.16E, %24.16E, %24.16E]]\n", A[6], A[7], A[8]);

    return;
}

//Lattice::Lattice( const Lattice& that ) :
//    data( that.data ),
//    components( that.components ),
//    volume( that.volume )
//{
//}
//
//Lattice& Lattice::operator=( const Lattice& that )
//{
//    if ( this != &that )
//    {
//        this-> data = that.data;
//        this-> components = that.components;
//        this->volume      =  that.volume;
//    }
//    return *this;
//}
//
//Lattice::Lattice( Lattice&& that ) noexcept :
//    data( std::move( that.data ) ), components( that.components ), volume( std::move( that.volume
//    ) )
//
//{
//    that.components.clear();
//}
//
//Lattice& Lattice::operator=( Lattice&& that )
//{
//    if ( this != &that )
//    {
//
//    }
//    return *this;
//}

const Vec1Real& Lattice::get_components() const
{
    return components;
}

void Lattice::convertUnits(const Real& factor)
{
    volume *= (factor * factor * factor);
    for (Real& x : components) x *= factor;
}
