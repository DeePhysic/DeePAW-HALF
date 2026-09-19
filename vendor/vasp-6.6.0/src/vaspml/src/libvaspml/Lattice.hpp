#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "types.hpp"

namespace vaspml
{
/***************************************************************************************************
 * Manage crystal lattices, use to rescale atoms from direct to Cartesian coordinates.
 *
 * This class implements a crystal lattice. It can be used to rescale atom positions from
 * Cartesian to direct coordinates and vice versa. The class can create its own inverse
 * which can then be used to invert the rescaling of the Direct coordinates to Cartesian.
 * The memory of the class can be managed from the outside by supplying a `shared_ptr` to
 * the constructor. If omitted, the memory can be managed by the class instance itself.
 *
 * Usage example:
 * @code
 * Lattice lattice;
 * // Set the elements with
 * lattice.set_components(ax, ay, az, bx, by, bz, cx, cy, cz);
 * // where a, b and c are the individual lattice vectors and x, y, z denote their respective
 * // components.
 * // Alternatively, provide a vector containing the lattice vector components:
 * lattice.set_components(vectorFormatData);
 * lattice.computeVolume();
 * Lattice inverseLattice;
 * inverseLattice.set_components( lattice.get_components() );
 * inverseLattice.invert();
 * @endcode
 **************************************************************************************************/
class Lattice
{
  public:
    /***********************************************************************************************
     * Creating Lattice where memory can be supplied.
     *
     * @param record Memory for the lattice components, 9 numbers required.
     *
     * @note When `record` is supplied the memory of the class can be managed and modified from
     * the outside. If `record` is not supplied memory will be managed by class itself.
     **********************************************************************************************/
    Lattice(ShRec record = nullptr);
    /***********************************************************************************************
     * Construct Lattice class with supplying elements directly in the constructor.
     *
     * @param record components Memory for the lattice components, 9 numbers required.
     * @param ax x-component of the first lattice vector a.
     * @param ay y-component of the first lattice vector a.
     * @param az z-component of the first lattice vector a.
     * @param bx x-component of the second lattice vector b.
     * @param by y-component of the second lattice vector b.
     * @param bz z-component of the second lattice vector b.
     * @param cx x-component of the third lattice vector c.
     * @param cy y-component of the third lattice vector c.
     * @param cz z-component of the third lattice vector c.
     *
     * @note When `record` is supplied the memory of the class can be managed and modified from
     * the outside. If `record` is not supplied memory will be managed by class itself.
     **********************************************************************************************/
    Lattice(Real  ax,
            Real  ay,
            Real  az,
            Real  bx,
            Real  by,
            Real  bz,
            Real  cx,
            Real  cy,
            Real  cz,
            ShRec record = nullptr);
    /***********************************************************************************************
     * Setting the components of the lattice by supplying 9 numbers.
     *
     * @param ax x-component of the first lattice vector a.
     * @param ay y-component of the first lattice vector a.
     * @param az z-component of the first lattice vector a.
     * @param bx x-component of the second lattice vector b.
     * @param by y-component of the second lattice vector b.
     * @param bz z-component of the second lattice vector b.
     * @param cx x-component of the third lattice vector c.
     * @param cy y-component of the third lattice vector c.
     * @param cz z-component of the third lattice vector c.
     **********************************************************************************************/
    void
    set_components(Real ax, Real ay, Real az, Real bx, Real by, Real bz, Real cx, Real cy, Real cz);
    // clang-format on
    /***********************************************************************************************
     * Set components of lattice by supplying a flat vector.
     *
     * @param components 9 element vector with lattice components in this order: `ax, ay, az, bx,
     *                   by, bz, cx, cy, cz`.
     **********************************************************************************************/
    void            set_components(const Vec1Real& componentsVector);
    const Vec1Real& get_components() const;
    /***********************************************************************************************
     * Set individual component of single lattice vector.
     *
     * @param i lattice vector index, 0, 1, 2 corresponds to a, b, c.
     * @param j lattice vector component, 0, 1, 2 corresponds to x, y, z.
     * @param value Numerical value of component.
     *
     * @note Whenever individual components are changed, the stored volume will not be automatically
     * recomputed. To update the volume, call computeVolume() after compoment change.
     **********************************************************************************************/
    void setComponent(Size i, Size j, Real value);
    /***********************************************************************************************
     * Get individual component of single lattice vector.
     *
     * @param i lattice vector index, 0, 1, 2 corresponds to a, b, c.
     * @param j lattice vector component, 0, 1, 2 corresponds to x, y, z.
     *
     * @return Numerical value of component.
     *
     **********************************************************************************************/
    Real component(Size i, Size j) const;
    /***********************************************************************************************
     * Multiply lattice matrix with given input vector.
     *
     * Typically used to transfer direct input coordinates @f$\vec{d}@f$ to Cartesian output
     * coordinates @f$\vec{r}@f$:
     *
     * @f[
     * \vec{r} = \begin{pmatrix}
     * \vec{a}^T \\
     * \vec{b}^T \\
     * \vec{c}^T
     * \end{pmatrix} \cdot \vec{d} =
     * \begin{pmatrix}
     * a_x & a_y & a_z \\
     * b_x & b_y & b_z \\
     * c_x & c_y & c_z
     * \end{pmatrix} \cdot
     * \begin{pmatrix} d_1 \\ d_2 \\ d_3 \end{pmatrix}
     * @f]
     *
     * Alternatively, if the inverse lattice is stored, e.g. in `inv` here:
     *
     * @code
     * Lattice inv = lattice.computeInverse();
     * @endcode
     *
     * then the direct coordinates can be recovered from the Cartesian coordinates:
     *
     * @f[
     * \vec{d} = \begin{pmatrix}
     * \vec{a}^T \\
     * \vec{b}^T \\
     * \vec{c}^T
     * \end{pmatrix}^{-1} \cdot \vec{r}
     * @f]
     *
     * @param vec Input vector to be multiplied with lattice.
     * @return Output vector containing product of lattice times input vector.
     **********************************************************************************************/
    Vec1Real timesVector(const Vec1Real& vec) const;
    /***********************************************************************************************
     * Multiply lattice matrix with given input vector (in-place version).
     *
     * See timesVector() for further details.
     *
     * @param vec Input vector to be multiplied with lattice, will be modified, contains the result
     *            on exit.
     **********************************************************************************************/
    void timesVectorInPlace(Vec1Real& vec) const;
    /***********************************************************************************************
     * Multiply lattice matrix with given input vector (in-place, component input version).
     *
     * See timesVector() for further details.
     *
     * @param v1 First component of input vector, contains result on exit.
     * @param v2 Second component of input vector, contains result on exit.
     * @param v3 Third component of input vector, contains result on exit.
     **********************************************************************************************/
    void timesVectorInPlace(Real& v1, Real& v2, Real& v3) const;
    /***********************************************************************************************
     * Invert the lattice matrix.
     *
     * If the lattice vectors @f$a, b, c@f$ are column vectors this function computes:
     *
     * @f[
     * \begin{pmatrix}
     * \vec{a}^T \\
     * \vec{b}^T \\
     * \vec{c}^T
     * \end{pmatrix}^{-1}
     * @f]
     *
     * The product of the lattice matrix and the inverse will result in the unit matrix.
     **********************************************************************************************/
    void invert();
    /***********************************************************************************************
     * Compute the volume of the box spanned by the lattice vectors.
     *
     * Computed volume will be stored to the member variable volume and can be retrieved by
     * get_volume().
     **********************************************************************************************/
    void computeVolume();
    /// Get precomputed value of volume spanned by this lattice.
    Real get_volume() const;
    /// Write lattice to screen.
    void writeToScreen() const;
    //Lattice ( const Lattice& that );
    //Lattice& operator=( const Lattice& that );
    //Lattice ( Lattice&& that ) noexcept;
    //Lattice& operator=( Lattice&& that );

    void convertUnits(const Real& factor);

  private:
    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /***********************************************************************************************
     * Memory for the lattice vector components in the following order: `ax, ay, az, bx, by, bz, cx,
     * cy, cz`.
     *
     * @note Can be managed from outside if provided in constructor.
     **********************************************************************************************/
    Vec1Real& components;
    /***********************************************************************************************
     * Volume of simulation box, will be recomputed if all components are reset at the same time.
     **********************************************************************************************/
    Real volume;
};

/**********************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **********************************************************************************************/
MapString dataMapLattice();

} // namespace vaspml

#endif
