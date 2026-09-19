#ifndef BASISFUNCTIONSRADIALSPLINE_HPP
#define BASISFUNCTIONSRADIALSPLINE_HPP

#include "BasisFunctions.hpp"
#include "cubic_spline.hpp"
#include "types.hpp"

namespace vaspml::math
{
enum class CutoffType;
}

namespace vaspml
{

/*******************************************************************************************
 * @class BasisFunctionsRadialSpline
 * @brief Implements radial basis function for description of pure 2 body interactions
 *
 * Class is used to express to body descriptors which only depend on the distance between
 * two objects. The distances between the 2 bodies will be smeared by a gaussian and then
 * projected ontto modified spherical Bessel functions. For a detailed descriptions
 * of the mathematics behind this two body basis functions please read the paper
 *
 * <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.014105">
 * On-the-fly machine learning force field generation: Application to melting points
 * Ryosuke Jinnouchi, Ferenc Karsai, and Georg Kresse
 * Phys. Rev. B 100, 014105 – Published 17 July 2019</a>
 *******************************************************************************************/
class BasisFunctionsRadialSpline : public BasisFunctions
{
  public:
    /// Intializes empty class and all parameters are set to zero.
    BasisFunctionsRadialSpline();
    /** Constructor with initialization, calls update(), i.e., spline coefficients etc. are
     * initialized.
     *
     * @param cutoffType_in Selects cutoff type, see math::CutoffType.
     * @param cutoff_in Selects the cutoff radius.
     * @param width_broadening_in Width for Gaussian broadening of atom positions.
     * @param nGrid_in Number of grid points used for radial grid.
     * @param maxOrderBessel_in Maximum order of the Bessel functions used.
     * @param nRootsBessel_in Maximum number of roots used per Bessel function.
     */
    BasisFunctionsRadialSpline(math::CutoffType  cutoffType_in,
                               Real              cutoff_in,
                               Real              widthBroadening_in,
                               Size              nGrid_in,
                               Size              maxOrderBessel_in,
                               Size              nRootsBessel_in,
                               BasisFunctionType type);

    //**********************************************************************************************
    // Regular member functions.
    //**********************************************************************************************
    void printValuesRadial() const;
    /** Updating the spline coefficients, has to be called when parameters changed by one/several of
     * the setters.
     */
    void     update();
    Vec2Real interpolate(const Real x) const;
    void     interpolate(const Real x, Vec2Real& values) const;
    void     interpolate(const Real x, Vec2Real& values, Vec2Real& derivatives) const;
    void     interpolate(const Real x, Vec1Real& values, Vec1Real& derivatives) const;
    void     interpolate(const Real  x,
                         const Size& offset,
                         Vec1Real&   values,
                         Vec1Real&   derivatives) const;
    void     interpolate(const Real& x,
                         const Size& order,
                         const Size& root,
                         Real&       values,
                         Real&       derivatives) const;
    void     computeAngularBasis(const Vec1Real& rhat,
                                 const Vec1Real& rnorm,
                                 Vec1Real&       ylm,
                                 Vec1Real&       ylmd) const;
    void     computeAngularBasis(const Real* rhat,
                                 const Real* rnorm,
                                 const Size  nPair,
                                 Real*       ylm,
                                 Real*       ylmd) const;
    void     computeAngularBasis(const Real& rhat_x,
                                 const Real& rhat_y,
                                 const Real& rhat_z,
                                 const Real& rnorm,
                                 const Size& l,
                                 const Size& m,
                                 Real&       ylm,
                                 Real&       ylmd_x,
                                 Real&       ylmd_y,
                                 Real&       ylmd_z) const;

    //**********************************************************************************************
    // Getters and setters.
    //**********************************************************************************************
    // Getters
    Size            get_maxOrder() const;
    Real            get_gridSpacing() const;
    Int             get_nGrid() const;
    Int             get_nRootsBessel() const;
    Size            get_totalNumberBasisFunctions() const;
    Size            get_totalNumberBesselFunctions() const;
    Size            get_ldim() const;
    const Vec1Size& get_nValidRoots() const;
    // Setters
    void set_nGrid(const Size in);
    void set_maxOrderBessel(const Size in);
    void set_nRootsBessel(const Size in);

  private:
    /** Maximum order of spherical Bessel functions used as radial basis functions.
     *
     *  Corresponds to maximum second index @f$l_{\max}@f$ of basis function @f$\chi_{nl}@f$.
     */
    Size maxOrderBessel_store;
    /// Maximum order of spherical Bessel functions used as radial basis functions + 1.
    Size maxOrderBessel;
    /** Maximum number of spherical Bessel function zeros used as radial basis functions.
     *
     *  Corresponds to maximum first index @f$n_{\max}@f$ of basis function @f$\chi_{nl}@f$.
     */
    Size nRootsBessel;
    /** Number of valid spherical Bessel function roots found by math::sphericalBesselRoots.
     *
     * Vector contains number of valid roots for each spherical Bessel function order.
     */
    Vec1Size nValidRoots;
    /// Cutoff radius / number points.
    Real gridSpacing;
    /// Total number of basis functions in radial descriptor.
    Size totalNumberBasisFunctions;
    Size totalNumberBesselFunctions;
    /// Number of grid points for spline interpolation.
    Size nGrid;
    /// Minimum of largest root found for each spherical Bessel function order, used for rescaling.
    Real rootMax = 1e12;
    /// Grid points for spline interpolation.
    Vec1Real grid;
    /// Cutoff function values computed on grid within cutoff sphere.
    Vec1Real cutoffFunction;
    /// Cutoff function derivatives computed on grid within cutoff sphere.
    Vec1Real cutoffFunction_derivative;
    /** Cubic spline interpolated versions of radial basis functions (spline coefficients).
     *
     * @note This array will only contain the proper spline coefficients when the function update()
     * finished successfully.
     */
    std::vector<std::vector<CubicSpline<Real>>> splines;

    //**********************************************************************************************
    // Private member functions.
    //**********************************************************************************************
    /// Set up grid for spline interpolation.
    void makeGrid();
    /// Compute cutoff function and derivatives on grid points.
    void computeCutoffFunction();
    /// Compute grid spacing when cutoff and nGrid are set.
    void set_gridSpacing();
    /** Compute the roots of a spherical Bessel basis function.
     *
     * @param rootsBesselBasis On output contains the roots of the spherical Bessel function.
     * @param actOrderBessel Order of the spherical Bessel function for which roots are computed.
     *
     * The first nRootsBessel are computed.
     */
    void computeSphericalBesselRoots(Vec2Real& rootsBesselBasis, const Int actOrderBessel);
    /** Check if root finding algorithm converged.
     *
     * @param rootsBesselBasis On output contains the roots of the spherical Bessel function with
     * invalid roots set to zero.
     *
     * If the algorithm did not converge for certain roots then those are set to zero. Count number
     * of valid roots and stores it in totalNumberBasisFunctions.
     */
    void checkSphericalBesselRootsConvergence(Vec2Real& rootsBesselBasis);
    /** Rescale roots by cutoff radius.
     *
     * @param rootsBesselBasis On output contains the roots of the spherical Bessel function scaled
     * via the cutoff radius.
     *
     * This guarantees that basis functions pass trough zero at the cutoff radius. rootMax is also
     * scaled.
     */
    void rescaleRoots(Vec2Real& rootsBesselBasis);
    /** Compute the spherical Bessel functions for all orders, all roots and all grid points.
     *
     * @param bessel 3D vector containing the spherical Bessel function values (see note).
     * @param rootsBesselBasis Precomuted, rescaled roots of the spherical Bessel functions.
     *
     * @note Memory allocation is performed inside this function (index 1). Results are stored as
     * three-dimensional vector:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     * - index 3: Grid point.
     */
    void computeSphericalBesselFunctions(Vec3Real& bessel, const Vec2Real& rootsBesselBasis);
    /** Compute the spherical Bessel functions for one order and all corresponding roots and grid
     * points.
     *
     * @param bessel 2D vector containing the spherical Bessel function values (see note).
     * @param roots First nBasis roots.
     * @param actOrderBessel Desired order of the spherical Bessel functions.
     * @param nBasis Number of spherical Bessel functions with given order.
     *
     * @note Memory allocation is performed inside this function (index 2). Results are stored as
     * two-dimensional vector:
     * - index 1: Spherical Bessel function root @f$n@f$.
     * - index 2: Grid point.
     */
    void computeSphericalBesselFunctionsSingleOrder(Vec2Real&       bessel,
                                                    const Vec1Real& roots,
                                                    const Size      actOrderBessel,
                                                    const Size      nBasis);
    /** Compute spherical Bessel function basis set for a single root and all grid points.
     *
     * @param bessel Vector containing the spherical Bessel function values for all grid points.
     * @param root Root which is used to rescale grid points.
     * @param actOrderBessel Order of spherical Bessel function.
     */
    void computeSphericalBesselFunctionsSingleRoot(Vec1Real& bessel,
                                                   Real      root,
                                                   const Int actOrderBessel);
    /** Normalize precomputed spherical Bessel function values for all orders and roots.
     *
     * @param bessel 3D vector containing spherical Bessel function values which are normalized on
     * output (see note).
     *
     * @note Results are stored as three-dimensional vector:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     * - index 3: Grid point.
     *
     * After calling this function radial basis functions are normalized with respect to the first
     * @f$n@f$ index only, i.e.,
     *
     * @f[
     * \chi_{nl} \chi_{n'l} = \delta_{nn'}.
     * @f]
     */
    void computeBesselNorms(Vec3Real& bessel);
    /** Normalize precomputed spherical Bessel function values for one order.
     *
     * @param bessel 2D vector containing spherical Bessel function values which are normalized on
     * output (see note).
     * @param nBasis Number of spherical Bessel functions for given order.
     *
     * @note Results are stored as two-dimensional vector:
     * - index 1: Spherical Bessel function root @f$n@f$.
     * - index 2: Grid point.
     *
     * See also computeBesselNorms().
     */
    void computeBesselNormsSingleOrder(Vec2Real& bessel, Size nBasis);

    /** Computes the radial part of the projection function onto the basis functions.
     *
     * This corresponds to the function @f$h_{nl}@f$ found in expression (19) in Phys. Rev. B 100,
     * 014105 (see also https://arxiv.org/pdf/1904.12961.pdf).
     *
     * @param sphericalBessel Precomputed spherical Bessel function values, normalized and root
     * scaled.
     * @param projectionRadial Contains values of the function @f$h_{nl}(r)@f$.
     *
     * This function computes the following integral:
     * @f[
     * h_{nl}\left(r\right) = \frac{4\pi}{\left(\sqrt{2\sigma_\mathrm{atom}^{2}\pi}\right)^{3}}
     * f_\mathrm{cut}\left(r\right)
     * \int  _{0}^{\infty} \chi_{nl}\left(r'\right)
     * \mathrm{exp}\left(-\frac{r^{2}+r'^{2}}{2\sigma_\mathrm{atom}^{2}}\right)
     * \iota_{l}\left(\frac{rr'}{\sigma_{\mathrm  {atom}}^{2}}\right)
     * r'^{2}dr'
     * @f]
     *
     * @note Input and output are stored as three-dimensional vector, output is allocated inside
     * this function:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     * - index 3: Grid point.
     */
    void computeProjectionRadial(const Vec3Real& sphericalBessel, Vec3Real& projectionRadial);
    /** Similar to computeProjectionRadial(), computes the derivative of the radial part of the
     * projection function at the grid origin.
     *
     * This corresponds to the derivative of the function @f$h_{nl}(r)@f$ at the grid origin found
     * in expression (19) in Phys. Rev. B 100, 014105 (see also
     * https://arxiv.org/pdf/1904.12961.pdf). The derivative is needed to get the correct boundary
     * behavior of the spline interpolation.
     *
     * @param sphericalBessel Precomputed spherical Bessel function values, normalized and root
     * scaled.
     * @param projectionRadialDerivative Contains values of the derivative
     * @f$\frac{d h_{nl}(r)}{d r}{\big|}_{r=r_0}@f$, where @f$r_0@f$ is the grid origin.
     *
     * This function computes the following expression:
     * @f[
     * \begin{split}
     * \frac{dh_{nl}}{dr}{\Big|}_{r=r_0} &= \frac{4\pi}{\left(2\sigma^{2}_{atom}\right)^{3/2}}
     * \left\{ \left[\frac{\partial f_{cut}(r))}{\partial r}{\Big|}_{r=r_0} -
     * \frac{r_0}{\sigma_{atom}^{2}}f_{cut}(r_0) \right] \right.
     * \int_{0}^{\infty}\chi_{nl}(r')e^{-\frac{{r'}^{2}+r_0^{2}}{2\sigma^{2}_{atom}}}
     * \iota_{l}(\frac{r_0r'}{\sigma^{2}}){r'}^{2}dr' \\
     * &\left.+ \int_{0}^{\infty}\chi_{nl}(r')e^{-\frac{{r'}^{2}+r_0^{2}}{2\sigma^{2}_{atom}}}
     * \frac{\partial \iota_{l}(\frac{rr'}{\sigma^{2}})}{\partial r}{\Big|}_{r=r_0}
     * \frac{{r'}^{3}}{\sigma^{2}_{atom}} dr'\right\}
     * \end{split}
     * @f]
     *
     * @note Input is a three-dimensional vector:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     * - index 3: Grid point.
     * @note Two-dimensional output vector is allocated inside this function:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     */
    void computeProjectionRadialDerivativeOrigin(const Vec3Real& sphericalBessel,
                                                 Vec2Real&       projectionRadialDerivative);
    /** Helper function for computeProjectionRadialDerivativeOrigin(), computes derivative for a
     * single Bessel function order and root.
     *
     * For details see formula in description of computeProjectionRadialDerivativeOrigin().
     *
     * @param factor1 Prefactor for term 1.
     * @param factor2 Prefactor for term 2.
     * @param factor2 Prefactor for term 3.
     * @param gauss Gaussian function values on the grid.
     * @param sphericalBessel Spherical Bessel function values on the grid.
     * @param modifiedBessel Modified spherical Bessel function values on the grid.
     * @param modifiedBesselDerivative Modified spherical Bessel function derivatives on the grid.
     */
    Real computeProjectionRadialDerivativeOriginSingle(
        const Real               factor1,
        const Real               factor2,
        const Real               factor3,
        const Vec1Real&          gauss,
        const Vec1Real&          sphericalBessel,
        Vec1Real::const_iterator modifiedBessel,
        Vec1Real::const_iterator modifiedBesselDerivative);

    /** Compute cubic spline coefficients for a given function.
     *
     * Here, this is used to create a spline representation for radial part of the projection
     * function @f$h_{nl}(r)@f$, see computeProjectionRadial();
     *
     * @param functions The functions to be splined (see note).
     * @param derivative The derivative of the function to be splined at the first (left-most) grid
     * point.
     *
     * @note The function input must be a three-dimensional vector:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     * - index 3: Grid point.
     * @note The derivative at the first grid point must be a two-dimensional vector:
     * - index 1: Spherical Bessel function order @f$l@f$.
     * - index 2: Spherical Bessel function root @f$n@f$.
     */
    void computeSplineCoefficients(const Vec3Real& functions, const Vec2Real& derivative);

    friend class BasisFunctionsAngular;
};

} // namespace vaspml

#endif
