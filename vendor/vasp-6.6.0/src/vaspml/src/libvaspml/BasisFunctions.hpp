#ifndef BASISFUNCTIONS_HPP
#define BASISFUNCTIONS_HPP

#include "types.hpp"

namespace vaspml
{

// Forward-declare cutoff type enum.
namespace math
{
enum class CutoffType;
}

enum class BasisFunctionType
{
    none = 0,
    bodyOrder2 = 1,
    bodyOrder3 = 2
};

/*******************************************************************************************
 * @class BasisFunctions
 * @brief Virtual base class from which any basis function class should inherit.
 *
 * This class is a purely virtual base class for basis functions. This class implements the
 * basic functionality which any basis function class will need. Class implements for
 * exmaple a cutoff radius, the number of used basis functions with angular paramter (l)
 * dependence. The virtual member functions are implemented as pure to ensure that all
 * derived classes inherit and implement these functions as a minimum level of functionality.
 *******************************************************************************************/
class BasisFunctions
{
  public:
    /// Initialize empty class with all parameters set to zero.
    BasisFunctions();
    /*******************************************************************************************
     * Construct class with desired parameters.
     *
     * @param cutoffType_in Desired cutoff type, see choices in math::CutoffType.
     * @param cutoff_in Desired cutoff radius.
     * @param width_broadening_in Width for Gaussian broadening of atom positions.
     *******************************************************************************************/
    BasisFunctions(const math::CutoffType cutoffType_in,
                   const Real             cutoff_in,
                   const Real             widthBroadening_in,
                   BasisFunctionType      type);

    //**********************************************************************************************
    // Regular member functions.
    //**********************************************************************************************
    void printValues();

    //**********************************************************************************************
    // Getters and setters (non-virtual).
    //**********************************************************************************************
    // Setters
    void set_cutoffType(const math::CutoffType in);
    void set_cutoff(const Real in);
    void set_widthBroadening(const Real in);
    // Getters
    math::CutoffType get_cutoffType() const;
    Real             get_cutoff() const;
    Real             get_widthBroadening() const;

    //**********************************************************************************************
    // Virtual functions to be implemented by derived types.
    //
    // Also check derived types for additional documentation.
    //**********************************************************************************************
    /// Get maximum spherical harmonic angular number (same as maximum spherical Bessel function
    /// order).
    Size get_maxOrder() const;
    /// Get grid spacing of radial grid.
    Real get_gridSpacing() const;
    /// Get number of grid points.
    Int get_nGrid() const;
    /// Get number of spherical Bessel function roots.
    Int get_nRootsBessel() const;
    /// Get the total number of basis functions.
    Size get_totalNumberBasisFunctions() const;
    Size get_totalNumberBesselFunctions() const;
    /// Get the total number of @f$l, m@f$ combinations.
    Size get_ldim() const;
    /*******************************************************************************************
     * Get number of valid roots found for radial basis.
     *
     * @returns Vector with number of valid roots. Index of vector determines the radial basis
     * function order.
     *******************************************************************************************/
    const Vec1Size& get_nValidRoots() const;
    /*******************************************************************************************
     * Set the number of grid points for radial grid.
     *
     *
     * @param in Number of radial grid points.
     *******************************************************************************************/
    void set_nGrid(const Size in);
    /*******************************************************************************************
     * Set maximum order (lmax) for radial functions.
     *
     * @param in Value to set for lmax.
     *******************************************************************************************/
    void set_maxOrderBessel(const Size in);
    /*******************************************************************************************
     * Set maximum number of roots in radial functions.
     *
     * @param in Value for desired maximum number of roots of radial basis functions.
     *******************************************************************************************/
    void set_nRootsBessel(const Size in);
    /// Print parameters of radial basis function to screen.
    void printValuesRadial() const;
    void update();
    /*******************************************************************************************
     * Compute interpolated radial basis functions for a given x-value.
     *
     * @param x-value at which functions are given back.
     * @return Radial basis function values at x.
     *******************************************************************************************/
    Vec2Real interpolate(const Real x) const;
    /*******************************************************************************************
     * Compute interpolated radial basis functions for a given x-value.
     *
     * @param x-value at which functions are given back.
     * @param result Radial basis function values at x (see note).
     *
     * @note Memory allocation of result vector is performed inside the function. Results are stored
     * as two-dimensional vector:
     * - index 1: Bessel function order @f$l@f$.
     * - index 2: Root number @f$n@f$.
     *******************************************************************************************/
    void interpolate(const Real x, Vec2Real& result) const;
    /*******************************************************************************************
     * Compute interpolated radial basis functions and derivative for a given x-value.
     *
     * @param x-value at which functions are given back.
     * @param result Radial basis function values at x (see note).
     * @param result_derivative Radial basis function derivatives at x (see note).
     *
     * @note Memory allocation of result vector is performed inside the function. Results are stored
     * as two-dimensional vector:
     * - index 1: Bessel function order @f$l@f$.
     * - index 2: Root number @f$n@f$.
     *******************************************************************************************/
    void interpolate(const Real x, Vec2Real& result, Vec2Real& result_derivative) const;
    /*******************************************************************************************
     * Compute interpolated radial basis functions and derivative for a given x-value (flattened
     * result vectors).
     *
     * @param x-value at which functions are given back.
     * @param result Radial basis function values at x (see note and warning).
     * @param result_derivative Radial basis function derivatives at x (see note and warning).
     *
     * @note Results are stored as flattened one-dimensional vector with the following memory
     * layout:
     * - Root number @f$n@f$ (fastest)
     * - Bessel function order @f$l@f$ (slowest)
     *
     * @warning The memory of return variables has to be managed outside of this function, the size
     * must be correct on entry. A size check is performed in debugging levels 1+.
     *******************************************************************************************/
    void interpolate(const Real x, Vec1Real& result, Vec1Real& result_derivative) const;
    void interpolate(const Real  x,
                     const Size& offest,
                     Vec1Real&   result,
                     Vec1Real&   result_derivative) const;
    void interpolate(const Real  x,
                     const Size& order,
                     const Size& root,
                     Real&       result,
                     Real&       result_derivative) const;
    /*******************************************************************************************
     * Computing angular part of basis functions.
     *
     * @param rhat Flattened vector of normalized 3D position vectors (see note).
     * @param rnorm Vector with norm of all position vectors (same order as position vectors).
     * @param ylm Flattened array with angular part of basis functions.
     * @param ylmd Flattened array with angular part gradient of basis functions.
     *
     * @note This function requires a vector of normalized position vectors @f$\hat{\vec{r}}@f$
     * (@f$|\hat{\vec{r}}| = 1@f$) as input. The multi-dimensional data must be passed in a
     * flattened vector with this memory layout from fastest to slowest dimension:
     *   - Cartesian coordinate x, y, z (fastest)
     *   - Position vector index (slowest)
     *
     * For output vector format look at the documentation of
     * math::SphericalHarmonics::computeSphericalHarmonicsAndGradients().
     *******************************************************************************************/
    void computeAngularBasis(const Vec1Real& rhat,
                             const Vec1Real& rnorm,
                             Vec1Real&       ylm,
                             Vec1Real&       ylmd) const;
    /*******************************************************************************************
     * Computing angular part of basis functions. Input parameters in raw pointers.
     *
     * @param rhat Flattened vector of normalized 3D position vectors (see note).
     * @param rnorm Vector with norm of all position vectors (same order as position vectors).
     * @param nPair parameter which defines the number of pairs between distance is given. Length
     * of rnorm input
     * @param ylm Flattened array with angular part of basis functions.
     * @param ylmd Flattened array with angular part gradient of basis functions.
     *
     * @note This function requires a vector of normalized position vectors @f$\hat{\vec{r}}@f$
     * (@f$|\hat{\vec{r}}| = 1@f$) as input. The multi-dimensional data must be passed in a
     * flattened vector with this memory layout from fastest to slowest dimension:
     *   - Cartesian coordinate x, y, z (fastest)
     *   - Position vector index (slowest)
     *
     * For output vector format look at the documentation of
     * math::SphericalHarmonics::computeSphericalHarmonicsAndGradients().
     *******************************************************************************************/
    void computeAngularBasis(const Real* rhat,
                             const Real* rnorm,
                             const Size  nPair,
                             Real*       ylm,
                             Real*       ylmd) const;
    /*******************************************************************************************
     * Computing angular part of basis functions.
     *
     * @param rhat_x x component of normalized input vector
     * @param rhat_y y component of normalized input vector
     * @param rhat_z z component of normalized input vector
     * @param rnorm norm of input vector
     * @param l angular quantum number for which angular basis is returned
     * @param m magnetic quantum number for which angular basis is returned
     * @param ylm component for given input parameters.
     * @param ylmd Flattened array with angular part gradient of basis functions.
     *
     * @note This function requires a vector of normalized position vectors @f$\hat{\vec{r}}@f$
     * (@f$|\hat{\vec{r}}| = 1@f$) as input. The multi-dimensional data must be passed in a
     * flattened vector with this memory layout from fastest to slowest dimension:
     *   - Cartesian coordinate x, y, z (fastest)
     *   - Position vector index (slowest)
     *
     * For output vector format look at the documentation of
     * math::SphericalHarmonics::computeSphericalHarmonicsAndGradients().
     *******************************************************************************************/
    void computeAngularBasis(const Real& rhat_x,
                             const Real& rhat_y,
                             const Real& rhat_z,
                             const Real& rnorm,
                             const Size& l,
                             const Size& m,
                             Real&       ylm,
                             Real&       ylmd_x,
                             Real&       ylmd_y,
                             Real&       ylmd_z) const;

  private:
    /// Cutoff function type for the basis functions.
    math::CutoffType cutoffType;
    /// Cutoff radius for the basis functions.
    Real cutoff;
    /// Width for Gaussian broadening of atom positions (variance).
    Real widthBroadening;

    BasisFunctionType functionType = BasisFunctionType::none;

    // List of friend classes. It is important to have some derived classes with access to private
    // members.  However, not all derived classes should have this privilege.
    friend class BasisFunctionsRadialSpline;
    friend class BasisFunctionsAngular;
};

using BasisFunctionMap = std::map<String, std::shared_ptr<BasisFunctions>>;

} // namespace vaspml

#endif
