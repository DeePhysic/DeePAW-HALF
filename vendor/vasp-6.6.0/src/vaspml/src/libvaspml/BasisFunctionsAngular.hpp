#ifndef BASISFUNCTIONSANGULAR_HPP
#define BASISFUNCTIONSANGULAR_HPP

#include "BasisFunctionsRadialSpline.hpp"
#include "types.hpp"

namespace vaspml
{

enum class BasisFunctionType;
namespace math
{
class SphericalHarmonics;
enum class CutoffType;
}

class BasisFunctionsAngular : public BasisFunctionsRadialSpline
{
  public:
    BasisFunctionsAngular();
    BasisFunctionsAngular(math::CutoffType  cutOffType_in,
                          Real              cutOff_in,
                          Real              widthBroadening_in,
                          Size              nGrid_in,
                          Size              maxOrderBessel_in,
                          Size              nRootsBessel_in,
                          BasisFunctionType type,
                          Size              angularFiltering_in = 9999999,
                          Real              filterScale_in = 0);

    //**********************************************************************************************
    // Regular member functions.
    //**********************************************************************************************
    void computeAngularBasis(const Vec1Real& rhat,
                             const Vec1Real& rnorm,
                             Vec1Real&       ylm,
                             Vec1Real&       ylmd) const;

    void computeAngularBasis(const Real* rhat,
                             const Real* rnorm,
                             const Size& nPairs,
                             Real*       ylm,
                             Real*       ylmd) const;

    void computeAngularBasis(const Real& rhat_x,
                             const Real& rhat_y,
                             const Real& rhat_z,
                             const Real& rnorm,
                             const Size& l,
                             const Int&  m,
                             Real&       ylm,
                             Real&       ylmd_x,
                             Real&       ylmd_y,
                             Real&       ylmd_z) const;

    Size get_ldim() const;

  private:
    void compute_totalNumberBasisFunctions();

    //**********************************************************************************************
    // Getters and setters.
    //**********************************************************************************************
    // Setters
    void set_lmax(Size lmax_in);
    // Getters
    const Vec1Real& get_filteringFactors() const;

    /// Maximum angular quantum number.
    Size lmax;
    /// Number of @f$l,m@f$ combinations for given @f$l_{\max}$, equals size of flattened
    /// SphericalHarmonics arrays.
    Size                                      ldim;
    std::shared_ptr<math::SphericalHarmonics> angularBasis;
    /** Angular filtering type.
     *
     * Angular filtering decreases importance of angular basis with increasing l.
     */
    Size angularFiltering;
    /// Scaling parameter for filtering function.
    Real filterScale;
    /// Computes filteringFactors.
    void computeAngularFilteringCoefficient();
    /// Filtering factors for each angular quantum number.
    Vec1Real filteringFactors;
};

} // namespace vaspml

#endif
