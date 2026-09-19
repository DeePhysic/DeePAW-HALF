#ifndef CUTOFF_HPP
#define CUTOFF_HPP

#include "types.hpp"

namespace vaspml::math
{

enum class CutoffType
{
    /** No cutoff.
     *
     * @f[
     * f_c(r) = 1
     * @f]
     */
    NONE = -1,
    /** Hard cutoff.
     * @f[
     * f_c(r) =
     * \begin{cases}
     * 1, & \text{for } 0 \le r < r_c, \\
     * 0, & \text{else.}
     * \end{cases}
     * @f]
     */
    HARD = 0,
    /** Cutoff used by Behler and Parrinello.
     *
     * @f[
     * f_c(r) = \frac{1}{2} \left[ \cos (\pi \frac{r}{r_c}) + 1\right]
     * @f]
     *
     * See https://doi.org/10.1103/PhysRevLett.98.146401.
     */
    BP = 1,
    /** Cutoff developed by Kazutoshi Miwa.
     *
     * @f[
     * f_c(r) =
     * \begin{cases}
     * 1, & \text{for } 0 \le r < \frac{r_c}{2} \\
     * \frac{1}{4} \left[x(r)^3 - 3x(r) + 2\right], & \text{for } \frac{r_c}{2} \le r < r_c
     * \end{cases}
     * @f]
     * where @f$ x(r) = \frac{4r}{r_c} - 3 @f$.
     */
    MIWA = 2,
    /** Cutoff developed by Ryosuke Jinnouchi (**not implemented**).
     *
     * Not converted to C++ yet.
     */
    JINNOUCHI = 3,
    /** Variation of cutoff used by Willatt, Musil and Ceriotti (**not implemented**).
     *
     * @f[
     * f_c(r) = \frac{c}{c + \left(r / r_0\right)^m}
     * @f]
     *
     * See https://doi.org/10.1039/c8cp05921g, assumes @f$c = 1@f$, @f$m = 7@f$
     * and @f$r_0 = 2\mathrm{\mathring{A}}@f$, not converted to C++ yet.
     */
    WMC = 4,
    /** Cutoff proposed by Singraber, Behler and Dellago.
     *
     * @f[
     * f(x) = ((15 - 6x)x - 10) x^3 + 1, \quad \text{where } x := \frac{r}{r_c}
     * @f]
     *
     * See https://doi.org/10.1021/acs.jctc.8b00770.
     */
    POLY2 = 5,
    /** Exponential bump function, smooth at cutoff radius.
     *
     * @f[
     * f(x) = e^{1 - \frac{1}{1 - \left(r/r_c\right)^2}}
     * @f]
     *
     * See https://en.wikipedia.org/wiki/Bump_function.
     */
    BUMP = 6
};

/** Computes cutoff function values for a given set of radii.
 *
 * @param r Vector with given input radii.
 * @param fc Output vector with cutoff function values.
 * @param rcut Cutoff radius.
 * @param cutoffType Cutoff function type to use, see #CutoffType.
 */
void cutoff(const Vec1Real& r, Vec1Real& fc, Real rcut, CutoffType cutoffType = CutoffType::BP);
/** Computes cutoff function value and derivatives for a given set of radii.
 *
 * @param r Vector with given input radii.
 * @param fc Output vector with cutoff function values.
 * @param dfc Output vector with cutoff function derivatives.
 * @param rcut Cutoff radius.
 * @param cutoffType Cutoff function type to use, see #CutoffType.
 */
void cutoffAndDerivative(const Vec1Real& r,
                         Vec1Real&       fc,
                         Vec1Real&       dfc,
                         Real            rcut,
                         CutoffType      cutoffType = CutoffType::BP);

} // namespace vaspml::math

#endif
