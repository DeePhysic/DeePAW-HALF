#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include "types.hpp"

// Avoid IWYU suggesting <stdlib.h> because of integer argument variant of std::abs (not used here).
// IWYU pragma: no_include <stdlib.h>
#include <algorithm>    // for max, transform, copy
#include <cmath>        // for sqrt, abs
#include <limits>       // for numeric_limits
#include <numeric>      // for inner_product
#include <type_traits>  // for is_floating_point

namespace vaspml
{

/*******************************************************************************************
 * InterpolLinesearch takes function and computes the next x position along the step direction
 * \[f
 * \mathbf{x}_{new} = \mathbf{x}_{old} + \alpha \mathbf{p}
 *\f]
 * where \f$ \alpha \f$ is the scale constant to be determined in this routine
 * and \f$ \mathbf{p} \f$ is the Newton direction along which the function is
 * minimized.
 * This InterpolLinesearch routine satisfies the 1st Wolfe
 * \[f
 *  f(\mathbf{x}_{new}) \leq f(\mathbf{x}_{old}) + c_{1}*\alpha \nabla f
 *(\mathbf{x}_{old})*\mathbf{p}
 *\f]
 *
 * Then a backtracking scheme is used. Therefore the function to be minimized is considered as
 * a function of alpha. In the first step a quadratic function is used to approximate the
 * function to be minimized as
 * \f[
 * \phi( \alpha ) = f( \mathbf{x}_{old} + \alpha \mathbf{p})
 * \f]
 * with gradient
 * \f[
 * \frac{d\phi( \alpha )}{d\alpha} = \nabla f_{old} \mathbf{p}
 * \f]
 * Then the quadratic approximation is written as
 * \f[
 * \phi( \alpha ) \approx [\phi(1)-\phi(0)-\nabla f_{old} \mathbf{p}]\alpha^{2}+\nabla f_{old}
 *\mathbf{p}\alpha+\phi(0)
 * \f]
 * Then the optimal new alpha value can be computed by finding the extremum of this function
 * \f[
 * \alpha = \frac{\nabla f_{old} \mathbf{p}}{2[\phi(1)-\phi(0)-\nabla f_{old} \mathbf{p}]}
 * \f]
 * If this new alpha does not satisfy the first Wolfe condition for sufficient decrease
 * a cubic function can be fitted since more values of \f$ \phi(\alpha) \f$ are available.
 * The cubic function can be written as
 * \f[
 * \phi(\alpha) = a * \alpha^{3}+b*\alpha^{2}+\nabla f_{old} \mathbf{p}\alpha+\phi(0)
 * \f]
 * where the coefficients a and b are given by
 * \f[
 * \begin{bmatrix}
 * a \\
 * b
 * \end{bmatrix} = \frac{1}{\alpha_{1}-\alpha_{2}}
 * \begin{bmatrix}
 *    \frac{1}{\alpha^{2}_{1}} & - \frac{1}{\alpha^{2}_{2}} \\
 *    -\frac{\alpha_{2}}{\alpha^{2}_{1}} & \frac{\alpha_{1}}{\alpha^{2}_{2}}
 * \end{bmatrix}
 * \begin{bmatrix}
 *    \phi(\alpha_{1})-\nabla f_{old} \mathbf{p}\alpha_{1}-\phi(0) \\
 *    \phi(\alpha_{2})-\nabla f_{old} \mathbf{p}\alpha_{2}-\phi(0) \\
 * \end{bmatrix}
 * \f]
 * The new optimized $\alpha$ is now given by
 * \f[
 *    \alpha = \frac{-b+\sqrt{b^{2}-3*a*\nabla f_{old} \mathbf{p}}}{3*a}
 * \f]
 * From then on the \f$ \alpha \f$ is updated again and again until the 1st Wolfe condition is
 * satiesfied.
 * To use this class create an instance
 * ```
 * InterpolLinesearch lnsearch;
 * lnsearch( xStart, fOld, gradient of function, Direction along which to mini, xNew, fNew,
 *maxStep,bool,function object)
 * ```
 *
 * For further information on the Bracketing method the following reference is suggested for
 * mathematical information:
 * Jorge Nocedal Stephen J. Wright Numerical Optimization Second Edition, Chapter 3.5
 * Step-length selection algorithms, Interpolation, p. 57
 * For information on the here implemented algorithm see
 * Numerical Recipes The Art of Scientific Computing Third Edition,
 * Chapter Line Searches and Bracketing, p. 479
 *******************************************************************************************/
template<typename T>
class InterpolLinesearch
{
    static_assert(std::is_floating_point<T>::value, "Wrong argument type in InterpolLinesearch\n");

  public:
    InterpolLinesearch();
    /*******************************************************************************************
     * Make a line search step with the Newton method. The method checks if the function
     * decreases sufficiently. If not the step length will be adapted.
     *
     * @param xOld x values at previous step
     * @param fOld old function value evaluated at xOld
     * @param grad gradient of function
     * @param newtDir newton direction obtained by J^{-1}F where is jacobian and F is function value
     * @param x vector which will be computed during the step
     * @param f function value which will be computed during step
     * @param maxStep is the maximal step width allowed to the algorithm
     * @param check variable to see if convergence was reached.
     * @param func class of the function to be minimized. Call has to look like f(x);
     *******************************************************************************************/
    template<typename U>
    void operator()(const std::vector<T>& xOld,
                    const T&              fOld,
                    const std::vector<T>& grad,
                    std::vector<T>&       newtDir,
                    std::vector<T>&       x,
                    T&                    f,
                    const T&              maxStep,
                    bool&                 check,
                    U&                    func);

  private:
    /*******************************************************************************************
     * parameter alpha to control step width
     *******************************************************************************************/
    const T c1 = (T)1.0e-4;
    /*******************************************************************************************
     * epsilon to control mini step width.
     *******************************************************************************************/
    const T TOLX = std::numeric_limits<T>::epsilon();
};

template<typename T>
InterpolLinesearch<T>::InterpolLinesearch()
{}

template<typename T>
template<typename U>
void InterpolLinesearch<T>::operator()(const std::vector<T>& xOld,
                                       const T&              fOld,
                                       const std::vector<T>& grad,
                                       std::vector<T>&       newtDir,
                                       std::vector<T>&       x,
                                       T&                    f,
                                       const T&              maxStep,
                                       bool&                 check,
                                       U&                    func)
{
    check = false;
    T sum = std::sqrt(std::inner_product(newtDir.cbegin(), newtDir.cend(), newtDir.cbegin(), (T)0));
    // rescale step when too large
    if (sum > maxStep)
        std::transform(newtDir.begin(),
                       newtDir.end(),
                       newtDir.begin(),
                       [&maxStep, &sum](T& value) { return value *= maxStep / sum; });
    T slope = std::inner_product(grad.cbegin(), grad.cend(), newtDir.cbegin(), (T)0);
    // this check is not so nice for the LBFGS algorithm since it often prevents the algorithm from
    // converging
    //if ( slope > 0 )
    //{
    //    std::cout << slope << std::endl;
    //    throw std::runtime_error("Problem in InterpolLinesearch<T>. Newton direction not along
    //    descent direction");
    //}

    //T test = std::transform_reduce( newtDir.cbegin(), newtDir.cend(), xOld.cbegin(),
    //                                (T)0,
    //                                []( const T& a, const T& b ) { return std::max(a, b); },
    //                                []( const T& dir, const T& xi ) {
    //                                    return std::abs( dir ) / std::max( std::abs( xi ), 1.0 );
    //                                });
    T test = (T)0;
    for (Size i = 0; i < newtDir.size(); ++i)
    {
        test = std::max(test, std::abs(newtDir[i]) / std::max(std::abs(xOld[i]), 1.0));
    }
    T alphaMin = TOLX / test;
    T alpha1 = (T)1;
    T alpha2 = (T)0;
    T f2 = (T)0;
    while (true)
    {
        // update x;
        std::transform(xOld.cbegin(),
                       xOld.cend(),
                       newtDir.cbegin(),
                       x.begin(),
                       [&alpha1](const T& xold_i, const T& dir) { return xold_i + alpha1 * dir; });
        f = func(x);
        // the newly determined alpha is less than the minimal alpha. Return starting x value
        if (alpha1 < alphaMin)
        {
            std::copy(xOld.cbegin(), xOld.cend(), x.begin());
            check = true;
            return;
        }
        // First Wolfe condition for sufficient decrease is satiesfied
        // return new x value and new function value. Everything went well.
        else if ((f <= fOld + c1 * alpha1 * slope)) { return; }
        else
        {
            T tmpAlpha;
            if (alpha1 == (T)1)
            {
                // in first round compute the new alpha with the quadtratic approx
                tmpAlpha = -slope / ((T)2 * (f - fOld - slope));
            }
            else
            {
                // do the alpha update with the cubic function approximation
                T rhs1 = f - fOld - alpha1 * slope;
                T rhs2 = f2 - fOld - alpha2 * slope;
                // compute coefficients for the cubic alpha computation
                T a = (rhs1 / (alpha1 * alpha1) - rhs2 / (alpha2 * alpha2)) / (alpha1 - alpha2);
                T b = (-alpha2 * rhs1 / (alpha1 * alpha1) + alpha1 * rhs2 / (alpha2 * alpha2))
                    / (alpha1 - alpha2);
                // check if a is zero because then the denominator can't be used.
                if (a == (T)0) tmpAlpha = -slope / ((T)2 * b);
                else
                {
                    T sqrtTerm = b * b - (T)3 * a * slope;
                    // check if root can be computed
                    if (sqrtTerm < (T)0) tmpAlpha = (T)0.5 * alpha1;
                    else if (b <= (T)0) tmpAlpha = (-b + std::sqrt(sqrtTerm)) / ((T)3 * a);
                    else tmpAlpha = -slope / (b + std::sqrt(sqrtTerm));
                }
                // if alpha gets too large reduce alpha1 by factor 2 instead of taking computed
                // alpha
                if (tmpAlpha > (T)0.5 * alpha1) tmpAlpha = (T)0.5 * alpha1;
            }
            alpha2 = alpha1;
            f2 = f;
            alpha1 = std::max(tmpAlpha, (T)0.1 * alpha1);
        }
    }
}

/*******************************************************************************************
 * Ths routine performs an Armijo-Wolfe line search. During the line search the paramater
 * \f$ alpha \f$ of the following equation is determined

 * \[f
 * \mathbf{x}_{new} = \mathbf{x}_{old} + \alpha \mathbf{p}.
 *\f]
 * The parameter \f$ \alpha \f$ is optimized until the Armijo condition is satisfied. The
 * Armijo condition is given by
 *\[f
 * f( x_{k} + \alpha_{k}p_{k}) \le f(x_{k}) +c_{1}+\alpha_{k}\nabla f^{T}_{k}p_{k},
 *\f]
 * where \f$ p_{k} \f$ is the Newton Direction determined by
 *\[f
 * p_{k} = -(\nabla^{2}f_{k})^{-1}\nabla f_{k}.
 *\f]
 * Additionally the second Wolfe contion has to be satisfied after determining
 * \f$ \alpha \f$. The second Wolfe condition is given by
 *\[f
 * \nabla f(x_{k}+\alpha_{k}p_{k})^{T}p_{k} \ge c_{2}\nabla f(x_{k})^{T}p_{k}.
 *\f]
 * The algorithm uses a bracketing approach where an upper \f$ \alpha_{max} \f$ and
 * lower \f$ \alpha_{min} \f$ bound to the optimal \f$ \alpha \f$ is defined. These
 * bracketing values are then adapted until the optimal \f$ \alpha \f$ is found.
 * If the Armijo condition is not satiesfied
 * \f$ \alpha_{max} \f$ is set to the actual \f$ \alpha \f$ because the function starts
 * increasing again if going further. Then the new \f$ \alpha \f$ is set to the average of
 * \f$ \alpha_{max} \f$ and \f$ \alpha_{min} \f$. If the Wolfe condition is not satiesfied
 * the new lower bound is set to the actual \f$ \alpha \f$ value. And then the algorithm
 * computes the new \f$ \alpha \f$ again as the average between of upper and lower bound.
 * Then a second phase of optimization is done in which \f$ \alpha \f$ is reduced until the
 * Armijo condition is satiesfied. In this phase a criterion for noisy gradients could
 * be used. For further information check
 * A noise-tolerant quasi-newton algorithm for unconstrained optimization H. J. M. Shi
 * Y. Xie, R. Byrd, Jorge Nocedal arxiv 09. September 2021
*******************************************************************************************/
template<typename T>
class ArmijoWolfeLinesearch
{
  public:
    ArmijoWolfeLinesearch(const T& c1 = 1e-4, const T& c2 = 0.9);
    /*******************************************************************************************
     * Function to update the x value along the Newton direction
     *
     * @param xOld current x value from which linesearch is started
     * @param fOld function value at xOld
     * @param grad gradient of function at xOld
     * @param newtDir Newton direction along which linesearch is done
     * @param x value after update is done along Newton direction
     * @param f function value at x.
     * @param maxStep maximal allowed step
     * @param check control variable to check for convergence
     * @param func function which has to be optimized. Function has to implement () operator and
     * df( x, grad ) member function.
     *******************************************************************************************/
    template<typename U>
    void operator()(const std::vector<T>& xOld,
                    const T&              fOld,
                    const std::vector<T>& grad,
                    std::vector<T>&       newtDir,
                    std::vector<T>&       x,
                    T&                    f,
                    const T&              maxStep,
                    bool&                 check,
                    U&                    func);

  private:
    /*******************************************************************************************
     * Reduce alpha unitl second Wolfe contion is satiesfied. Here a noise condition could be
     * added if the gradient would contain noise.
     *
     * @param xOld current x value from which linesearch is started
     * @param fOld function value at xOld
     * @param grad gradient of function at xOld
     * @param newtDir Newton direction along which linesearch is done
     * @param x value after update is done along Newton direction
     * @param f function value at x.
     * @param maxStep maximal allowed step
     * @param check control variable to check for convergence
     * @param func function which has to be optimized. Function has to implement () operator and
     * df( x, grad ) member function.
     *******************************************************************************************/
    template<typename U>
    void splitPhase(const std::vector<T>& xOld,
                    const T&              fOld,
                    const std::vector<T>& grad,
                    std::vector<T>&       newtDir,
                    std::vector<T>&       x,
                    T&                    f,
                    const T&              maxStep,
                    bool&                 check,
                    U&                    func);
    /*******************************************************************************************
     * Constant for determining the Armijo condtion. Typical value is 1E-4
     *******************************************************************************************/
    T c1;
    /*******************************************************************************************
     * Constant for the second Wolfe condition. Typical value is 0.9
     *******************************************************************************************/
    T c2;
    /*******************************************************************************************
     * Variable which determines how far along the Newton direction the current x
     * value will be changed
     *******************************************************************************************/
    T alpha;
    /*******************************************************************************************
     * Maximal number of iterations
     *******************************************************************************************/
    Size maxIter = 20;
    /*******************************************************************************************
     * Variable stores Armijo condition at input x value xOld.
     *******************************************************************************************/
    T armijoCheck;
};

template<typename T>
ArmijoWolfeLinesearch<T>::ArmijoWolfeLinesearch(const T& c1, const T& c2) : c1(c1), c2(c2)
{}

template<typename T>
template<typename U>
void ArmijoWolfeLinesearch<T>::operator()(const std::vector<T>& xOld,
                                          const T&              fOld,
                                          const std::vector<T>& grad,
                                          std::vector<T>&       newtDir,
                                          std::vector<T>&       x,
                                          T&                    f,
                                          const T&              /* maxStep */,
                                          bool&                 check,
                                          U&                    func)
{
    T lower = 0;
    T upper = std::numeric_limits<T>::max();
    alpha = (T)1;
    std::vector<T> xi = xOld;
    f = fOld;

    T gp = std::inner_product(grad.cbegin(), grad.cend(), newtDir.cbegin(), (T)0);
    armijoCheck = fOld + c1 * alpha * gp;
    T wolfeCheck = c2 * gp;

    std::vector<T> gradNew = grad;

    for (Size iter = 0; iter < maxIter; iter++)
    {
        // update x with x = xold + alpha * newtonDir;
        std::transform(xOld.cbegin(),
                       xOld.cend(),
                       newtDir.cbegin(),
                       x.begin(),
                       [&](const T& xoi, const T& dir) { return xoi + alpha * dir; });
        //std::cout << "X values " << x[0] << "   " << x[1] << std::endl;
        f = func(x);
        //        func.df( x, gradNew );
        if (f > armijoCheck) // check if Armijo condition is violated.
        {
            // do backtracking
            upper = alpha;
            alpha = (T)0.5 * (upper + lower);
        }
        else if (func.df(x, gradNew);
                 std::inner_product(gradNew.cbegin(), gradNew.cend(), newtDir.cend(), (T)0)
                 < wolfeCheck)
        //check Wolfe condition
        {
            //std::cout << std::inner_product( gradNew.cbegin(), gradNew.cend(), newtDir.cend(),
            //(T)0 ) << "   "
            //          << wolfeCheck << std::endl;
            lower = alpha;
            if (upper == std::numeric_limits<T>::max()) { alpha = (T)2 * alpha; }
            else { alpha = (T)0.5 * (upper + lower); }
        }
        else // everything went fine. The algo succeded
        {
            check = true;
            return;
        }
    }
    //std::cout << "x end "<< x[0] << "   "<< x[1]<< std::endl;
    //throw std::runtime_error( "Too many iterations in ArmijoWolfeLinesearch<T>" );
}

template<typename T>
template<typename U>
void ArmijoWolfeLinesearch<T>::splitPhase(const std::vector<T>& xOld,
                                          const T&              /* fOld */,
                                          const std::vector<T>& /* grad */,
                                          std::vector<T>&       newtDir,
                                          std::vector<T>&       x,
                                          T&                    f,
                                          const T&              /* maxStep */,
                                          bool&                 /* check */,
                                          U&                    func)
{
    while (f > armijoCheck)
    {
        alpha = alpha / 10;
        // update x with x = xold + alpha * newtonDir;
        std::transform(xOld.cbegin(),
                       xOld.cend(),
                       newtDir.cbegin(),
                       x.begin(),
                       [&](const T& xoi, const T& dir) { return xoi + alpha * dir; });
        f = func(x);
    }
}

} //namespace vaspml

#endif
