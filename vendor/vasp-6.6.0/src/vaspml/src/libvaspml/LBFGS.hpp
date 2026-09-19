#ifndef LBFGS_HPP
#define LBFGS_HPP

#include "Linesearch.hpp"
#include "types.hpp"

// Avoid IWYU suggesting <stdlib.h> because of integer argument variant of std::abs (not used here).
// IWYU pragma: no_include <stdlib.h>
#include <algorithm>   // for transform, max, copy, for_each
#include <cmath>       // for abs, sqrt
#include <limits>      // for numeric_limits
#include <numeric>     // for inner_product
#include <stdexcept>   // for runtime_error
#include <type_traits> // for is_floating_point

namespace vaspml
{

template<typename T>
inline T sign(const T& value)
{
    return (T(0) < value) - (value < T(0));
}

/*******************************************************************************************
 *  Minimize the supplied function with a quasi-Newton algorithm in LBFGS form.
 *
 * The objective function will be approximated by a convex function of the form
 *\f[
 * m_{k}(\mathbf{p})=&f(\mathbf{x}_{k})+\nabla \mathbf{f}^{T}(\mathbf{x}_{k})\mathbf{p}+
 *        \frac{1}{2}\mathbf{p}^{T}\nabla^{2}f(\mathbf{x}_{k})\mathbf{p}\\
 *                =&f_{k}+\nabla f^{T}_{k}\mathbf{p} +
 *\frac{1}{2}\mathbf{p}^{T}\mathbf{B}_{k}\mathbf{p}.
 *\f]
 * The optimization of the function is done along the Newton direction \f$ \mathbf{p} \f$, given by
 *\f[
 * \mathbf{p}_{k}=-\mathbf{B}^{-1}_{k}\nabla f_{k}
 *\f]
 * The new x position is then found by optimizing the paramater \f$ \alpha \f$ defining how long
 * the step will be along Newton's direction
 *\f[
 * \mathbf{x}_{k+1}=\mathbf{x}_{k}+\alpha \mathbf{p}_{k}.
 *\f]
 * The step width can be found by any linesearch algorithm. The linesearch algorithm
 * has to satisfy Wolfe conditions too keep the approximate Hessian positive definite.
 * The update of the approximate Hessian is achieved by the following optimization
 * criterion
 *\f[
 * \min_{\mathbf{B}} ||\mathbf{B}-\mathbf{B}_{k} || \\
 * \text{ subject to } \mathbf{B}=\mathbf{B}^{T} \text{ , } \mathbf{B}\mathbf{s}_{k}=\mathbf{y}_{k}
 *\f]
 * In the LBFGS the approximate to the inverse Hessian is computed by a two loop
 * recursive form. For this recursive form curvature pairs
 * \f$ \{\mathbf{s}_{k},\mathbf{y}_{k}\} \f$. The number of curvature pairs is an input
 * parameter supplied by the user. Given the curvature pairs the inverse Hessian can be
 * computed by the following equation
 *\f[
 * For easier notation the quantity \f$ \mathbf{V}_{k} \f$ has to be defined. \f$ \mathbf{V}_{k} \f$
 *is a rank 1 matrix and given by
 * \f[
 * \mathbf{V}_{k} = I - \rho_{k}\mathbf{y}_{k}\mathbf{s}^{T}_{k}.
 * \f]
 * With this the iterative scheme to compute $H_{k}$ is given by
 *\f[
 * \mathbf{H}_{k}&=(\mathbf{V}_{k-1}^{T}\cdots
 *\mathbf{V}^{T}_{k-m})\mathbf{H}_{k}^{0}(\mathbf{V}_{k-m}\cdots
 * \mathbf{V}_{k-1})\\
 *            &+\rho_{k-m}(\mathbf{V}_{k-1}^{T}\cdots
 *\mathbf{V}^{T}_{k-m+1})\mathbf{s}_{k-m}\mathbf{s}^{T}_{k-m}
 *                        (\mathbf{V}_{k-m+1}\cdots \mathbf{V}_{k-1})\\
 *            &+\rho_{k-m+1}(\mathbf{V}_{k-1}^{T}\cdots
 *\mathbf{V}^{T}_{k-m+2})\mathbf{s}_{k-m+1}\mathbf{s}^{T}_{k-m+1}
 *                        (\mathbf{V}_{k-m+2}\cdots \mathbf{V}_{k-1})\\
 *            &+\cdots \\
 *            &+\rho_{k-1}\mathbf{s}_{k-1}\mathbf{s}^{T}_{k-1}.
 * Note that this recursive scheme represents a sum of rank-1 matrices. This recursive
 * equation is not used directly but in a form such that the product of the approximate
 * inverse hessian times the gradient is readily obtained.
 *******************************************************************************************/
template<typename T>
class LBFGS
{
    static_assert(std::is_floating_point<T>::value, "Wrong argument type in LBFGS\n");

  public:
    /*******************************************************************************************
     * Find minimum of supplied function func.
     *
     * @param xPos[inout] on intput this is the initial gues for the minimum on output the optimized
     *value
     * @param gradTol[in] accepted tolerance in the gradient.
     * @param nIter[out]  stores number of iterations needed until BFGS converges
     * @param fReturned[out] function value on output
     * @param func function to be optimized. Function has to implement func( xPos ) to obtain
     * function value and func.df( xPos, grad ) for gradint computation at xPos.
     *******************************************************************************************/
    template<typename U>
    void minimize(const Int&      m,
                  std::vector<T>& xPos,
                  const T&        gradTol,
                  Size&           nIter,
                  T&              fReturned,
                  U&              func);

  private:
    /*******************************************************************************************
     * Initialize Local arrays needed for minimize function. Initializes Hessian matrix as
     * unit matrix.
     *******************************************************************************************/
    void init(const Size& dim, const Size& m);
    /*******************************************************************************************
     * Compute new Newton directon with the two loop recursion.
     *******************************************************************************************/
    void compute_newtDirNew(const Size& iter);
    /*******************************************************************************************
     * Function which collects the curvature pairs sk,yk.
     *
     * If number of iterations is higher than sizeHistory the oldest pair is removed.
     * All the other pairs are moved to the back once. The new pair is put into the front of the
     *list.
     *******************************************************************************************/
    void collectCurvaturePairs(const Size&           iter,
                               const std::vector<T>& deltaGrad,
                               const std::vector<T>& deltaX);

    /*******************************************************************************************
     * Maximal number of iterations to find function minimum
     *******************************************************************************************/
    const Size MAXITER = 100;
    /*******************************************************************************************
     * Smallest numeric value of T which can be distinguished from zero. The factor ten
     * is a trial and error thing and should be considered as a heuristic to
     * improve convergence rate
     *******************************************************************************************/
    const T EPS = std::numeric_limits<T>::epsilon() * 10;
    /*******************************************************************************************
     * Check if change in x is converged.
     *******************************************************************************************/
    const T TOLX = (T)4 * EPS;
    /*******************************************************************************************
     * Scaling constant to determine maximal step width
     *******************************************************************************************/
    const T STPMX = (T)100;
    /*******************************************************************************************
     * Difference between actual function gradient and previous function gradient
     *******************************************************************************************/
    std::vector<T> deltaGrad;
    /*******************************************************************************************
     * Stores gradient obtained during previous iteration.
     *******************************************************************************************/
    std::vector<T> gradOld;
    /*******************************************************************************************
     * Stores gradient obtained during current iteration.
     *******************************************************************************************/
    std::vector<T> grad;
    /*******************************************************************************************
     * Stores new x-value obtained during current iteration.
     *******************************************************************************************/
    std::vector<T> xNew;
    /*******************************************************************************************
     * stores difference between xOld and xNew. Step which was done on the independent variable.
     *******************************************************************************************/
    std::vector<T> deltaX;
    /*******************************************************************************************
     * Storing history of \f$ \rho = \frac{1}{s_{k}^{T}y_{k}} \f$
     *******************************************************************************************/
    std::vector<T> rhos;
    /*******************************************************************************************
     * Stores the history of the differences of x-values in consecutive iterations.
     * Dimension of the array is sizeHistory x nDim.
     *******************************************************************************************/
    std::vector<std::vector<T>> sk;
    /*******************************************************************************************
     * Stores the history of the differences of the gradients in consecutive iterations.
     * Dimension of the array is sizeHistory x nDim.
     *******************************************************************************************/
    std::vector<std::vector<T>> yk;
    /*******************************************************************************************
    *******************************************************************************************/
    std::vector<T> scales;
    /*******************************************************************************************
     *Dimension of minimization problem
     *******************************************************************************************/
    Size nDim;
    /*******************************************************************************************
     * Determines the number of curvature pairs which are used for computing the approximate
     * Hessian with the two loop recursion.
     *******************************************************************************************/
    Size sizeHistory;
    /*******************************************************************************************
     * Newton direction which is computed with two loop recursion.
     *******************************************************************************************/
    std::vector<T> q;
    //InterpolLinesearch<T> linesearch;
    ArmijoWolfeLinesearch<T> linesearch;
};

template<typename T>
void LBFGS<T>::init(const Size& dim, const Size& m)
{
    if (dim != nDim)
    {
        nDim = dim;
        sizeHistory = m;
        deltaGrad.resize(nDim);
        grad.resize(nDim);
        xNew.resize(nDim);
        gradOld.resize(nDim);
        deltaX.resize(nDim);
        // allocate curvature pairs
        sk.resize(m);
        rhos.resize(m);
        std::for_each(sk.begin(), sk.end(), [&](Vec1Real& s) { s.resize(nDim); });
        yk.resize(m);
        std::for_each(yk.begin(), yk.end(), [&](Vec1Real& y) { y.resize(nDim); });
        q.resize(nDim);
        scales.resize(m);
    }
}

template<typename T>
void LBFGS<T>::compute_newtDirNew(const Size& iter)
{
    std::transform(grad.cbegin(), grad.cend(), q.begin(), [](const T& gi) { return gi; });
    Int nSamples = iter < sizeHistory ? iter : sizeHistory;
    for (Int i = 0; i < nSamples; i++)
    {
        T dotp = std::inner_product(yk[i].cbegin(), yk[i].cend(), sk[i].cbegin(), (T)0);
        rhos[i] = T(1) / ((std::abs(dotp) > EPS) ? dotp : sign(dotp) * EPS);
        scales[i] = rhos[i] * std::inner_product(sk[i].cbegin(), sk[i].cend(), q.cbegin(), (T)0);
        std::transform(yk[i].cbegin(),
                       yk[i].cend(),
                       q.cbegin(),
                       q.begin(),
                       [s = scales[i]](const T& yki, const T& qi) { return qi - s * yki; });
    }
    // multiply by initial hessian matrix
    if (nSamples > 1)
    {
        T gamma = std::inner_product(sk[1].cbegin(), sk[1].cend(), yk[1].cend(), (T)0);
        T denom = std::inner_product(yk[1].cbegin(), yk[1].cend(), yk[1].begin(), (T)0);
        gamma /= (std::abs(denom) > EPS) ? denom : sign(denom) * EPS;
        std::transform(q.cbegin(),
                       q.cend(),
                       q.begin(),
                       [gamma](const T& qi) { return gamma * qi; });
    }

    for (Int i = nSamples - 1; i >= 0; i--)
    {
        T beta = rhos[i] * std::inner_product(yk[i].cbegin(), yk[i].cend(), q.cbegin(), (T)0);
        //q = q - ( scales[i] - beta ) * sk[i];
        std::transform(sk[i].cbegin(),
                       sk[i].cend(),
                       q.cbegin(),
                       q.begin(),
                       [s = scales[i] - beta](const T& ski, const T& qi) { return qi + s * ski; });
    }
    std::transform(q.cbegin(), q.cend(), q.begin(), [](const T& qi) { return -qi; });
}

template<typename T>
void LBFGS<T>::collectCurvaturePairs(const Size&           iter,
                                     const std::vector<T>& deltaGrad,
                                     const std::vector<T>& deltaX)
{
    if (iter < sizeHistory)
    {
        yk[iter] = deltaGrad;
        sk[iter] = deltaX;
    }
    else
    {
        for (Size i = sizeHistory - 1; i > 0; i--)
        {
            std::copy(yk[i - 1].begin(), yk[i - 1].end(), yk[i].begin());
            std::copy(sk[i - 1].begin(), sk[i - 1].end(), sk[i].begin());
        }
        yk[0] = deltaGrad;
    }
}

template<typename T>
template<typename U>
void LBFGS<T>::minimize(const Int&      m,
                        std::vector<T>& xPos,
                        const T&        gradTol,
                        Size&           nIter,
                        T&              fReturned,
                        U&              func)
{
    // compute pk new with compute_newtDirNew
    // make line search
    // collect curvature pairs
    // do until convergence
    init(xPos.size(), m);
    Real fp = func(xPos);
    func.df(xPos, grad);

    T    sum = std::inner_product(xPos.cbegin(), xPos.cend(), xPos.cbegin(), (T)0);
    T    stpMax = STPMX * std::max(std::sqrt(sum), (T)nDim);
    bool check = false;
    for (Size iter = 0; iter < MAXITER; iter++)
    {
        nIter = iter;
        compute_newtDirNew(iter);
        linesearch(xPos, fp, grad, q, xNew, fReturned, stpMax, check, func);
        // compute deltaX
        std::transform(xPos.cbegin(),
                       xPos.cend(),
                       xNew.cbegin(),
                       deltaX.begin(),
                       [](const T& xold, const T& xnew) { return xnew - xold; });
        T qTip = std::inner_product(deltaX.cbegin(), deltaX.cend(), deltaX.cend(), (T)0);
        //update x
        xPos = xNew;
        // update function value
        fp = fReturned;
        // store old gradient
        gradOld = grad;
        //compute new gradient
        func.df(xPos, grad);
        // do another convergence check where max(abs(nabla f / xNew )/scale ) is determined
        T gradScale = std::max(fReturned, (T)1);
        T gradNorm = std::inner_product(grad.cbegin(), grad.cend(), grad.cbegin(), T(0));
        if (gradNorm < gradTol) return;
        // compute delta gradient for curvature pair
        std::transform(grad.cbegin(),
                       grad.cend(),
                       gradOld.cbegin(),
                       deltaGrad.begin(),
                       [](const T& gnew, const T& gold) { return gnew - gold; });
        collectCurvaturePairs(iter, deltaGrad, deltaX);
    }

    throw std::runtime_error("Too many iterations in LBFGS<T>::minimize");
}

} //namespace vaspml

#endif
