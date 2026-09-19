#ifndef BFGS_HPP
#define BFGS_HPP

#include "Linesearch.hpp"
#include "types.hpp"

// Avoid IWYU suggesting <stdlib.h> because of integer argument variant of std::abs (not used here).
// IWYU pragma: no_include <stdlib.h>
#include <algorithm>   // for max, transform, for_each, copy, fill
#include <cmath>       // for abs, sqrt
#include <functional>  // for minus
#include <limits>      // for numeric_limits
#include <numeric>     // for inner_product, iota
#include <stdexcept>   // for runtime_error
#include <type_traits> // for is_floating_point

namespace vaspml
{

/*******************************************************************************************
 * Minimize the supplied function with a quasi-Newton algorithm in BFGS form.
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
 * where \f$ \mathbf{y}_{k} = \nabla f_{k+1} - \nabla f_{k}\f$ and \f$ s_{k} = x_{k+1} - s_{k}\f$.
 *
 * With this constraints a solution to the update of the Hessian matrix is given by and
 *\f[
 *    \mathbf{H}_{k+1}=\mathbf{H}_{k}+(I-\rho_{k}\mathbf{s}_{k}\mathbf{y}^{T}_{k})\mathbf{H}_{k}
 *    (I-\rho_{k}\mathbf{s}_{k}\mathbf{y}^{T}_{k})+\rho_{k}\mathbf{s}_{k}\mathbf{s}^{T}_{k},
 *\f]
 * where \f$ \rho_{k}=\frac{1/}{\mathbf{s}^{T}_{k}\mathbf{y}_{k}}\f$. For computation
 * of the inverse Hessian the following update function can be used
 *\f[
 *  \mathbf{B}_{k+1}=\mathbf{B}_{k}-
 *  \frac{\mathbf{B}_{k}\mathbf{s}_{k}\mathbf{s}^{T}_{k}\mathbf{B}}{\mathbf{s}^{T}_{k}\mathbf{s}_{k}}+
 *  \frac{\mathbf{y}_{k}\mathbf{y}_{k}^{T}}{\mathbf{y}^{T}_{k}\mathbf{s}_{k}}
 *\f]
 *The implementation here uses the update of inverse matrix directly because a matrix
 *inversion can be safed.
 *******************************************************************************************/
template<typename T>
class BFGS
{
    static_assert(std::is_floating_point<T>::value, "Unsupported template type for BFGS class\n");

  public:
    BFGS();
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
    void minimize(std::vector<T>& xPos, const T& gradTol, Size& nIter, T& fReturned, U& func);
    /*******************************************************************************************
     * Initialize Local arrays needed for minimize function. Initializes Hessian matrix as
     * unit matrix.
     *******************************************************************************************/
    void init(const Size& dim);
    /*******************************************************************************************
     * Return approximate to hessian matrix.
     *******************************************************************************************/
    const std::vector<std::vector<T>>& get_hessian() const;

  private:
    /*******************************************************************************************
     * Maximal number of iterations to find function minimum
     *******************************************************************************************/
    const Size MAXITER = 1000;
    /*******************************************************************************************
     * Smallest numeric value of T which can be distinguished from zero
     *******************************************************************************************/
    const T EPS = std::numeric_limits<T>::epsilon();
    /*******************************************************************************************
     * Check if change in x is converged.
     *******************************************************************************************/
    const T TOLX = (T)4 * EPS;
    /*******************************************************************************************
     * Scaling constant to determine maximal step width
     *******************************************************************************************/
    const T STPMX = (T)100;
    /*******************************************************************************************
     * Hessian matrix which is updated during the quasi newton iterations.
     *******************************************************************************************/
    std::vector<std::vector<T>> hessian;
    /*******************************************************************************************
     * Difference between actual function gradient and previous function gradient
     *******************************************************************************************/
    std::vector<T> deltaGrad;
    /*******************************************************************************************
     * Actual gradient of function.
     *******************************************************************************************/
    std::vector<T> grad;
    /*******************************************************************************************
     * Stores function gradient of previous step.
     *******************************************************************************************/
    std::vector<T> gradOld;
    /*******************************************************************************************
     * Hessian times deltaGrad
     *******************************************************************************************/
    std::vector<T> hDeltaGrad;
    /*******************************************************************************************
     * Newton direction
     *******************************************************************************************/
    std::vector<T> newtDir;
    /*******************************************************************************************
     * Stores x positions during the minimization
     *******************************************************************************************/
    std::vector<T> xi;
    /*******************************************************************************************
     * Number of dimension of the gradient. Number of independet variables of the function to
     * be minimized.
     *******************************************************************************************/
    Size nDim = 0;
    /*******************************************************************************************
     * Index array which is defined as an iota to fill the Hessian as unit matrix.
     *******************************************************************************************/
    Vec1Size index;
    /*******************************************************************************************
     * Line search method to optimize the function value along the new Newton direction.
     *******************************************************************************************/
    InterpolLinesearch<T> linesearch;
    //ArmijoWolfeLinesearch<T> linesearch;
};

template<typename T>
BFGS<T>::BFGS()
{}

template<typename T>
void BFGS<T>::init(const Size& dim)
{
    if (dim != nDim)
    {
        nDim = dim;
        hessian.resize(nDim);
        std::for_each(hessian.begin(),
                      hessian.end(),
                      [&](std::vector<T>& slice) { slice.resize(nDim); });
        deltaGrad.resize(nDim);
        grad.resize(nDim);
        gradOld.resize(nDim);
        hDeltaGrad.resize(nDim);
        newtDir.resize(nDim);
        xi.resize(nDim);
        index.resize(nDim);
        std::iota(index.begin(), index.end(), 0);
        std::for_each(index.cbegin(),
                      index.cend(),
                      [&](const Size& indx)
                      {
                          std::fill(hessian[indx].begin(), hessian[indx].end(), (T)0);
                          hessian[indx][indx] = (T)1;
                      });
    }
    else
    {
        std::for_each(index.cbegin(),
                      index.cend(),
                      [&](const Size& indx)
                      {
                          std::fill(hessian[indx].begin(), hessian[indx].end(), (T)0);
                          hessian[indx][indx] = (T)1;
                      });
    }
}

template<typename T>
template<typename U>
void BFGS<T>::minimize(std::vector<T>& xPos, const T& gradTol, Size& nIter, T& fReturned, U& func)
{
    bool check;
    init(xPos.size());
    T fp = func(xPos);
    func.df(xPos, grad);

    // initial search direction is set to the gradient
    std::transform(grad.cbegin(), grad.cend(), newtDir.begin(), [](const T& gi) { return -gi; });
    // compute the length of the input vector x to determine the maximal step length
    T sum = std::inner_product(xPos.cbegin(), xPos.cend(), xPos.cbegin(), (T)0);
    T stpMax = STPMX * std::max(std::sqrt(sum), (T)nDim);
    for (Size iter = 0; iter < MAXITER; iter++)
    {
        nIter = iter;
        linesearch(xPos, fp, grad, newtDir, xi, fReturned, stpMax, check, func);
        // after line search xi contains xNew
        // and newtDir contains the Newton direction
        // update function value
        fp = fReturned;
        // compute the new search direction p = xNew - xOld
        std::transform(xi.cbegin(), xi.cend(), xPos.cbegin(), newtDir.begin(), std::minus<>());
        // update current position
        std::copy(xi.cbegin(), xi.cend(), xPos.begin());
        // check for convergence max(abs(p/xNew)) it is not the max newton direction tested
        // but a a version relative to the xvalues of the function input to avoid scaling problems
        // TODO: Are the arguments of the reduce lambda here in the wrong order?
        //T convergenceTest = std::transform_reduce( newtDir.cbegin(), newtDir.cend(),
        //                                xPos.cbegin(), (T) 0,
        //                                []( const T& a, const T& b)
        //                                {
        //                                    return std::max(a, b);
        //                                },
        //                                [](const T& x, const T& dir)
        //                                {
        //                                    return std::abs(x) / std::max( std::abs(dir), (T)1);
        //                                });
        T convergenceTest = (T)0;
        for (Size i = 0; i < newtDir.size(); ++i)
        {
            convergenceTest =
                std::max(convergenceTest, std::abs(newtDir[i]) / std::max(std::abs(xPos[i]), (T)1));
        }
        // return if convergence was reached
        if (convergenceTest < TOLX) return;
        // store gradient
        std::copy(grad.cbegin(), grad.cend(), gradOld.begin());
        // compute new gradient
        func.df(xPos, grad);
        convergenceTest = (T)0;
        // do another convergence check where max(abs(nabla f / xNew )/scale ) is determined
        T gradScale = std::max(fReturned, (T)1);
        //convergenceTest = std::transform_reduce( grad.cbegin(), grad.cend(),
        //                              xPos.cbegin(), (T)0,
        //                              [](const T& a, const T& b )
        //                              {
        //                                  return std::max(a, b);
        //                              },
        //                             [&gradScale](const T& gradi, const T& dir)
        //                             {
        //                                 return std::abs(gradi) * std::max(std::abs( dir ), 1.0) /
        //                                 gradScale;
        //                             });
        for (Size i = 0; i < newtDir.size(); ++i)
        {
            convergenceTest =
                std::max(convergenceTest,
                         std::abs(grad[i]) * std::max(std::abs(xPos[i]), 1.0) / gradScale);
        }
        // return if gradient at new position x is small enough
        if (convergenceTest < gradTol) return;
        // compute the differnce of the new and old gradient in BFGS literature notation this is
        // y_{k}
        std::transform(grad.cbegin(),
                       grad.cend(),
                       gradOld.begin(),
                       deltaGrad.begin(),
                       std::minus<>());
        // Here the updating part for the hessian matrix starts.
        // compute B_{k} * y_{k}
        std::transform(
            hessian.cbegin(),
            hessian.cend(),
            hDeltaGrad.begin(),
            [&](const std::vector<T>& row)
            { return std::inner_product(row.cbegin(), row.cend(), deltaGrad.cbegin(), (T)0); });
        // computing the scalar denominators for the BFGS updating scheme
        // compute (x_{i+1}-x_{i})*(\nabla f_{k+1}-\nabla_{f_{k}})=y_{k}*p_{k}
        T rho = std::inner_product(deltaGrad.cbegin(), deltaGrad.cend(), newtDir.cbegin(), (T)0);
        // compute (\nabla f_{k+1}-\nabla_{f_{k}})*H*(\nabla f_{k+1}-\nabla_{f_{k}})=y_{k}H_{k}y_{k}
        T yHy = std::inner_product(deltaGrad.cbegin(), deltaGrad.cend(), hDeltaGrad.cbegin(), (T)0);
        //compute the square of y_{k}*y_{k}
        T yk2 = std::inner_product(deltaGrad.cbegin(), deltaGrad.cend(), deltaGrad.cbegin(), (T)0);
        // compute the length of the used update direction. If the update is too small and hence x
        // does not change much then one can skip the update of the hessian since the old
        // approximate is still good enough
        T sumnewtDir = std::inner_product(newtDir.cbegin(), newtDir.cend(), newtDir.cbegin(), T(0));
        if (rho > std::sqrt(EPS * yk2
                            * sumnewtDir)) // skip hessian update if not enough decrease of function
        {
            rho = (T)1 / rho;
            T fad = (T)1 / yHy;
            // compute p_{k} / y_{k}*p_{k} - H*y_{k} / yHy this is Eq 10.9.10 in numerical recipes
            // for c++
            std::transform(newtDir.cbegin(),
                           newtDir.cend(),
                           hDeltaGrad.cbegin(),
                           deltaGrad.begin(),
                           [&rho, &fad](const T& dir, const T& sk)
                           { return rho * dir - fad * sk; });
            // update the hessian matrix with the BFGS formulae shown in class description
            // this is shown in equation 6.19 in Nocedal numerical Optimization. But implementation
            // equation in this for each is from numerical recipes in c++
            std::for_each(hessian.begin(),
                          hessian.end(),
                          [&](std::vector<T>& row)
                          {
                              Size i = &row - &hessian[0];
                              std::for_each(row.begin() + i,
                                            row.end(),
                                            [&](T& element)
                                            {
                                                Size j = &element - &row[0];
                                                element +=
                                                    rho * newtDir[i] * newtDir[j]
                                                    - fad * hDeltaGrad[i] * hDeltaGrad[j] +
                                                    // this term is taken from numerical recipes.
                                                    // This is not part of the derivation of BFGS
                                                    // but is added to improve numerical stability
                                                    yHy * deltaGrad[i] * deltaGrad[j];
                                                if (i != j) { hessian[j][i] = element; }
                                            });
                          });
        }
        // compute next direction by p_{k} = - (H_{k})^{-1} * g_{k}
        std::transform(
            hessian.cbegin(),
            hessian.cend(),
            newtDir.begin(),
            [&](const std::vector<T>& row)
            { return -std::inner_product(row.cbegin(), row.cend(), grad.cbegin(), (T)0); });
    }

    throw std::runtime_error("!!Too many iterations in BFGS<T>::minimize!!");
}

template<typename T>
const std::vector<std::vector<T>>& BFGS<T>::get_hessian() const
{
    return hessian;
}

} //namespace vaspml

#endif
