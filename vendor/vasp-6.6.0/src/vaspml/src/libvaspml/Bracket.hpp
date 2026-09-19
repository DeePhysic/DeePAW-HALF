#ifndef BRACKET_HPP
#define BRACKET_HPP

// Avoid IWYU suggesting <stdlib.h> because of integer argument variant of std::abs (not used here).
// IWYU pragma: no_include <stdlib.h>
#include <algorithm>    // for max
#include <cmath>        // for abs
#include <tuple>        // for swap, tuple, make_tuple
#include <type_traits>  // for is_floating_point

namespace vaspml
{

template<typename T>
T sign(const T& value)
{
    return (T(0) < value) - (value < T(0));
}

template<typename T>
class BracketMinimum
{
    static_assert(std::is_floating_point<T>::value,
                  "Unsupported template type for BracketMinimum class\n");

  public:
    /*******************************************************************************************
     * Construct class and initialize all values by zero.
     *******************************************************************************************/
    BracketMinimum();
    /*******************************************************************************************
     * Compute the bracketing values for a splied function.
     *
     * @param a First starting guess x value
     * @param b Second starting guess x value
     *******************************************************************************************/
    template<typename U>
    void bracket(const T& a, const T& b, U& func);
    /*******************************************************************************************
     * Return left bracketing value after bracketing was carried out.
     *******************************************************************************************/
    T get_ax() const;
    /*******************************************************************************************
     * Return central x value after bracketing which is closest to the actual minimum in interval
     *******************************************************************************************/
    T get_bx() const;
    /*******************************************************************************************
     * Return right bracketing value after bracketing was carried out.
     *******************************************************************************************/
    T get_cx() const;
    /*******************************************************************************************
     * Return function value at xa
     *******************************************************************************************/
    T get_fa() const;
    /*******************************************************************************************
     * Return function value at xb
     *******************************************************************************************/
    T get_fb() const;
    /*******************************************************************************************
     * Return function value at xc
     *******************************************************************************************/
    T get_fc() const;
    /*******************************************************************************************
     * Return all x values and function values. The values will be returned as tuple.
     *
     * The return tuple is ordered as xa,xb,xc,fa,fb,fc
     *******************************************************************************************/
    std::tuple<T, T, T, T, T, T> get_result() const;

  private:
    /*******************************************************************************************
     * Function to shift values
     * @param a input will be overwitten by b
     * @param b input will be overwitten by c
     * @param c input will be overwitten by d
     * @param d will over write c
     *******************************************************************************************/
    inline void shft3(T& a, T& b, T& c, const T& d);
    /*******************************************************************************************
     * Left bracketing x value
     *******************************************************************************************/
    T ax;
    /*******************************************************************************************
     * Central value within intervall (a,c). Is closest to actual minimum of function
     *******************************************************************************************/
    T bx;
    /*******************************************************************************************
     * Right bracketing x value
     *******************************************************************************************/
    T cx;
    /*******************************************************************************************
     * Function value at ax.
     *******************************************************************************************/
    T fa;
    /*******************************************************************************************
     * Function value at bx.
     *******************************************************************************************/
    T fb;
    /*******************************************************************************************
     * Function value at cx.
     *******************************************************************************************/
    T fc;
    /*******************************************************************************************
     * Golden ration needed for x update.
     *******************************************************************************************/
    const T GOLD = 1.618034;
    /*******************************************************************************************
     * Constant to compute outer limit of
     *******************************************************************************************/
    const T GLIMIT = 100.0;
    /*******************************************************************************************
     * Tiny value to avoid divison by zero.
     *******************************************************************************************/
    const T TINY = 1.0e-20;
};

template<typename T>
BracketMinimum<T>::BracketMinimum() : ax((T)0), bx((T)0), fa((T)0), fb((T)0), fc((T)0)
{}

template<typename T>
template<typename U>
void BracketMinimum<T>::bracket(const T& a, const T& b, U& func)
{
    T fu;
    ax = a;
    bx = b;
    fa = func(ax);
    fb = func(bx);
    if (fb > fa)
    {
        std::swap(fa, fb);
        std::swap(ax, bx);
    }

    cx = bx + GOLD * (bx - ax);
    fc = func(cx);

    while (fb > fc)
    {
        // compute parabola and it's extremum for x update
        T r = (bx - ax) * (fb - fc);
        T q = (bx - cx) * (fb - fa);
        T u = bx
            - ((bx - ax) * q - (bx - ax) * r)
                  / ((T)2.0 * sign(q - r) * std::max(std::abs(q - r), TINY));
        T ulim = bx + GLIMIT * (cx - bx);
        if ((bx - u) * (u - cx) > 0)
        {
            fu = func(u);
            if (fu < fc)
            {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            }
            else if (fu > fb)
            {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + GOLD * (cx - bx);
        }
        else if ((cx - u) * (u - ulim) > (T)0)
        {
            fu = func(u);
            if (fu < fc)
            {
                shft3(bx, cx, u, u + GOLD * (u - cx));
                shft3(fb, fc, fu, func(u));
            }
        }
        else if ((u - ulim) * (ulim - cx) >= (T)0)
        {
            u = ulim;
            fu = func(u);
        }
        else
        {
            u = cx + GOLD * (cx - bx);
            fu = func(u);
            shft3(ax, bx, cx, u);
            shft3(fa, fb, fc, fu);
        }
    }
}

template<typename T>
inline void BracketMinimum<T>::shft3(T& a, T& b, T& c, const T& d)
{
    a = b;
    b = c;
    c = d;
}

template<typename T>
T BracketMinimum<T>::get_ax() const
{
    return ax;
}

template<typename T>
T BracketMinimum<T>::get_bx() const
{
    return bx;
}

template<typename T>
T BracketMinimum<T>::get_cx() const
{
    return cx;
}

template<typename T>
T BracketMinimum<T>::get_fa() const
{
    return fa;
}

template<typename T>
T BracketMinimum<T>::get_fb() const
{
    return fb;
}

template<typename T>
T BracketMinimum<T>::get_fc() const
{
    return fc;
}

template<typename T>
std::tuple<T, T, T, T, T, T> BracketMinimum<T>::get_result() const
{
    return std::make_tuple(ax, bx, cx, fa, fb, fc);
}

} // namespace vaspml

#endif
