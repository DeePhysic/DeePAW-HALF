#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "types.hpp"

#include <chrono>       // for duration, high_resolution_clock
#include <cmath>        // for log, sqrt
#include <limits>       // for numeric_limits
#include <random>       // for mt19937
#include <type_traits>  // for is_floating_point_v

namespace vaspml
{

/*******************************************************************************************
 * @class UniqueRandomizer
 * @brief A utility class for generating unique random numbers within a specified range.
 *
 * Function draws count number of random numbers in a supplied range without repetition.
 * If a next draw is done where count2 random numbers are choosen. The set will have
 * not have repetitions within itself but can have to the previously drawn set.
 *******************************************************************************************/
class UniqueRandomizer
{
  public:
    UniqueRandomizer();
    /*******************************************************************************************
     * @brief Constructs a UniqueRandomizer with a single seed value.
     * @param seed The seed value used to initialize the random number generator.
     *******************************************************************************************/
    UniqueRandomizer(const Size seed);
    /*******************************************************************************************
     * @brief Constructs a UniqueRandomizer with a sequence of seed values.
     * @param seed A vector of integers used to initialize the random number generator.
     *******************************************************************************************/
    UniqueRandomizer(const Vec1Int& seed);
    /*******************************************************************************************
     * @brief Reinitializes the random number generator with a single seed value.
     * @param seed The new seed value used to reinitialize the random number generator.
     *******************************************************************************************/
    void init(const Size seed);
    /*******************************************************************************************
     * @brief Reinitializes the random number generator with a sequence of seed values.
     * @param seed A vector of integers used to reinitialize the random number generator.
     *******************************************************************************************/
    void init(const Vec1Int& seed);
    /*******************************************************************************************
     * @brief Constructs a UniqueRandomizer with a specified range and a random seed.
     * @param min The minimum value of the range (inclusive).
     * @param max The maximum value of the range (inclusive).
     *******************************************************************************************/
    UniqueRandomizer(const Int min, const Int max);
    /*******************************************************************************************
     * @brief Constructs a UniqueRandomizer with a specified range and a user-provided seed.
     * @param min The minimum value of the range (inclusive).
     * @param max The maximum value of the range (inclusive).
     * @param seed The seed value for the random number generator.
     *******************************************************************************************/
    UniqueRandomizer(const Int min, const Int max, const Size seed);
    /*******************************************************************************************
     * @brief Reinitializes the range of the randomizer with new minimum and maximum values.
     * @param newMin The new minimum value of the range (inclusive).
     * @param newMax The new maximum value of the range (inclusive).
     * @note If `newMin` is greater than `newMax`, the range is not updated, and a warning is
     *logged.
     *******************************************************************************************/
    void reinitializeRange(const Int newMin, const Int newMax);
    /*******************************************************************************************
     * @brief Draws a specified number of unique random numbers from the current range.
     * @param count The number of unique random numbers to draw.
     * @return A vector of unique random numbers. Returns an empty vector if `count` exceeds the
     *range size.
     * @note The order of the numbers in the returned vector is randomized.
     *******************************************************************************************/
    Vec1Int draw(const Int count);

  private:
    /*******************************************************************************************
     * @brief Initializes the internal number sequence based on the specified range.
     * @param min The minimum value of the range (inclusive).
     * @param max The maximum value of the range (inclusive).
     * @note This method clears and repopulates the internal vector with the sequence `[min, ...,
     *max]`.
     *******************************************************************************************/
    void initializeNumbers(const Int min, const Int max);
    /*******************************************************************************************
     * minimum value of the range (inclusive).
     *******************************************************************************************/
    Int minVal;
    /*******************************************************************************************
     * maximum value of the range.
     *******************************************************************************************/
    Int maxVal;
    /*******************************************************************************************
     * The vector containing the complete range of numbers
     *******************************************************************************************/
    Vec1Int numbers;
    /*******************************************************************************************
     * @brief Pseudo-random number generator based on the Mersenne Twister algorithm.
     *
     * This instance of `std::mt19937` is used to generate high-quality random numbers.
     * It is a 32-bit Mersenne Twister engine.
     *******************************************************************************************/
    std::mt19937 generator;
};

/*******************************************************************************************
 * @class UniformRandomIntGenerator
 * @brief A utility class for generating uniformly distributed random integers.
 *******************************************************************************************/
class UniformRandomIntGenerator
{
  public:
    /*******************************************************************************************
     * @brief Constructs a UniformRandomIntGenerator and seeds the random number generator.
     *
     * Initializes the random number generator with a seed obtained from a random device.
     *******************************************************************************************/
    UniformRandomIntGenerator();
    /*******************************************************************************************
     * @brief Constructs a UniformRandomIntGenerator with a specified seed.
     *
     * @param seed The seed value used to initialize the random number generator.
     *******************************************************************************************/
    UniformRandomIntGenerator(const Size seed);
    /*******************************************************************************************
     * @param seed A vector of integers used to initialize the random number generator.
     *******************************************************************************************/
    UniformRandomIntGenerator(const Vec1Int& seed);
    /*******************************************************************************************
     * @param seed The seed value to reinitialize the random number generator.
     *******************************************************************************************/
    void init(const Size seed);
    /*******************************************************************************************
     * @param seed A vector of integers used to reinitialize the random number generator.
     *******************************************************************************************/
    void init(const Vec1Int& seed);
    /*******************************************************************************************
     * @brief Generates a random integer within a specified range.
     *
     * @param low The lower bound of the range (inclusive).
     * @param high The upper bound of the range (inclusive).
     * @return A random integer between low and high.
     *
     * @note Assumes that a random number generator `generator` is globally defined.
     *******************************************************************************************/
    Int uniformInt(Int low, Int high);

  private:
    /*******************************************************************************************
     * @brief Pseudo-random number generator based on the Mersenne Twister algorithm.
     *
     * This instance of `std::mt19937` is used to generate high-quality random numbers.
     * It is a 32-bit Mersenne Twister engine.
     *******************************************************************************************/
    std::mt19937 generator;
};

typedef unsigned long long int Ullong;

/*******************************************************************************************
 * @class Random
 * @brief A random number generator class based on Numerical Recipes for C++.
 *
 * This class generates random numbers with a period of approximately ~3.138 * 10^57.
 * It supports generating 64-bit integers, 32-bit integers, uniform floating-point numbers,
 * Gaussian-distributed numbers, exponential-distributed numbers, and gamma-distributed numbers.
 *
 * @tparam T Floating-point type (e.g., float, double). Non-floating-point types are not allowed.
 *
 * The routines where adapted from: Numerical Recipes in C++: The Art of Scientific Computing
 * by William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery
 *******************************************************************************************/
template<typename T>
class Random
{
    // Excluding non-floating-point types.
    static_assert(std::is_floating_point_v<T>, "Random requires a floating-point template type.");

  public:
    /*******************************************************************************************
     * @brief Constructor to initialize the random number generator with a seed.
     *
     * @param j Seed value (unsigned long long).
     *******************************************************************************************/
    Random(Ullong j) : v(4101942887655102017LL), w(1)
    {
        u = j ^ v;
        int64();
        v = u;
        int64();
        w = v;
        int64();
    };
    /*******************************************************************************************
     * @brief Constructor for the Random class.
     *
     * Initializes the random number generator using the current high-resolution
     * time as a seed. Sets up internal state variables `u`, `v`, and `w` for
     * generating random numbers.
     *
     * @details
     * - `v` is initialized to a fixed constant.
     * - `w` is initialized to 1.
     * - The current time in nanoseconds is used to seed `u`.
     * - Calls `int64()` three times to further randomize the internal state.
     *******************************************************************************************/
    Random() : v(4101942887655102017LL), w(1)
    {
        const auto now = std::chrono::high_resolution_clock::now();
        Ullong     j = now.time_since_epoch().count();
        u = j ^ v;
        int64();
        v = u;
        int64();
        w = v;
        int64();
    };
    /*******************************************************************************************
     * @brief Generates a 64-bit random integer.
     *
     * @return A 64-bit unsigned integer.
     *
     * @note This function updates the internal state variables to produce a new random number.
     *******************************************************************************************/
    inline Ullong int64()
    {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17;
        v ^= v << 31;
        v ^= v >> 8;
        w = 429495766U * (w & 0xffffffff) + (w >> 32);
        Ullong x = u ^ (u << 21);
        x ^= x >> 35;
        x ^= x << 4;
        return (x + v) ^ w;
    };
    /*******************************************************************************************
     * @brief Generates a uniformly distributed floating-point random number in the interval [0, 1).
     *
     * @return A floating-point random number.
     *******************************************************************************************/
    inline T uniformFloat() { return 5.42101086242752217E-20 * int64(); };
    /*******************************************************************************************
     * @brief Generates a 32-bit random integer.
     *
     * @return A 32-bit unsigned integer.
     *******************************************************************************************/
    inline unsigned int int32() { return (unsigned int)int64(); };
    /*******************************************************************************************
     * @brief Generates a random integer uniformly distributed within a specified range.
     *
     * This function returns a random integer in the range [nlow, nhigh], inclusive.
     * If `nlow` is greater than `nhigh`, their values are swapped to ensure the range is valid.
     * The function avoids bias by rejecting out-of-range random values.
     *
     * @param nlow The lower bound of the range (inclusive).
     * @param nhigh The upper bound of the range (inclusive).
     * @return A random integer uniformly distributed in the range [nlow, nhigh].
     *******************************************************************************************/
    inline Int uniformInt(Int nlow, Int nhigh)
    {
        // Ensure nlow <= nhigh
        if (nlow > nhigh) { std::swap(nlow, nhigh); }

        // Calculate the range size
        Int range = nhigh - nlow + 1;

        // Calculate the maximum acceptable value to avoid bias
        Size maxAcceptable = (std::numeric_limits<Size>::max() / range) * range;

        Size randomValue;
        do {
            randomValue = int32(); // Assume int32() generates a random 32-bit unsigned integer
        }
        while (randomValue >= maxAcceptable); // Reject out-of-range values

        // Scale the random value to the desired range
        Int randomInt = nlow + static_cast<Int>(randomValue % range);
        return randomInt;
    }
    /*******************************************************************************************
     * @brief Generates a Gaussian-distributed random number.
     *
     * @param mu Mean of the distribution.
     * @param sigma Standard deviation of the distribution.
     * @return A Gaussian-distributed random number.
     *******************************************************************************************/
    inline T generateGaussDist(const T mu, const T sigma)
    {
        T                 u1;
        T                 u2;
        T                 z1;
        static T          z2;
        T                 s;
        thread_local bool generate_new;
        generate_new = !generate_new;
        //sigma *= ( double ) 2549.876663134 ;
        if (!generate_new) { return z2 * sigma + mu; }

        do {
            u1 = uniformFloat() * 2.0 - 1.0;
            u2 = uniformFloat() * 2.0 - 1.0;
            s = u1 * u1 + u2 * u2;
        }
        while ((s >= 1.0) || (s == 0.0));

        s = std::sqrt(-2.0 * std::log(s) / s);

        z1 = u1 * s;
        z2 = u2 * s;

        return z1 * sigma + mu;
    }
    /*******************************************************************************************
     * @brief Generates an exponentially distributed random number.
     *
     * @param lambda Rate parameter of the exponential distribution.
     * @return An exponentially distributed random number.
     *******************************************************************************************/
    inline T generateExpDist(const T lambda)
    {
        T u;
        u = uniformFloat();
        return -std::log(1.0 - u) / lambda;
    }
    /*******************************************************************************************
     * @brief Generates a gamma-distributed random number.
     *
     * This method is implemented based on the paper:
     * "A Simple Method for Generating Gamma Variables" by George Marsaglia and Wai Wan Tsang,
     * ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000, Pages 363–372.
     *
     * @param a Shape parameter of the gamma distribution.
     * @return A gamma-distributed random number.
     *******************************************************************************************/
    inline T generateGammaDist(const T a)
    {
        T c;
        T d;
        T x;
        T v;
        T u;
        T value;

        d = a - 1.0 / 3.0;
        c = 1.0 / std::sqrt(9.0 * d);

        // initialize
        u = uniformFloat();
        x = generateGaussDist(0.0, 1.0);
        v = computeVGamma(c, x);
        do {
            x = generateGaussDist(0.0, 1.0);
            v = computeVGamma(c, x);
        }
        while (v <= 0);

        value = v * d;
        do {
            u = uniformFloat();
            do {
                x = generateGaussDist(0.0, 1.0);
                v = computeVGamma(c, x);
            }
            while (v <= 0);
            value = v * d;
        }
        while (u > computeTerm1Gamma(x) || std::log(u) > computeTerm2Gamma(d, v, x));
        return value;
    }

  private:
    Ullong u; ///< Internal state variable.
    Ullong v; ///< Internal state variable.
    Ullong w; ///< Internal state variable.

    /*******************************************************************************************
     * @brief Helper function to compute intermediate value v for gamma distribution.
     *
     * @param c Constant derived from the shape parameter.
     * @param x Gaussian random variable.
     * @return Computed value of v.
     *******************************************************************************************/
    T computeVGamma(const T c, const T x)
    {
        T results = 1.0 + c * x;
        return results * results * results;
    }
    /*******************************************************************************************
     * @brief Helper function to compute the first acceptance condition for gamma distribution.
     *
     * @param x Gaussian random variable.
     * @return Computed value for the first acceptance condition.
     *******************************************************************************************/
    T computeTerm1Gamma(const T x) { return 1.0 - 0.331 * (x * x) * (x * x); }
    /*******************************************************************************************
     * @brief Helper function to compute the second acceptance condition for gamma distribution.
     *
     * @param d Derived constant from the shape parameter.
     * @param v Intermediate value v.
     * @param x Gaussian random variable.
     * @return Computed value for the second acceptance condition.
     *******************************************************************************************/
    T computeTerm2Gamma(const T d, const T v, const T x)
    {
        return 0.5 * (x * x) + d - d * v + d * log(v);
    }
};

} //namespace vaspml

#endif
