#ifndef SEQUENTIALLIKELIHOODCOVERAGE_HPP
#define SEQUENTIALLIKELIHOODCOVERAGE_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm>    // for max, remove
#include <cmath>        // for exp, pow, sqrt
#include <limits>       // for numeric_limits
#include <numeric>      // for iota
#include <type_traits>  // for is_arithmetic_v

namespace vaspml::bayes_select
{

template<typename T, typename MetricFunction>
T estimateSilvermanSigmaWithMetric(const std::vector<std::vector<T>>& data,
                                   MetricFunction                     metricFunction)
{
    const Size n = data.size();
    if (n < 2) return 1.0;
    // For the scaling factor, we need the dimensionality
    const Size d = data[0].size();
    // 1. Collect a sample of distances to estimate the spread of the metric space
    // We'll sample O(n) distances instead of O(n^2) to keep it fast
    std::vector<T> distances(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size i = 0; i < n; i++)
    {
        // Measure distance from each point to a "neighbor" or a random reference
        // Using i and (i + 1) % n provides a simple linear chain of distances
        // just trick to get a cheap estimate for variance O(n)
        distances[i] = metricFunction(data[i], data[(i + 1) % n]);
    }

    // 2. Calculate Mean and Variance of these distances in parallel
    T sum = 0.0;
    T sumSq = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum, sumSq)
#endif
    for (Size i = 0; i < n; i++)
    {
        sum += distances[i];
        sumSq += distances[i] * distances[i];
    }

    T meanDist = sum / (T)n;
    T varianceDist = (sumSq / (T)n) - (meanDist * meanDist);
    T sigmaDist = std::sqrt(std::max<T>(0.0, varianceDist));

    // If all distances are the same, fallback to mean distance or a default
    if (sigmaDist <= 1e-9) sigmaDist = (meanDist > 0) ? meanDist : 1.0;

    // 3. Apply Silverman's Scaling
    // This adjusts the "spread" of your metric to the optimal kernel bandwidth
    T exponent = 1.0 / (static_cast<T>(d) + 4.0);
    T base = 4.0 / (static_cast<T>(n) * (static_cast<T>(d) + 2.0));
    T factor = std::pow(base, exponent);

    return factor * sigmaDist;
}

// Precompute symmetric kernel matrix using Gaussian kernel
template<typename T, typename NormFunction>
std::vector<T> computeKernelMatrix(const std::vector<std::vector<T>>& data,
                                   const T                            sigmaIn,
                                   NormFunction                       computeMetric)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );
    const Size     n = data.size();
    T              sigma = sigmaIn;
    std::vector<T> K(n * n, 0.0);
    if (sigma <= 0.0) sigma = estimateSilvermanSigmaWithMetric(data, computeMetric);
    const T denom = 2.0 * sigma * sigma;
    for (Size i = 0; i < n; i++)
    {
        K[i * n + i] = 1.0;
        for (Size j = i + 1; j < n; j++)
        {
            T d2 = computeMetric(data[i], data[j]);
            d2 *= d2;
            const T val = std::exp(-d2 / denom);
            K[i * n + j] = val;
            K[j * n + i] = val;
        }
    }

    return K;
}

// Convenience function to retrieve selected vectors
template<typename T>
std::vector<std::vector<T>> extractSelectedVectors(const std::vector<std::vector<T>>& data,
                                                   const Vec1Int&                     indices)
{
    std::vector<std::vector<T>> subset;
    subset.reserve(indices.size());
    for (Size idx : indices) subset.push_back(data[idx]);

    return subset;
}

// Sequential likelihood coverage selection.
// coverage[i] maintained as: coverage[i] = 1 - product_{s in S} (1 - K[i,s])
// Marginal gain from adding candidate c: sum_i (newCoverage[i] - coverage[i])
template<typename T>
struct SequentialCoverageResult
{
    Vec1Int selectedIndices;
    T       sigmaUsed;
    T       finalCoverageSum;
    Size    iterations;
};

/*******************************************************************************************
 * @brief Selects a subset of elements from the input data to maximize coverage likelihood
 *        using a sequential greedy algorithm.
 *
 * This function implements a sequential likelihood coverage selection algorithm. It selects
 * elements from the input data to maximize a coverage metric, which is computed using a
 * kernel matrix derived from the input data and a user-defined metric function. The selection
 * process stops based on specified criteria, such as minimum marginal gain or relative gain drop.
 *
 * @tparam T Type of the data elements (must be an arithmetic type).
 * @tparam MetricFunction Type of the metric function used to compute pairwise similarity.
 *
 * @param data A 2D vector representing the input data. Each inner vector corresponds to a data
 *point.
 * @param metricFunction A callable object (e.g., function, lambda) that computes the similarity
 *        between two data points. It should accept two data points of type `T` and return a
 *similarity score.
 * @param sigma The kernel width parameter. If set to 0.0, it will be estimated automatically
 *        based on the input data. Default is 0.0.
 * @param minMarginalGain The minimum marginal gain required to continue the selection process.
 *        If the gain falls below this value, the algorithm stops. Default is 1e-6.
 * @param relativeGainDrop The relative drop in gain that triggers stopping. If the gain in the
 *        current iteration is less than this fraction of the previous gain, the algorithm stops.
 *Default is 0.01.
 * @param maxSelectionsIn The maximum number of elements to select. If set to 0, the algorithm
 *        allows selecting up to the total number of elements in the input data. Default is 0.
 *
 * @return A vector of indices representing the selected elements from the input data.
 *
 * @note The function assumes that the input data is non-empty and that the metric function is
 *       well-defined for the given data type. If the data type `T` is not arithmetic, the function
 *       will trigger a debug error.
 *
 * @warning The function modifies the internal state of the algorithm during execution, including
 *          the kernel matrix and coverage probabilities. Ensure that the input data and metric
 *          function are consistent and valid.
 *
 * @details
 * The algorithm works as follows:
 * 1. Compute a kernel matrix using the input data and the metric function.
 * 2. Iteratively select the element that provides the highest marginal gain in coverage.
 * 3. Update the coverage probabilities and remove the selected element from the candidate list.
 * 4. Stop the selection process based on the specified stopping criteria.
 *
 * The kernel matrix is computed using the Gaussian kernel:
 * \f[
 * K(i, j) = \exp\left(-\frac{\text{metricFunction}(data[i], data[j])^2}{2 \sigma^2}\right)
 * \f]
 * where \f$\sigma\f$ is the kernel width parameter.
 *
 * The algorithm ensures that the selected subset maximizes the overall coverage likelihood
 * while adhering to the specified constraints.
 *******************************************************************************************/
template<typename T, typename MetricFunction>
Vec1Int sequentialLikelihoodCoverageSelect(const std::vector<std::vector<T>>& data,
                                           MetricFunction                     metricFunction,
                                           const Size                         maxSelectionsIn = 0,
                                           T                                  sigma = 0.0,
                                           T minMarginalGain = 1e-6,
                                           T relativeGainDrop = 0.01)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );
    SequentialCoverageResult<T> result;
    Size                        n = data.size();
    if (n == 0)
    {
        result.finalCoverageSum = 0.0;
        result.sigmaUsed = sigma;
        result.iterations = 0;
        return {};
    }
    // use some factor here. Otherwise sigma is to small and one get's a Diract delta function as
    // kernel
    if (sigma <= 0.0) sigma = 1000.0 * estimateSilvermanSigmaWithMetric(data, metricFunction);
    //std::cout << "SIGMA SIGMA NO OMP " << sigma << std::endl;
    result.sigmaUsed = sigma;

    std::vector<T> K = computeKernelMatrix(data, sigma, metricFunction);

    std::vector<T> coverage(n, 0.0); // current coverage probabilities
    Vec1Size       candidates(n);
    std::iota(candidates.begin(), candidates.end(), 0);

    T    lastGain = std::numeric_limits<Real>::infinity();
    Size iteration = 0;
    T    maxSelections = maxSelectionsIn;
    if (maxSelectionsIn == 0) maxSelections = n; // allow full if not specified

    while (!candidates.empty() && iteration < maxSelections)
    {
        T    bestGain = -1.0;
        Size bestCandidate = std::numeric_limits<size_t>::max();

        // Evaluate marginal gain for each candidate
        for (Size cIdx = 0; cIdx < candidates.size(); cIdx++)
        {
            Real gain = 0.0;
            Size c = candidates[cIdx];
            for (Size i = 0; i < n; i++)
            {
                T kVal = K[i * n + c];
                //std::cout << kVal << "   " << sigma << std::endl;
                T newCoverage = 1.0 - (1.0 - coverage[i]) * (1.0 - kVal);
                gain += (newCoverage - coverage[i]);
            }
            if (gain > bestGain)
            {
                bestGain = gain;
                bestCandidate = c;
            }
        }
        if (bestCandidate == std::numeric_limits<Size>::max()) break;

        //std::cout << "Best gain = " << bestGain << " minMarginalGain =  " << minMarginalGain << "
        //" << lastGain
        //          << "   " << relativeGainDrop * lastGain << std::endl;
        // Stopping criteria
        if (bestGain < minMarginalGain) { break; }
        if (lastGain < std::numeric_limits<T>::infinity() && bestGain < relativeGainDrop * lastGain)
        {
            break;
        }

        // Accept bestCandidate
        Size c = bestCandidate;
        Size colOffset = c;
        for (Size i = 0; i < n; i++)
        {
            T kVal = K[i * n + colOffset];
            coverage[i] = 1.0 - (1.0 - coverage[i]) * (1.0 - kVal);
        }

        result.selectedIndices.push_back(c);
        lastGain = bestGain;
        iteration++;

        // Remove candidate from list
        candidates.erase(std::remove(candidates.begin(), candidates.end(), c), candidates.end());
    }

    T coverageSum = 0.0;
    for (T v : coverage) coverageSum += v;
    result.finalCoverageSum = coverageSum;
    result.iterations = iteration;

    return result.selectedIndices;
}

template<typename T, typename MetricFunction>
Vec1Int sequentialLikelihoodCoverageSelectOMP(const std::vector<std::vector<T>>& data,
                                              MetricFunction                     metricFunction,
                                              const Size maxSelectionsIn = 0,
                                              T          sigma = 0.0,
                                              T          minMarginalGain = 1e-6,
                                              T          relativeGainDrop = 0.01)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );

    SequentialCoverageResult<T> result;
    const Size                  n = data.size();
    if (n == 0)
    {
        result.finalCoverageSum = 0.0;
        result.sigmaUsed = sigma;
        result.iterations = 0;
        return {};
    }

    if (sigma <= 0.0) sigma = estimateSilvermanSigmaWithMetric(data, metricFunction);
    result.sigmaUsed = sigma;

    //std::cout << "SIGMA SIGMA " << sigma << std::endl;
    // 1. Parallel Kernel Computation (Assumed computeKernelMatrix can be OMP enabled)
    std::vector<T> K = computeKernelMatrix(data, sigma, metricFunction);

    std::vector<T>    coverage(n, 0.0);
    std::vector<char> isSelected(n, 0);
    Vec1Size          candidates(n);
    std::iota(candidates.begin(), candidates.end(), 0);

    T    lastGain = std::numeric_limits<T>::infinity();
    Size iteration = 0;
    Size maxSelections = (maxSelectionsIn == 0) ? n : maxSelectionsIn;

    while (!candidates.empty() && iteration < maxSelections)
    {
        T    bestGain = -1.0;
        Size bestCandidate = std::numeric_limits<Size>::max();

// 2. Parallel Evaluation of Marginal Gain
// We use a custom reduction to find the index and value of the maximum gain
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            T    localBestGain = -1.0;
            Size localBestCandidate = std::numeric_limits<Size>::max();

#ifdef _OPENMP
#pragma omp for nowait
#endif
            for (Size cIdx = 0; cIdx < candidates.size(); cIdx++)
            {
                T    gain = 0.0;
                Size c = candidates[cIdx];

                // Inner loop optimization: compute gain for this candidate
                for (Size i = 0; i < n; i++)
                {
                    if (isSelected[i] && i == c) continue;
                    T kVal = K[i * n + c];
                    T newCoverage = 1.0 - (1.0 - coverage[i]) * (1.0 - kVal);
                    gain += (newCoverage - coverage[i]);
                }

                if (gain > localBestGain)
                {
                    localBestGain = gain;
                    localBestCandidate = c;
                }
            }

// Critical section to find the global best across threads
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (localBestGain > bestGain)
                {
                    bestGain = localBestGain;
                    bestCandidate = localBestCandidate;
                }
            }
        }

        if (bestCandidate == std::numeric_limits<Size>::max()) break;

        //std::cout << "Best gain = " << bestGain << " minMarginalGain =  " << minMarginalGain << "
        //" << lastGain
        //          << "   " << relativeGainDrop * lastGain << std::endl;
        // Stopping criteria
        if (bestGain < minMarginalGain)
        {
            //std::cout << "DOING THE FIRST BREAK " << std::endl;
            break;
        }
        if (lastGain < std::numeric_limits<T>::infinity() && bestGain < relativeGainDrop * lastGain)
        {
            //std::cout << "DOING THE SECOND BREAK " << std::endl;
            break;
        }

        // Accept bestCandidate
        isSelected[bestCandidate] = 1;

// 3. Parallel Update of Coverage vector
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (Size i = 0; i < n; i++)
        {
            T kVal = K[i * n + bestCandidate];
            coverage[i] = 1.0 - (1.0 - coverage[i]) * (1.0 - kVal);
        }

        result.selectedIndices.push_back(bestCandidate);
        lastGain = bestGain;
        iteration++;

        // Remove candidate (Sequential)
        candidates.erase(std::remove(candidates.begin(), candidates.end(), bestCandidate),
                         candidates.end());
    }

    // 4. Parallel Reduction for Final Sum
    T coverageSum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : coverageSum)
#endif
    for (Size i = 0; i < n; i++) { coverageSum += coverage[i]; }

    result.finalCoverageSum = coverageSum;
    result.iterations = iteration;

    return result.selectedIndices;
}

} //namespace vaspml::bayes_select

#endif
