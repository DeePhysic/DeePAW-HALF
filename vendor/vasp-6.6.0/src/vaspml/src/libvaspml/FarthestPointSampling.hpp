#ifndef FARTHESTPOINTSAMPLING_HPP
#define FARTHESTPOINTSAMPLING_HPP

#include "MlMPI.hpp"

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm>    // for min
#include <limits>       // for numeric_limits
#include <numeric>      // for iota
#include <type_traits>  // for is_arithmetic_v


namespace vaspml
{
template<class T>
class ShmemArray2D;
}

namespace vaspml::fps
{

/*******************************************************************************************
 * @brief Extracts a subset of vectors from a 2D vector based on specified indexes.
 *
 * This function takes a list of indexes and a 2D vector, and returns a new 2D vector
 * containing only the vectors at the specified indexes.
 *
 * @tparam T The type of elements in the 2D vector.
 * @param indexes A vector of integers specifying the indexes of the vectors to extract.
 * @param data A 2D vector from which the vectors will be extracted.
 * @return A 2D vector containing the extracted vectors in the order specified by `indexes`.
 *******************************************************************************************/
template<typename T>
std::vector<std::vector<T>> extractVectors(const Vec1Int&                     indexes,
                                           const std::vector<std::vector<T>>& data)
{
    std::vector<std::vector<T>> extract(indexes.size());
    for (Size i = 0; i < indexes.size(); i++) extract[i] = data[indexes[i]];
    return extract;
}
/*******************************************************************************************
 * @brief Selects the first index randomly using OpenMP.
 *
 * This function generates a random index in the range [0, nPoints-1] using a
 * uniform random integer generator. If the provided seed is 0, a random seed
 * is generated using `std::random_device`.
 *
 * @param nPoints The total number of points (upper bound for the random index).
 * @param seedIn The seed value for the random number generator. If 0, a random
 *               seed is generated.
 * @return A randomly chosen index in the range [0, nPoints-1].
 *******************************************************************************************/
Int chooseFirstIndexOMP(const Size nPoints, const Size seedIn = 0);
/*******************************************************************************************
 * @brief Selects the first index randomly using MPI for distributed environments.
 *
 * This function generates a random index in the range [0, nPoints-1] using a
 * uniform random integer generator. If the provided seed is 0 and the current
 * MPI rank is 0, a random seed is generated using `std::random_device`. The
 * chosen index is then broadcasted to all MPI ranks.
 *
 * @param nPoints The total number of points (upper bound for the random index).
 * @param seedIn The seed value for the random number generator. If 0, a random
 *               seed is generated on rank 0.
 * @param mlmpi A shared pointer to an MlMPI object for MPI communication. If
 *              provided, the chosen index is broadcasted to all ranks.
 * @return A randomly chosen index in the range [0, nPoints-1].
 *******************************************************************************************/
Int chooseFirstIndexMPI(const Size                          nPoints,
                        const Size                          seedIn = 0,
                        const std::shared_ptr<const MlMPI>& mlmpi = nullptr);
/*******************************************************************************************
 * @brief Performs Farthest Point Sampling (FPS) on a set of points.
 *
 * This function selects a subset of points from the input set such that the selected points
 * are as far apart as possible in terms of Euclidean distance. The algorithm starts with a
 * randomly chosen point and iteratively selects the farthest point from the already selected
 * set until the desired number of samples is reached.
 *
 * @tparam T The type of the point coordinates. Must be an arithmetic type (e.g., float, double).
 * @param points A 2D vector representing the input points. Each inner vector corresponds to a
 *point.
 * @param numSamples The number of points to sample. Must be less than or equal to the size of
 *`points`.
 * @return A tuple containing:
 *         - A 2D vector of the sampled points.
 *         - A vector of indices corresponding to the sampled points in the original input.
 *
 * @note If `numSamples` is 0 or `points` is empty, the function returns an empty result.
 * @note If `numSamples` is greater than or equal to the size of `points`, the function returns all
 *points.
 * @note The function only supports arithmetic types for the point coordinates. Non-arithmetic types
 *will trigger an error.
 *
 * @throws std::runtime_error If the input type `T` is not arithmetic.
 *
 * @warning Ensure that the `euclideanDistance` function is defined and computes the Euclidean
 *distance between two points of type `std::vector<T>`.
 *
 * @example
 * std::vector<std::vector<float>> points = {{0.0f, 0.0f}, {1.0f, 1.0f}, {2.0f, 2.0f},
 *{3.0f, 3.0f}}; auto [sampledPoints, indices] = farthestPointSampling(points, 2);
 * // sampledPoints will contain 2 points farthest apart, and indices will contain their original
 *indices.
 *******************************************************************************************/
template<typename T, typename NormFunction>
Vec1Int farthestPointSamplingKPointsOMP(const std::vector<std::vector<T>>& points,
                                        Size                               numSamples,
                                        NormFunction                       computeDist,
                                        const Size                         seed = 0)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );
    Vec1Int sampledIndices;
    sampledIndices.reserve(numSamples);

    Int firstIndex = chooseFirstIndexOMP(points.size(), seed);
    sampledIndices.push_back(firstIndex);
    std::vector<T> minDistances(points.size(), std::numeric_limits<T>::max());
    minDistances[firstIndex] = -1.0;

// Compute initial distances to the first selected point
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size i = 0; i < points.size(); i++)
    {
        if (minDistances[i] >= 0.0)
        {
            minDistances[i] = computeDist(points[i], points[firstIndex]);
        }
    }
    for (Size sampleIdx = 1; sampleIdx < numSamples; sampleIdx++)
    {
        T   maxMinDist = -1;
        Int nextFarthestIdx = -1;
// Find the farthest point (and update minDistances in-place)
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            T           threadMaxDist = -1;
            Int         threadIdx = -1;
            const auto& lastPointAdded = points[sampledIndices.back()];
// Shared view for the next point's distances
#ifdef _OPENMP
#pragma omp for
#endif
            for (Size j = 0; j < points.size(); j++)
            {
                if (minDistances[j] >= 0)
                {
                    T dist = computeDist(points[j], lastPointAdded);

                    if (dist < minDistances[j]) minDistances[j] = dist;
                    if (minDistances[j] > threadMaxDist)
                    {
                        threadMaxDist = minDistances[j];
                        threadIdx = static_cast<Int>(j);
                    }
                }
            }
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (threadMaxDist > maxMinDist)
                {
                    maxMinDist = threadMaxDist;
                    nextFarthestIdx = threadIdx;
                }
            }
        }
        if (nextFarthestIdx == -1) break;
        sampledIndices.push_back(nextFarthestIdx);
        minDistances[nextFarthestIdx] = -1.0;
    }

    return sampledIndices;
}

/*******************************************************************************************
 * @struct MaxCandidate for finding MPI maxDistance and associated nextFarthestIndex in
 * Allreduce
 * @brief Represents a candidate with the maximum minimum distance and its associated index.
 *******************************************************************************************/
struct MaxCandidate
{
    Real maxMinDist;
    Int  nextFarthestIdx;
};
/*******************************************************************************************
 * @brief Performs farthest point sampling (FPS) to select a specified number of points from a
 *distributed dataset.
 *
 * This function implements the farthest point sampling algorithm in a distributed memory
 *environment using MPI. It selects `numSamples` points from the input dataset `points` such that
 *the selected points are as far apart as possible based on a given distance metric. The algorithm
 *is parallelized across MPI ranks.
 *
 * @tparam T The data type of the points. Must be an arithmetic type.
 * @tparam NormFunction A callable object or function that computes the distance between two points.
 *
 * @param points A 2D array of points, where each row represents a point in the dataset.
 * @param numSamples The number of points to sample. If `numSamples` is 0, an empty result is
 *returned. If `numSamples` is greater than or equal to the number of points, all points are
 *returned.
 * @param computeDist A callable object or function that computes the distance between two points.
 *                    It should take two pointers to the data of the points and the dimensionality
 *as arguments.
 * @param mlmpi A shared pointer to an MPI wrapper object (`MlMPI`) for managing MPI communication.
 *              If `nullptr`, the function assumes a single-process execution.
 * @param seed An optional seed for random number generation. If `0`, a random seed is generated on
 *rank 0 and broadcasted to all ranks.
 *
 * @return A vector of indices representing the sampled points. The indices correspond to the rows
 *in the input `points`.
 *
 * @note The function assumes that the input dataset is distributed across MPI ranks. Each rank is
 *responsible for a subset of the points. The algorithm uses MPI collective operations to
 *synchronize and determine the farthest points globally.
 *
 * @warning This function only supports arithmetic types for the data type `T`. If `T` is not
 *arithmetic, a runtime error is triggered.
 *
 * @details
 * The algorithm proceeds as follows:
 * 1. Randomly select the first point.
 * 2. Iteratively select the farthest point from the already selected points.
 * 3. Use MPI to synchronize and determine the globally farthest point at each iteration.
 * 4. Repeat until the desired number of points is selected.
 *
 * Example usage:
 * @code
 * ShmemArray2D<float> points = ...; // Initialize distributed dataset
 * Size numSamples = 10;
 * auto computeDist = [](const float* p1, const float* p2, Size dims) {
 *     float dist = 0.0;
 *     for (Size i = 0; i < dims; ++i) {
 *         dist += (p1[i] - p2[i]) * (p1[i] - p2[i]);
 *     }
 *     return dist;
 * };
 * auto mlmpi = std::make_shared<MlMPI>();
 * Vec1Int sampledIndices = farthestPointSamplingKPointsMPI(points, numSamples, computeDist, mlmpi);
 * @endcode
 *******************************************************************************************/
template<typename T, typename NormFunction>
Vec1Int farthestPointSamplingKPointsMPI(const ShmemArray2D<T>&              points,
                                        Size                                numSamples,
                                        NormFunction                        computeDist,
                                        const std::shared_ptr<const MlMPI>& mlmpi = nullptr,
                                        const Size                          seed = 0)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );
    // handling base cases
    if (numSamples == 0 || points.get_size0() == 0) { return {}; }
    if (numSamples >= points.get_size0())
    {
        Vec1Int indexes(points.get_size0());
        std::iota(indexes.begin(), indexes.end(), 0);
        return indexes;
    }
    // Stores the indices of the sampled points
    Vec1Int sampledIndices;
    sampledIndices.reserve(numSamples);

    Int        firstIndex = chooseFirstIndexMPI(points.get_size0(), seed, mlmpi);
    Int        nextFarthestIdx = firstIndex;
    const Size dims = points.get_size1();

    Vec1Int        distribution = mlmpi::getBlockDistributionIndexes(points.get_size0(),
                                                              mlmpi->get_numberRanks(),
                                                              mlmpi->get_rank());
    Size           localSize = distribution.size();
    std::vector<T> localMinDistances(localSize, std::numeric_limits<T>::max());

    for (Size i = 0; i < numSamples; i++)
    {
        sampledIndices.push_back(nextFarthestIdx);

        const T* currentFarthestSlice = points.get_slice(nextFarthestIdx);
        T        maxMinDistLocal = -1.0;
        Int      nextFarthestIdxLocal = -1;

        for (Size k = 0; k < localSize; k++)
        {
            const Size globalIdx = distribution[k];
            if (globalIdx == (Size)nextFarthestIdx) localMinDistances[k] = -1.0;
            if (localMinDistances[k] >= 0.0)
            {
                T dist = computeDist(points.get_slice(globalIdx), currentFarthestSlice, dims);
                if (dist < localMinDistances[k]) localMinDistances[k] = dist;
                if (localMinDistances[k] > maxMinDistLocal)
                {
                    maxMinDistLocal = localMinDistances[k];
                    nextFarthestIdxLocal = static_cast<Int>(globalIdx);
                }
            }
        }
        // Global Sync
        MaxCandidate maxCandidateLocal = {maxMinDistLocal, nextFarthestIdxLocal};
        MaxCandidate maxCandidate = {-1.0, -1};
        MPI_Allreduce(&maxCandidateLocal,
                      &maxCandidate,
                      1,
                      MPI_DOUBLE_INT,
                      MPI_MAXLOC,
                      mlmpi->get_communicator());
        nextFarthestIdx = maxCandidate.nextFarthestIdx;
        if (nextFarthestIdx == -1) break;
    }

    return sampledIndices;
}
/*******************************************************************************************
 * @brief Performs farthest point sampling on a set of points with a distance threshold.
 *
 * This function selects a subset of points from the input set such that the selected points
 * are as far apart as possible, subject to a minimum distance threshold. The distance between
 * points is computed using a user-provided distance function.
 *
 * @tparam T The type of the elements in the points (e.g., float, double). Must be an arithmetic
 *type.
 * @tparam NormFunction A callable type that computes the distance between two points.
 *
 * @param points A 2D vector representing the input points. Each inner vector represents a point.
 * @param distanceThreshold The minimum distance threshold for selecting points.
 * @param computeDist A callable object (e.g., lambda, function pointer) that computes the distance
 *                    between two points. It should take two points as input and return a distance
 *of type T.
 *
 * @return A vector of indices representing the selected points. The indices correspond to the
 *positions of the selected points in the input `points` vector. Returns an empty vector if the
 *input is empty or no points satisfy the distance threshold.
 *
 * @warning This function only supports arithmetic types for the template parameter `T`. If a
 *non-arithmetic type is used, a runtime error will be triggered in debug mode.
 *******************************************************************************************/
template<typename T, typename NormFunction>
Vec1Int farthestPointSamplingDistTreshOMP(const std::vector<std::vector<T>>& points,
                                          const T                            distanceThreshold,
                                          NormFunction                       computeDist,
                                          const Size                         seed = 0)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );
    if (points.empty()) return {};
    Size nPoints = points.size();

    Int     firstIndex = chooseFirstIndexOMP(points.size(), seed);
    Vec1Int sampledIndices;
    sampledIndices.push_back(firstIndex);
    std::vector<T> minDistances(points.size(), std::numeric_limits<T>::max());
    minDistances[firstIndex] = -1.0;
    // 2. Iteratively select points
    while (true)
    {
        T   maxMinDist = -1.0;
        Int nextIndex = -1;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            T   localMaxDist = -1.0;
            Int localIndex = -1;
#ifdef _OPENMP
#pragma omp for nowait
#endif
            for (Size i = 0; i < nPoints; i++)
            {
                if (minDistances[i] < 0.0) continue;
                T minDistI = std::numeric_limits<T>::max();
                for (const Int selIdx : sampledIndices)
                {
                    T dist = computeDist(points[i], points[selIdx]);
                    minDistI = std::min(minDistI, dist);
                }
                if (minDistI > localMaxDist)
                {
                    localMaxDist = minDistI;
                    localIndex = static_cast<Int>(i);
                }
            }
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (localMaxDist > maxMinDist)
                {
                    maxMinDist = localMaxDist;
                    nextIndex = localIndex;
                }
            }
        }
        if (nextIndex != -1 && maxMinDist > distanceThreshold)
        {
            sampledIndices.push_back(nextIndex);
            minDistances[nextIndex] = -1.0;
        }
        else
        {
            break; // Stop when the distance gain is too small
        }
    }

    return sampledIndices;
}
/*******************************************************************************************
 * @brief Performs farthest point sampling with a distance threshold using MPI for distributed
 *computation.
 *
 * This function selects a subset of points from the input dataset such that the selected points are
 * as far apart as possible, subject to a distance threshold. The computation is distributed across
 * multiple MPI ranks for scalability.
 *
 * @tparam T Type of the data points. Must be an arithmetic type.
 * @tparam NormFunction Callable type for computing the distance between two points.
 *
 * @param points           A 2D array of points to sample from. Each row represents a point.
 * @param distanceThreshold The minimum distance threshold. Sampling stops when no point is farther
 *                          than this threshold from the already selected points.
 * @param computeDist       A callable object that computes the distance between two points.
 *                          It should have the signature:
 *                          `T computeDist(const T* pointA, const T* pointB, Size dims)`.
 * @param mlmpi             Shared pointer to an MPI wrapper object for distributed computation.
 *                          Defaults to `nullptr`.
 * @param seed              Seed for random number generation. If `0`, a random seed is generated
 *                          on rank 0 and broadcasted to other ranks. Defaults to `0`.
 *
 * @return A vector of indices representing the sampled points.
 *
 * @note This function assumes that the input points are distributed across MPI ranks. The
 *computation involves global synchronization using MPI to determine the farthest point at each
 *iteration. If the input dataset is empty, an empty vector is returned.
 *
 * @warning Only arithmetic types are supported for the template parameter `T`. If a non-arithmetic
 *          type is used, the function will trigger a runtime error.
 *
 * @details
 * The algorithm proceeds as follows:
 * 1. Randomly select the first point.
 * 2. Iteratively select the farthest point from the already selected points, ensuring that the
 *    distance is greater than the specified threshold.
 * 3. Use MPI to synchronize and determine the globally farthest point at each iteration.
 * 4. Stop when no point is farther than the threshold.
 *
 * The distance computation is delegated to the user-provided `computeDist` function, allowing
 * flexibility in defining the distance metric.
 *******************************************************************************************/
template<typename T, typename NormFunction>
Vec1Int farthestPointSamplingDistTreshMPI(const ShmemArray2D<T>&              points,
                                          const T                             distanceThreshold,
                                          NormFunction                        computeDist,
                                          const std::shared_ptr<const MlMPI>& mlmpi = nullptr,
                                          const Size                          seed = 0)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_arithmetic_v<T>)
        {
            global::tutor.bug("ERROR:" + flf(VASPML_FLF) + " only supports arithmetic types.");
        }
    );

    // handling base cases
    if (points.get_size0() == 0) { return {}; }
    // Stores the indices of the sampled points
    Vec1Int    sampledIndices;
    Int        firstIndex = chooseFirstIndexMPI(points.get_size0(), seed, mlmpi);
    Int        nextFarthestIdx = firstIndex;
    const Size dims = points.get_size1();

    Vec1Int distribution = mlmpi::getBlockDistributionIndexes(points.get_size0(),
                                                              mlmpi->get_numberRanks(),
                                                              mlmpi->get_rank());
    Size    localSize = distribution.size();

    std::vector<T> localMinDistances(localSize, std::numeric_limits<T>::max());
    // 2. Iteratively select points
    while (true)
    {
        sampledIndices.push_back(nextFarthestIdx);

        const T* currentFarthestSlice = points.get_slice(nextFarthestIdx);
        T        maxMinDistLocal = -1.0;
        Int      nextFarthestIdxLocal = -1;

        for (Size k = 0; k < localSize; k++)
        {
            const Size globalIdx = distribution[k];

            // If point is already selected, mark it so it's ignored
            if (globalIdx == (Size)nextFarthestIdx) localMinDistances[k] = -1.0;

            if (localMinDistances[k] >= 0.0)
            {
                T dist = computeDist(points.get_slice(globalIdx), currentFarthestSlice, dims);
                if (dist < localMinDistances[k]) localMinDistances[k] = dist;

                if (localMinDistances[k] > maxMinDistLocal)
                {
                    maxMinDistLocal = localMinDistances[k];
                    nextFarthestIdxLocal = static_cast<Int>(globalIdx);
                }
            }
        }

        // Global Sync to find the max of the local minimum distances
        MaxCandidate maxCandidateLocal = {static_cast<double>(maxMinDistLocal),
                                          nextFarthestIdxLocal};
        MaxCandidate maxCandidate = {-1.0, -1};

        MPI_Allreduce(&maxCandidateLocal,
                      &maxCandidate,
                      1,
                      MPI_DOUBLE_INT,
                      MPI_MAXLOC,
                      mlmpi->get_communicator());

        // Threshold Logic: Stop if no point is further away than the threshold
        if (maxCandidate.nextFarthestIdx != -1
            && maxCandidate.maxMinDist > static_cast<double>(distanceThreshold))
        {
            nextFarthestIdx = maxCandidate.nextFarthestIdx;
        }
        else { break; }
    }

    return sampledIndices;
}

} // namespace vaspml::fps

#endif
