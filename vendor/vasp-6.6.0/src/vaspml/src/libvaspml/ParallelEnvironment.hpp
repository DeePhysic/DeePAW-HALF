#ifndef PARALLELENVIRONMENT_HPP
#define PARALLELENVIRONMENT_HPP

#include "debug.hpp"
#include "types.hpp"

#if VASPML_DEBUG_LEVEL > 2
#include <iostream>
#endif

// Macro definitions inside the VASPML_DOXYGEN block are only for documentation purposes.
#ifdef VASPML_DOXYGEN
/**************************************************************************************************
 * Enables multithreading via C++ Parallel Algorithms code path.
 **************************************************************************************************/
#define VASPML_PALGO_THREADS
/**************************************************************************************************
 * Enables offloading to GPU via C++ Parallel Algorithms code path.
 **************************************************************************************************/
#define VASPML_PALGO_GPU
#endif // VASPML_DOXYGEN

/*================================================================================================+
 |
 | SECTION Preprocessor macros, compiler-specific settings.
 |
 +================================================================================================*/

/*================================================================================================+
 | Invert default for using std::execution.
 +================================================================================================*/
#ifndef VASPML_USE_STDEXECUTION
#define VASPML_NO_STDEXECUTION
#endif

/*================================================================================================+
 | Compile without any usage of std::execution and other selected C++17 features.
 |
 | This build should avoid most compiler troubles related to C++17 features: no std::execution
 | policies are used for C++ algorithms like std::for_each(). Typically, the functions without
 | execution policy argument were present since C++98 or C++11. However, a few were introduced
 | together with policies in C++17 and are also disabled with this build:
 |
 | * std::reduce()
 | * std::transform_reduce()
 |
 | Explanation
 | -----------
 | Those cause further problems in combination with policies and old system compilers: some
 | compilers depend on an underlying system compiler and its headers. Specifically, Intel's oneAPI
 | and NVIDIA's HPC SDK may rely on a potentially old system GCC. Then, GCC headers may not contain
 | std::reduce() and std::transform_reduce():
 |
 | * with    execution policy first argument: present since GCC 9.1.0
 | * without execution policy first argument: present since GCC 9.3.0
 |
 | The std::(transform_)reduce() function WITH execution policy is typically defined in a header
 | from Intel or NVIDIA but for the version WITHOUT policy they are relying on the system compiler's
 | headers. However, those will not contain the declarations if GCC is older than 9.3.0 and cause
 | compilation errors. Theoretically, in an effort to keep std::(transform_)reduce() we would have
 | to actually USE execution policies just for std::(transform_)reduce() while removing them for all
 | other functions. To avoid this weirdness we do not use std::(transform_)reduce() when C++17
 | execution policies are disabled.
 +================================================================================================*/
#ifdef VASPML_NO_STDEXECUTION

#define VASPML_POLICY
#define VASPML_POLICY_SEQ
#define VASPML_POLICY_PAR
#define VASPML_POLICY_PAR_UNSEQ
#define VASPML_EXEC_SPACE_SPECIFIER

/*================================================================================================+
 | Standard build, assume existence of std::execution.
 +================================================================================================*/
#else // VASPML_NO_STDEXECUTION

#define VASPML_POLICY policy,
#define VASPML_POLICY_SEQ seq,
#define VASPML_POLICY_PAR par,
#define VASPML_POLICY_PAR_UNSEQ par_unseq,

/*================================================================================================+
 | NVIDIA-specific implementation details.
 |
 | Available compiler flags:
 | * -stdpar=multicore (uses OpenMP threads via Thrust)
 | * -stdpar=gpu
 | * -stdpar=gpu:acc
 | Combinations of "multicore" and "gpu(:acc)" are allowed: -stdpar=multicore,gpu.
 |
 | If GPUs are used we also need to prepare the use of cuBLAS (include headers, prepare handle to
 | library context).
 +================================================================================================*/
#if defined(__NVCOMPILER)

#include <execution>

#ifdef __NVCOMPILER_STDPAR_MULTICORE
#define VASPML_PALGO_THREADS
#endif

#ifdef __NVCOMPILER_STDPAR_GPU
// clang-format off
#include <cuda_runtime.h>
#include <cublas_v2.h>
// clang-format on
#define VASPML_PALGO_GPU "nvidia"
#define VASPML_PALGO_GPU_NVIDIA
#define VASPML_EXEC_SPACE_SPECIFIER __host__ __device__
#endif

#ifdef __NVCOMPILER_STDPAR_OPENACC_GPU
// clang-format off
#include <cuda_runtime.h>
#include <cublas_v2.h>
// clang-format on
#define VASPML_PALGO_GPU "nvidia-acc"
#define VASPML_PALGO_GPU_NVIDIA
#define VASPML_EXEC_SPACE_SPECIFIER __host__ __device__
#endif

#ifndef VASPML_EXEC_SPACE_SPECIFIER
#define VASPML_EXEC_SPACE_SPECIFIER
#endif

#ifdef VASPML_PALGO_GPU_NVIDIA
#if VASPML_DEBUG_LEVEL > 2
#include <nvtx3/nvtx3.hpp>
#endif
#endif

namespace vaspml
{
constexpr auto seq = std::execution::seq;
constexpr auto par = std::execution::par;
constexpr auto par_unseq = std::execution::par_unseq;
} //namespace vaspml

/*================================================================================================+
 | Intel oneAPI-specific implementation details.
 |
 | TODO: This does not work yet, Intel requires using their special headers like
 | <oneapi/dpl/algorithm> instead of standard headers (<algorithm>) every time an execution policy
 | is passed.
 +================================================================================================*/
//#elif defined(__INTEL_LLVM_COMPILER)
//
//// clang-format off
//#include <oneapi/dpl/execution>
//#include <execution>
//// clang-format on
//namespace vaspml
//{
//constexpr auto seq = oneapi::dpl::execution::seq;
//constexpr auto par = oneapi::dpl::execution::par;
//constexpr auto par_unseq = oneapi::dpl::execution::par_unseq;
//} //namespace vaspml
//#define VASPML_EXEC_SPACE_SPECIFIER

/*================================================================================================+
 | Compiler-agnostic setup.
 +================================================================================================*/
#else

#include <execution> // IWYU pragma: export
namespace vaspml
{
constexpr auto seq = std::execution::seq;
constexpr auto par = std::execution::par;
constexpr auto par_unseq = std::execution::par_unseq;
} //namespace vaspml
#define VASPML_PALGO_THREADS
#define VASPML_EXEC_SPACE_SPECIFIER

#endif // Compiler-specific implementation details.

#endif // VASPML_NO_STDEXECUTION

/*================================================================================================+
 |
 | SECTION ParallelEnvironment class.
 |
 +================================================================================================*/

namespace vaspml
{

#ifdef VASPML_NO_STDEXECUTION
class DummyPolicy
{};
#endif

class ParallelEnvironment
{
  public:
    explicit ParallelEnvironment(ShRec record = nullptr);
    ~ParallelEnvironment();
    /**********************************************************************************************
     * Determine requested parallelization mode from string (one of palgos).
     **********************************************************************************************/
    void init(const String palgoRequested);
    /**********************************************************************************************
     * Call given lambda with execution policy corresponding to selected parallelization mode.
     **********************************************************************************************/
    template<class F>
    auto run(F f, [[maybe_unused]] String hint = "") const
    {
#ifdef VASPML_NO_STDEXECUTION
        DummyPolicy dummyPolicy;
        return f(dummyPolicy);
#else
        bool parallel = true;
        if (palgo == "off" || palgo == "serial" || hint == "seq") parallel = false;
        VASPML_DEBUG_L3(
            std::cout << "PALGO (" << palgo << ") call: " << hint << "\n";
        );
// NVTX scoped range is active until it goes out of scope, therefore it cannot be wrapped by
// VASPML_DEBUG_L3(). 
#if VASPML_DEBUG_LEVEL > 2
#ifdef VASPML_PALGO_GPU_NVIDIA
        nvtx3::scoped_range r{hint};
#endif
#endif
        if (parallel) { return f(par_unseq); }
        else { return f(seq); }
#endif
    }
    /**********************************************************************************************
     * Check if init() was ever called (strictly required only if some handle needs to be set up).
     **********************************************************************************************/
    bool isInitialized() const;
    /**********************************************************************************************
     * Check if currently selected parallelization mode does not use C++ Parallel Algorithms.
     **********************************************************************************************/
    bool off() const;
    /**********************************************************************************************
     * Check if currently selected parallelization mode makes use of GPU.
     **********************************************************************************************/
    bool gpu() const;
    /**********************************************************************************************
     * Get currently selected parallelization algorithm.
     **********************************************************************************************/
    String selected() const;
    /**********************************************************************************************
     * Get comma-separated list of supported parallelization algorithms.
     **********************************************************************************************/
    String supported() const;
#ifdef VASPML_PALGO_GPU_NVIDIA
    /**********************************************************************************************
     * Convert status enum to string.
     **********************************************************************************************/
    String convert(cublasStatus_t status) const;
    /**********************************************************************************************
     * Getter for handle to cuBLAS library context.
     **********************************************************************************************/
    const cublasHandle_t& get_handle() const;
    /**********************************************************************************************
     * Destroy cuBLAS handle.
     **********************************************************************************************/
    void destroyHandle();
#endif // VASPML_PALGO_GPU_NVIDIA

  private:
    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    bool  initialized;
    /**********************************************************************************************
     * Actually selected parallelization options derived from ML_PALGO.
     **********************************************************************************************/
    String& palgo;
    /**********************************************************************************************
     * Further specification of parallelization option (if necessary).
     **********************************************************************************************/
    String& flavor;
    /**********************************************************************************************
     * List of available parallelization options (depends on build options).
     **********************************************************************************************/
    Vec1String& palgos;
#ifdef VASPML_PALGO_GPU_NVIDIA
    /**********************************************************************************************
     * Handle to cuBLAS library context.
     **********************************************************************************************/
    cublasHandle_t handle;
#endif // VASPML_PALGO_GPU_NVIDIA
};

/**************************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **************************************************************************************************/
MapString dataMapParallelEnvironment();

namespace global
{
/**************************************************************************************************
 * Global instance of ParallelEnvironment class, makes selected parallelization available
 * everywhere.
 **************************************************************************************************/
inline ParallelEnvironment parallel;
} //namespace global

} //namespace vaspml

#endif // PARALLEL_ENVIRONMENT_HPP
