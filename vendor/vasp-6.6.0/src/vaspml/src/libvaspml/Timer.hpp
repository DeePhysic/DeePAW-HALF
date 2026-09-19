#ifndef TIMER_HPP
#define TIMER_HPP

#include "types.hpp"

#include <chrono>

#if PROFILING
#define VASPML_PROFILING_START(keyword) do { \
      global::timer.start(keyword); \
   } while(0)
#else
#define VASPML_PROFILING_START(X)
#endif

#if PROFILING
#define VASPML_PROFILING_STOP(keyword) do { \
      global::timer.stop(keyword); \
   } while(0)
#else
#define VASPML_PROFILING_STOP(X)
#endif

#if PROFILING
#define VASPML_PROFILING_WRITE() do { \
      global::timer.writeToScreen(); \
   } while(0)
#else
#define VASPML_PROFILING_WRITE()
#endif

namespace vaspml
{

/*******************************************************************************************
 * @class Timer
 * @brief implements a class which can be used for code timings
 *
 * This class can be used to time code segements. The timings will be averaged over all
 * calls of the function calls and will also write the absolute time.
 *
 * Usage example:\n
 * @code
 * Timer time; \n
 * time.start("MyTimerID"); \n
 * some code segment which will be timed \n
 * time.stop("MyTimerID"); \n
 * time.writeToScreen(); \n
 * @endcode
 *******************************************************************************************/
class Timer
{
  public:
    explicit Timer(ShRec record = nullptr);
    /*******************************************************************************************
     * Start timing with given ID.
     *
     * @param id Representative name for the timer, e.g., function name that is timed.
     *******************************************************************************************/
    void start(const String& id);
    /*******************************************************************************************
     * Stop timing with given ID.
     *
     * @param id Representative name for the timer, e.g., function name that is timed.
     *******************************************************************************************/
    void stop(const String& id);
    /*******************************************************************************************
     * writing the timing data to screen.
     *
     * The write to screen function will print the total time for evert set timer, the average
     * time for every set timer and the factor between the current timer and the fastest timer
     *******************************************************************************************/
    void writeToScreen();

  private:
    /**********************************************************************************************
     * Record with memory which may be set up externally.
     **********************************************************************************************/
    ShRec data;
    /**********************************************************************************************
     * Accumulated time in seconds.
     **********************************************************************************************/
    Vec1String& ids;
    /**********************************************************************************************
     * Accumulated time in seconds.
     **********************************************************************************************/
    Vec1Real& tSum;
    /**********************************************************************************************
     * Minimum time in seconds.
     **********************************************************************************************/
    Vec1Real& tMin;
    /**********************************************************************************************
     * Maximum time in seconds.
     **********************************************************************************************/
    Vec1Real& tMax;
    /**********************************************************************************************
     * Number of timer starts.
     **********************************************************************************************/
    Vec1Real& nStarts;
    /**********************************************************************************************
     * Number of timer starts.
     **********************************************************************************************/
    Vec1Int& status;
    /**********************************************************************************************
     * Start time point (not saved in data record).
     **********************************************************************************************/
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> tStart;
    /**********************************************************************************************
     * Stop time point (not saved in data record).
     **********************************************************************************************/
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> tStop;
};

/**************************************************************************************************
 * Create data map with (key, type)-pairs for setting up ::data member.
 **************************************************************************************************/
MapString dataMapTimer();

namespace global
{
inline Timer timer;
}

} //namespace vaspml

#endif
