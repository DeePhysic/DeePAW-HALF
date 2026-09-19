/*****************************************************************
      FORTRAN interface to get user-time and real time on
      AIX-system might work on other UNIX systems as well
      (must work on all BSDish systems and SYSVish systems
      with BSD-compatibility timing routines ...)!?
*****************************************************************/

#include <sys/resource.h>  // getrusage()
#include <time.h>          // clock_gettime()

void vtime_C(double *vputim, double *cputim)
{
    struct rusage ppt, cpt;
    struct timespec now;
    int ierr;

    // Get process user and system time
    ierr = getrusage(RUSAGE_SELF, &ppt);
    // Get child process user and system time
    ierr = getrusage(RUSAGE_CHILDREN, &cpt);

    double tpu = (double)ppt.ru_utime.tv_sec + (double)ppt.ru_utime.tv_usec / 1e6;
    double tps = (double)ppt.ru_stime.tv_sec + (double)ppt.ru_stime.tv_usec / 1e6;
    double tcu = (double)cpt.ru_utime.tv_sec + (double)cpt.ru_utime.tv_usec / 1e6;
    double tcs = (double)cpt.ru_stime.tv_sec + (double)cpt.ru_stime.tv_usec / 1e6;

    *vputim = tpu + tps + tcu + tcs;

    ierr = clock_gettime(CLOCK_MONOTONIC, &now);
    *cputim = (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

void vtime2_C(double *cputim)
{
    struct timespec now;
    int ierr = clock_gettime(CLOCK_MONOTONIC, &now);
    *cputim = (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

/* Get current epoch from syscall */
long get_current_epoch()
{
    struct timespec current_time;
    long error = clock_gettime(CLOCK_REALTIME, &current_time);
    if (error != 0)
    {
        return -8;
    }
    return (long)current_time.tv_sec;
}

