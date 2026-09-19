/* provide some timing and resource usage data (CPU, wall clock, memory, ...) */

#include <sys/resource.h>  // for getrusage(), struct rusage
#include <time.h>          // for clock_gettime(), struct timespec


void timing_C(int *mode, double *utime, double *stime, double *now,
              int *minpgf, int *majpgf, double *maxrsize, double *avsize,
              int *swaps, int *ios, int *cswitch, int *ierr)
{
    struct rusage rudata;
    struct timespec realtime;
    double intsize, totaltime;
    int dumerr;

    if (*mode == 0) {
        *ierr = getrusage(RUSAGE_SELF, &rudata);
    } else {
        *ierr = getrusage(RUSAGE_CHILDREN, &rudata);
    }

    *utime = (double)rudata.ru_utime.tv_sec + (double)rudata.ru_utime.tv_usec / 1.0e6;
    *stime = (double)rudata.ru_stime.tv_sec + (double)rudata.ru_stime.tv_usec / 1.0e6;

    dumerr = clock_gettime(CLOCK_REALTIME, &realtime);
    (void)dumerr;  // silence unused variable warning if not used

    *now = (double) realtime.tv_sec + (double) realtime.tv_nsec / 1.0e9;

    totaltime = *utime + *stime;
    if (totaltime <= 0.0) totaltime = 1.0;  // prevent division by zero

    *minpgf   = rudata.ru_minflt;
    *majpgf   = rudata.ru_majflt;

    intsize = (double)rudata.ru_ixrss + (double)rudata.ru_idrss + (double)rudata.ru_isrss;
    *maxrsize = (double)rudata.ru_maxrss;
    *avsize   = intsize / totaltime / 100.0;

    *swaps    = rudata.ru_nswap;
    *ios      = rudata.ru_inblock + rudata.ru_oublock;
    *cswitch  = rudata.ru_nvcsw;
}

