/*****************************************************************
      FORTRAN interface fuer ERRF
*****************************************************************/

#include <math.h>


double errf_C(const double *x) {
    return erf(*x);
}

double errfc_C(const double *x) {
    return erfc(*x);
}

