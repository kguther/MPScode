#ifdef USE_MKL
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif
#include <omp.h>
#include <complex>
#include <mkl.h>
#endif
#ifndef USE_MKL
#include <lapacke.h>
#include <cblas.h>
#include <lapacke_utils.h>
#endif
