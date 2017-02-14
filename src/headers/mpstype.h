#ifndef TYPEDEF_FOR_MPS_ENTRIES
#define TYPEDEF_FOR_MPS_ENTRIES

#ifdef REAL_MPS_ENTRIES
typedef double mpsEntryType;
//some hack to formulate everything for complex entries
inline double conj(double x){
  return x;
}
inline double real(double x){
  return x;
}
#else
#include <complex>
typedef std::complex<double> mpsEntryType;
#endif

void nancheck(int dim, mpsEntryType *array);
#endif
