#ifndef UNCACHED_ITERATIVE_MEASUREMENT
#define UNCACHED_ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"
#include <vector>

class uncachedMeasurement: public baseMeasurement{
 public:
  uncachedMeasurement();
  uncachedMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  void update();
  void getLctr(lapack_complex_double *&target){target=&(Lctr[0]);}
  void getRctr(lapack_complex_double *&target){target=&(Rctr[0]);}
 private:
  std::vector<lapack_complex_double> Rctr, Lctr;
};

#endif
