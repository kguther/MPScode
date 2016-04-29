#ifndef UNCACHED_ITERATIVE_MEASUREMENT
#define UNCACHED_ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"
#include "impBase.h"
#include <vector>

//---------------------------------------------------------------------------------------------------//
// A variant of the iterative measurement without a global cache. Only the last update of the partial
// contractions is stored, drastically reducing memory usage. Applicable in iDMRG, where no global
// cache is needed.
//---------------------------------------------------------------------------------------------------//


class uncachedMeasurement{
 public:
  uncachedMeasurement();
  uncachedMeasurement(mpo<lapack_complex_double> *const MPOperatorIn, impBase *const MPStateIn);
  void update();
  void getLctr(lapack_complex_double *&target){target=&(Lctr[0]);}
  void getRctr(lapack_complex_double *&target){target=&(Rctr[0]);}
  void getLeftCtr(int i);
  void getRightCtr(int i);
  void getContractions(int site);
 private:
  impBase *MPState;
  mpo<lapack_complex_double> *MPOperator;
  std::vector<lapack_complex_double> Rctr, Lctr;
};

#endif
