#ifndef UNCACHED_ITERATIVE_MEASUREMENT
#define UNCACHED_ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"
#include "impBase.h"
#include "contractor.h"
#include <vector>

//---------------------------------------------------------------------------------------------------//
// A variant of the iterative measurement without a global cache. Only the last update of the partial
// contractions is stored, drastically reducing memory usage. Applicable in iDMRG, where no global
// cache is needed.
//---------------------------------------------------------------------------------------------------//


class uncachedMeasurement{
 public:
  uncachedMeasurement();
  uncachedMeasurement(mpo<arcomplex<double> > *const MPOperatorIn, impBase *const MPStateIn);
  void update();
  void getLctr(arcomplex<double>  *&target){target=&(Lctr[0]);}
  void getRctr(arcomplex<double>  *&target){target=&(Rctr[0]);}
  void setContractions(std::vector<arcomplex<double> > const &R, std::vector<arcomplex<double> > const &L);
  void getLeftCtr();
  void getRightCtr();
 private:
  impBase *MPState;
  mpo<arcomplex<double> > *MPOperator;
  contractor calcer;
  std::vector<arcomplex<double> > Rctr, Lctr;
};

#endif
