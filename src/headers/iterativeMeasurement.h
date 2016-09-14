#ifndef ITERATIVE_MEASUREMENT
#define ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"
#include "templates/pContraction.h"

//---------------------------------------------------------------------------------------------------//
// This measurement class is used for computing the expectation value of some operator during the 
// optimization. It therefore comes with a left-side partial contraction in addition to the default
// right-side one. Also, the iterative calculations automatically use the next subarray of the partial
// contraction as target.
//---------------------------------------------------------------------------------------------------//

class iterativeMeasurement: protected baseMeasurement{
 public:
  iterativeMeasurement();
  iterativeMeasurement(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
  void initialize(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
  int calcCtrFull(int const direction);
  virtual void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  virtual void calcCtrIterRightBase(int i, mpsEntryType *const target);
  void calcOuterContainerLeft(int const i, tmpContainer<mpsEntryType > &outerContainer);
  void calcOuterContainerRight(int const i, tmpContainer<mpsEntryType > &outerContainer);
  //Containers for caching of partial contractions of the expectation value. This is what distinguishes the iterativeMeasurement
  pContraction<mpsEntryType > Lctr;
  pContraction<mpsEntryType > Rctr;
};

#endif
