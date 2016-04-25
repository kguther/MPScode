#ifndef ITERATIVE_MEASUREMENT
#define ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"

//---------------------------------------------------------------------------------------------------//
// This measurement class is used for computing the expectation value of some operator during the 
// optimization. It therefore comes with a left-side partial contraction in addition to the default
// right-side one. Also, the iterative calculations automatically use the next subarray of the partial
// contraction as target.
//---------------------------------------------------------------------------------------------------//

class iterativeMeasurement: protected baseMeasurement{
 public:
  iterativeMeasurement();
  iterativeMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  void initialize(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  int calcCtrFull(int const direction);
  virtual void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  virtual void calcCtrIterRightBase(int i, lapack_complex_double *const target);
  void calcOuterContainerLeft(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  void calcOuterContainerRight(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  //Containers for caching of partial contractions of the expectation value. This is what distinguishes the iterativeMeasurement
  pContraction<lapack_complex_double> Lctr;
  pContraction<lapack_complex_double> Rctr;
};

#endif
