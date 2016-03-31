#ifndef ITERATIVE_MEASUREMENT
#define ITERATIVE_MEASUREMENT

#include "baseMeasurement.h"

//---------------------------------------------------------------------------------------------------//
// This measurement class is used for computing the expectation value of some operator during the 
// optimization. It therefore comes with a left-side partial contraction in addition to the default
// right-side one. Also, the iterative calculations automatically use the next subarray of the partial
// contraction as target.
//---------------------------------------------------------------------------------------------------//

class iterativeMeasurement: public baseMeasurement{
 public:
  iterativeMeasurement();
  iterativeMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  void initialize(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  int calcCtrFull(int const direction);
  void calcCtrIterRightBase(int const i, lapack_complex_double *target);
  void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  void calcOuterContainerLeft(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  void calcOuterContainerRight(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  //Containers for caching of partial contractions of the expectation value. This is what distinguishes the iterativeMeasurement
  pContraction<lapack_complex_double> Lctr;
  pContraction<lapack_complex_double> Rctr;
 protected:
  void calcCtrIterRightBaseQNOpt(int const i, lapack_complex_double *target);
  void calcOuterContainerRightQNOpt(int const i, tmpContainer<lapack_complex_double> &outerContainer);
};

#endif
