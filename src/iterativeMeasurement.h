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
  iterativeMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  void initialize(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  int calcCtrFull(int const direction);
  void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  void calcOuterContainerRight(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  pContraction<lapack_complex_double> Rctr;
 protected:
  void calcCtrIterLeftQNOpt(int const i, lapack_complex_double *targetPctr);
  void calcOuterContainerLeftQNOpt(int const i, tmpContainer<lapack_complex_double> &outerContainer);
};

#endif
