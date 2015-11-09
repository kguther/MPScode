#ifndef BASE_MEASUREMENT
#define BASE_MEASUREMENT

#include <lapacke.h>
#include "mpo.h"
#include "mps.h"
#include "pContraction.h"

class baseMeasurement{
 public: 
  void calcCtrIterRightBase(int const i, lapack_complex_double *target);
  pContraction<lapack_complex_double> Rctr;
 protected:
  baseMeasurement();
  baseMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  mpo<lapack_complex_double> *MPOperator;
  mps *MPState;
  void getLocalDimensions(int const i);
  void initializeBase(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  int lDwL, lDwR, lDL, lDR, ld, D, Dw;
  int pctrIndex(int const ai, int const bi, int const aip){return aip+bi*D+ai*D*Dw;}
};

#endif
