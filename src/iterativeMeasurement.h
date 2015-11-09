#ifndef ITERATIVE_MEASUREMENT
#define ITERATIVE_MEASUREMENT

#include <lapacke.h>
#include "mpo.h"
#include "mps.h"
#include "pContraction.h"
#include "baseMeasurement.h"

class iterativeMeasurement: public baseMeasurement{
 public:
  iterativeMeasurement();
  void initialize(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  int calcCtrFull(int const direction);
  void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  pContraction<lapack_complex_double> Lctr;
};

#endif
