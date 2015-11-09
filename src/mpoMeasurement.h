#ifndef MPO_MEASURE
#define MPO_MEASURE

#include <lapacke.h>
#include "mpo.h"
#include "mps.h"
#include "baseMeasurement.h"

class mpoMeasurement: public baseMeasurement{
 public:
  mpoMeasurement();
  mpoMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  double measureFull();
};

#endif
