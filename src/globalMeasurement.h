#ifndef MPO_MEASURE
#define MPO_MEASURE

#include <lapacke.h>
#include "mpo.h"
#include "mps.h"
#include "baseMeasurement.h"

class globalMeasurement: public baseMeasurement{
 public:
  globalMeasurement();
  globalMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  double measureFull();
};

#endif
