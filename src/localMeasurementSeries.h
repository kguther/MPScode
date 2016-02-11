#ifndef LOCAL_MEASUREMENT_SERIES
#define LOCAL_MEASUREMENT_SERIES

#include <vector>
#include "iterativeMeasurement.h"
#include "localMpo.h"

class localMeasurementSeries: public iterativeMeasurement{
 public:
  localMeasurementSeries(localMpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  void measureFull(std::vector<lapack_complex_double> &lambda);
 private:
  localMpo<lapack_complex_double> *localMPOperator;
  void getCurrentValue(std::vector<lapack_complex_double> &lambda, int const i);
};

#endif
