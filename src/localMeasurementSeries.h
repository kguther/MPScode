#ifndef LOCAL_MEASUREMENT_SERIES
#define LOCAL_MEASUREMENT_SERIES

#include <vector>
#include "iterativeMeasurement.h"
#include "localMpo.h"

//---------------------------------------------------------------------------------------------------//
// This important class computes the expectation value of some site-dependant ('local') MPO for all 
// sites right to its current site. The operator is not changed.
//---------------------------------------------------------------------------------------------------//

class localMeasurementSeries: private iterativeMeasurement{
 public:
  localMeasurementSeries(localMpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  void measureFull(std::vector<lapack_complex_double> &lambda);
 private:
  localMpo<lapack_complex_double> *localMPOperator;
  void getCurrentValue(std::vector<lapack_complex_double> &lambda, int const i);
};

#endif
