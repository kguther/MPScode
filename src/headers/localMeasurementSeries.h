#ifndef LOCAL_MEASUREMENT_SERIES
#define LOCAL_MEASUREMENT_SERIES

#include <vector>
#include "iterativeMeasurement.h"
#include "templates/localMpo.h"

//---------------------------------------------------------------------------------------------------//
// This important class computes the expectation value of some site-dependant ('local') MPO for all 
// sites right to its current site. The operator is not changed.
//---------------------------------------------------------------------------------------------------//

class localMeasurementSeries: private iterativeMeasurement{
 public:
  localMeasurementSeries(localMpo<std::complex<double> > *const MPOperator, mps *const MPState);
  void measureFull(std::vector<std::complex<double> > &lambda);
 private:
  localMpo<std::complex<double> > *localMPOperator;
  void getCurrentValue(std::vector<std::complex<double> > &lambda, int const i);
};

#endif
