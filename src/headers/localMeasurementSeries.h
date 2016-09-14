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
  localMeasurementSeries(localMpo<mpsEntryType > *const MPOperator, mps *const MPState);
  void measureFull(std::vector<mpsEntryType > &lambda);
 private:
  localMpo<mpsEntryType > *localMPOperator;
  void getCurrentValue(std::vector<mpsEntryType > &lambda, int const i);
};

#endif
