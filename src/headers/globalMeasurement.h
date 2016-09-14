#ifndef MPO_MEASURE
#define MPO_MEASURE

#include "baseMeasurement.h"

//---------------------------------------------------------------------------------------------------//
// The global measurement can be used to compute the expectation value of some operator MPOperator
// when in some state MPState. It is quite straightforward. The main applciation is to check for
// fixed QNs and to get the variance of energy.
//---------------------------------------------------------------------------------------------------//

class globalMeasurement: public baseMeasurement{
 public:
  globalMeasurement();
  globalMeasurement(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
  void setupMeasurement(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
  void measureFull(double &lambda);
};

#endif
