#ifndef MPO_MEASURE
#define MPO_MEASURE

#include "baseMeasurement.h"

//---------------------------------------------------------------------------------------------------//
// The global measurement can be used to compute the expectation value of some operator MPOperator
// when in some state MPState. It is quite straightforward.
//---------------------------------------------------------------------------------------------------//

class globalMeasurement: public baseMeasurement{
 public:
  globalMeasurement();
  globalMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  double measureFull();
};

#endif
