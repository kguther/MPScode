#ifndef MPO_MEASURE
#define MPO_MEASURE

#include "iterativeMeasurement.h"

//---------------------------------------------------------------------------------------------------//
// The global measurement can be used to compute the expectation value of some operator MPOperator
// when in some state MPState. It is quite straightforward.
//---------------------------------------------------------------------------------------------------//

class globalMeasurement: public iterativeMeasurement{
 public:
  globalMeasurement();
  globalMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  void setupMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  void measureFull(double &lambda);
};

#endif
