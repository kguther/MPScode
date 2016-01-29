#include <iostream>
#include "globalMeasurement.h"

globalMeasurement::globalMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  iterativeMeasurement(MPOperatorIn,MPStateIn)
{}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::setupMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  baseMeasurement::setupMeasurement(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::measureFull(double &lambda){
  calcCtrFull(-1);
  lapack_complex_double result;
  calcCtrIterLeft(MPOperator->length(),&result);
  lambda=real(result);
}


