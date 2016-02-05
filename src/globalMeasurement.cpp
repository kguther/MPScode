#include <iostream>
#include "globalMeasurement.h"

globalMeasurement::globalMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::setupMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  baseMeasurement::setupMeasurement(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::measureFull(double &lambda){
  Lctr.global_access(0,0,0,0)=1.0;
  lapack_complex_double *targetPctr;
  for(int i=1;i<MPOperator->length();++i){
    Lctr.subContractionStart(targetPctr,i);
    calcCtrIterLeftBase(i,targetPctr);
  }
  lapack_complex_double result;
  calcCtrIterLeftBase(MPOperator->length(),&result);
  lambda=real(result);
}


