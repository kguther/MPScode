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
  lapack_complex_double *targetPtr;
  lapack_complex_double result;
  Rctr.global_access((*MPOperator).length()-1,0,0,0)=1.0;
  for(int i=(*MPOperator).length()-2;i>=-1;--i){
    if(i==-1){
      targetPtr=&result;
    }
    else{
      Rctr.subContractionStart(targetPtr,i);
    }
    calcCtrIterRightBase(i,targetPtr);
  }
  lambda=real(result);
}


