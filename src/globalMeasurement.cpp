#include <iostream>
#include "globalMeasurement.h"

globalMeasurement::globalMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{}

//---------------------------------------------------------------------------------------------------//

double globalMeasurement::measureFull(){
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
  return real(result);
}


