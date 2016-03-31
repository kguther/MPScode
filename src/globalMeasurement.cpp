#include <iostream>
#include "globalMeasurement.h"

globalMeasurement::globalMeasurement():
  baseMeasurement()
{}

//---------------------------------------------------------------------------------------------------//

globalMeasurement::globalMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::setupMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  baseMeasurement::setupMeasurement(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//
// The global measurement class does not have a full pContraction object to cache intermediate results
// since this is not necessary. We only compute one full contraction of the network. This is done
// iteratively, using the results of the last step in the current step.
//---------------------------------------------------------------------------------------------------//

void globalMeasurement::measureFull(double &lambda){
  lapack_complex_double *targetPctr=new lapack_complex_double[MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim()];
  targetPctr[0]=1.0;
  for(int i=1;i<=MPOperator->length();++i){
    calcCtrIterLeftBase(i,targetPctr,targetPctr);
  }
  lambda=real(targetPctr[0]);
  delete[] targetPctr;
}


