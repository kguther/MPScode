#include "globalMeasurement.h"
#include <memory>
#include "mkl_complex_defined.h"

globalMeasurement::globalMeasurement():
  baseMeasurement()
{}

//---------------------------------------------------------------------------------------------------//

globalMeasurement::globalMeasurement(mpo<mpsEntryType > *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{}

//---------------------------------------------------------------------------------------------------//

void globalMeasurement::setupMeasurement(mpo<mpsEntryType > *MPOperatorIn, mps *MPStateIn){
  baseMeasurement::setupMeasurement(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//
// The global measurement class does not have a full pContraction object to cache intermediate results
// since this is not necessary. We only compute one full contraction of the network. This is done
// iteratively, using the results of the last step in the current step.
//---------------------------------------------------------------------------------------------------//

void globalMeasurement::measureFull(double &lambda){
  std::unique_ptr<mpsEntryType >targetPctrP(new mpsEntryType[MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim()]);
  mpsEntryType *targetPctr=targetPctrP.get();
  targetPctr[0]=1.0;
  /*
  for(int i=MPOperator->length()-2;i>=-1;--i){
    calcCtrIterRightBase(i,targetPctr,targetPctr);
  }
  */
  
  for(int i=1;i<=MPOperator->length();++i){
    calcCtrIterLeftBase(i,targetPctr,targetPctr);
  }
  
  lambda=real(targetPctr[0]);
}


