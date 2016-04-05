#include "uncachedMeasurement.h"

uncachedMeasurement::uncachedMeasurement():
  baseMeasurement()
{}

//---------------------------------------------------------------------------------------------------//

uncachedMeasurement::uncachedMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState):
  baseMeasurement(MPOperator,MPState)
{
  Lctr.resize(MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim());
  Rctr.resize(MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::update(){
  int i=(MPState->length())/2;
  calcCtrIterLeftBase(i+1,&(Lctr[0]),&(Lctr[0]));
  calcCtrIterRightBase(i-1,&(Rctr[0]),&(Rctr[0]));
}
