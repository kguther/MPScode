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
  Lctr[0]=1.0;
  Rctr[0]=1.0;
}

//---------------------------------------------------------------------------------------------------//
// Functions for updating the L/R-terms
//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::update(){
  int i=(MPState->length())/2;
  getLeftCtr(i);
  getRightCtr(i);
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getContractions(int site){
  int const i=site+1;
  for(int j=1;j<=i;++j){
    getLeftCtr(j);
  }
  for(int j=MPState->dimInfo.L()-1;j>=i;++j){
    getRightCtr(j);
  }
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getLeftCtr(int i){
  calcCtrIterLeftBase(i,&(Lctr[0]),&(Lctr[0]));
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getRightCtr(int i){
  calcCtrIterRightBase(i-1,&(Rctr[0]),&(Rctr[0]));
}
