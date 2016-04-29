#include "uncachedMeasurement.h"

uncachedMeasurement::uncachedMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

uncachedMeasurement::uncachedMeasurement(mpo<lapack_complex_double> *const MPOperatorIn, impBase *const MPStateIn):
  MPState(MPStateIn),
  MPOperator(MPOperatorIn)
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
  getLeftCtr();
  getRightCtr();
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getLeftCtr(){
  int iGlobal=MPState->currentSite();
  int iInternal=MPState->internalSite();
  lapack_complex_double *container;
  //calcCtrIterLeftBase(i,&(Lctr[0]),&(Lctr[0]));
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getRightCtr(){
  //calcCtrIterRightBase(i-1,&(Rctr[0]),&(Rctr[0]));
}

//---------------------------------------------------------------------------------------------------//
// Load Lctr and Rctr from some external source (necessary when starting with extended system)
//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::setContractions(std::vector<lapack_complex_double> const &R, std::vector<lapack_complex_double> const &L){
  Rctr=R;
  Lctr=L;
}


