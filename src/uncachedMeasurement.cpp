#include "uncachedMeasurement.h"

uncachedMeasurement::uncachedMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

uncachedMeasurement::uncachedMeasurement(mpo<arcomplex<double> > *const MPOperatorIn, impBase *const MPStateIn):
  MPState(MPStateIn),
  MPOperator(MPOperatorIn),
  calcer(contractor(MPOperatorIn->maxDim(),MPStateIn->getDimInfo(),MPStateIn->indexTable()))
{
  Lctr.resize(MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim());
  Rctr.resize(MPState->maxDim()*MPState->maxDim()*MPOperator->maxDim());
  Lctr[0]=1.0;
  Rctr[0]=1.0;
  MPOperator->setUpSparse();
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
  arcomplex<double> *container;
  MPState->subMatrixStart(container,iGlobal);
  calcer.calcLeftContraction(iGlobal+1,iInternal+1,container,MPOperator->getSiteTensor(iInternal),&(Lctr[0]),&(Lctr[0]));
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getRightCtr(){
  int iGlobal=MPState->currentSite();
  int iInternal=MPState->internalSite();
  arcomplex<double> *container;
  MPState->subMatrixStart(container,iGlobal+1);
  calcer.calcRightContraction(iGlobal,iInternal,container,MPOperator->getSiteTensor(iInternal+1),&(Rctr[0]),&(Rctr[0]));
}

//---------------------------------------------------------------------------------------------------//
// Load Lctr and Rctr from some external source (necessary when starting with extended system)
//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::setContractions(std::vector<arcomplex<double> > const &R, std::vector<arcomplex<double> > const &L){
  Rctr=R;
  Lctr=L;
}


