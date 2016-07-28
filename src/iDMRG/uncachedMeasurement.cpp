#include "uncachedMeasurement.h"

uncachedMeasurement::uncachedMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

uncachedMeasurement::uncachedMeasurement(mpo<arcomplex<double> > *const MPOperatorIn, impBase *const MPStateIn):
  MPState(MPStateIn),
  MPOperator(MPOperatorIn),
  calcer(contractor(MPOperatorIn->maxDim(),MPStateIn->getDimInfo()))
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
  getRightCtr();
  getLeftCtr();
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getLeftCtr(){
  int iGlobal=MPState->currentSite();
  int const iH=(iGlobal==0)?0:1;
  arcomplex<double> *container;
  MPState->subMatrixStart(container,iGlobal);
  //Here, we have to differ between left- and right index table
  calcer.calcLeftContraction(iGlobal+1,container,MPOperator->getSiteTensor(iH),&(Lctr[0]),&(Lctr[0]));
}

//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::getRightCtr(){
  int iGlobal=MPState->currentSite();
  int const iH=(iGlobal==0)?(MPOperator->length()-1):1;
  arcomplex<double> *container;
  MPState->subMatrixStart(container,iGlobal+1);
  calcer.calcRightContraction(iGlobal,container,MPOperator->getSiteTensor(iH),&(Rctr[0]),&(Rctr[0]));
}

//---------------------------------------------------------------------------------------------------//
// Load Lctr and Rctr from some external source (necessary when starting with extended system)
//---------------------------------------------------------------------------------------------------//

void uncachedMeasurement::setContractions(std::vector<arcomplex<double> > const &R, std::vector<arcomplex<double> > const &L){
  Rctr=R;
  Lctr=L;
}


