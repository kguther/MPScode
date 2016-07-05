#include "baseMeasurement.h"
#include "mpo.h"
#include "mps.h"
#include "tmpContainer.h"
#include "contractor.h"

//---------------------------------------------------------------------------------------------------//
// Constructor and initialize() function for the baseMeasurement class.
//---------------------------------------------------------------------------------------------------//

baseMeasurement::baseMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

baseMeasurement::baseMeasurement(mpo<arcomplex<double> > *const MPOperatorIn, mps *const MPStateIn):
  MPOperator(MPOperatorIn),
  MPState(MPStateIn)
{
  initializeBase();  
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::initializeBase(){
  Dw=MPOperator->maxDim();
  D=MPState->maxDim();
  MPOperator->setUpSparse();
  calcer=contractor(MPOperator->maxDim(),MPState->getDimInfo());
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::setupMeasurement(mpo<arcomplex<double> > *const MPOperatorIn, mps *const MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase();
}

//---------------------------------------------------------------------------------------------------//  
// This does the same as the getLocalDimensions of the network
//---------------------------------------------------------------------------------------------------//

void baseMeasurement::getLocalDimensions(int i){
  lDL=MPState->locDimL(i);
  lDR=MPState->locDimR(i);
  ld=MPState->locd(i);
  lDwR=MPOperator->locDimR(i);
  lDwL=MPOperator->locDimL(i);
}

//---------------------------------------------------------------------------------------------------//
// Here we compute the partial contraction of the network consisting of the mps and the mpo from
// the right up to site i from that up to site i+1. The result is stored in the targetPctr array,
// which will usually be the next subcontraction of Rctr.
//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterLeftBase(int const i, arcomplex<double>  *const sourcePctr, arcomplex<double>  *const targetPctr){
  arcomplex<double> *siteMatrixState;
  MPState->subMatrixStart(siteMatrixState,i-1);
  if(MPState->indexTable().nQNs()){
      calcer.calcLeftContraction(i,MPState->indexTable().getLocalIndexTable(i-1),siteMatrixState,MPOperator->getSiteTensor(i-1),sourcePctr,targetPctr);
  }
  else{
    calcer.calcLeftContraction(i,siteMatrixState,MPOperator->getSiteTensor(i-1),sourcePctr,targetPctr);
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterRightBase(int const i, arcomplex<double>  *const sourcePctr, arcomplex<double>  *const targetPctr){
  arcomplex<double> *siteMatrixState;
  MPState->subMatrixStart(siteMatrixState,i+1);
  if(MPState->indexTable().nQNs()){
    //Only this is usually relevant
      calcer.calcRightContraction(i,MPState->indexTable().getLocalIndexTable(i+1),siteMatrixState,MPOperator->getSiteTensor(i+1),sourcePctr,targetPctr);
  }
  else{
    calcer.calcRightContraction(i,siteMatrixState,MPOperator->getSiteTensor(i+1),sourcePctr,targetPctr);
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcOuterContainerLeft(int const i, arcomplex<double>  *const source, tmpContainer<arcomplex<double> > &outercontainer){
  arcomplex<double> *siteMatrixState;
  MPState->subMatrixStart(siteMatrixState,i-1);
  if(MPState->indexTable().nQNs()){
    calcer.calcLeftOuterContainer(i,MPState->indexTable().getLocalIndexTable(i-1),siteMatrixState,MPOperator->getSiteTensor(i-1),source,outercontainer);
  }
  else{
    calcer.calcLeftOuterContainer(i,siteMatrixState,MPOperator->getSiteTensor(i-1),source,outercontainer);
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcOuterContainerRight(int const i, arcomplex<double>  *const source, tmpContainer<arcomplex<double> > &outercontainer){
  arcomplex<double> *siteMatrixState;
  MPState->subMatrixStart(siteMatrixState,i+1);
  if(MPState->indexTable().nQNs()){
    calcer.calcRightOuterContainer(i,MPState->indexTable().getLocalIndexTable(i+1),siteMatrixState,MPOperator->getSiteTensor(i+1),source,outercontainer);
  }
  else{
    calcer.calcRightOuterContainer(i,siteMatrixState,MPOperator->getSiteTensor(i+1),source,outercontainer);
  }
}
