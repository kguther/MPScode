#include "iterativeMeasurement.h"

iterativeMeasurement::iterativeMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

iterativeMeasurement::iterativeMeasurement(mpo<mpsEntryType > *const MPOperatorIn, mps *const MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{
  Rctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
  Lctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::initialize(mpo<mpsEntryType > *const MPOperatorIn, mps *const MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase();
  Rctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
  Lctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state.
//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterLeft(int const i){
  mpsEntryType *targetPctr, *sourcePctr;
  Lctr.subContractionStart(sourcePctr,i-1);
  Lctr.subContractionStart(targetPctr,i);
  calcCtrIterLeftBase(i,sourcePctr,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterRight(int const i){
  mpsEntryType *targetPctr;
  Rctr.subContractionStart(targetPctr,i);
  calcCtrIterRightBase(i,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterRightBase(int i, mpsEntryType  *const targetPctr){
  //only required to get expectation values - this uses the right-sided contraction
  mpsEntryType *sourcePctr;
  Rctr.subContractionStart(sourcePctr,i+1);
  baseMeasurement::calcCtrIterRightBase(i,sourcePctr,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcOuterContainerLeft(int const i, tmpContainer<mpsEntryType > &outercontainer){
  mpsEntryType *sourcePctr;
  Lctr.subContractionStart(sourcePctr,i-1);
  baseMeasurement::calcOuterContainerLeft(i,sourcePctr,outercontainer);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcOuterContainerRight(int const i, tmpContainer<mpsEntryType > &outercontainer){
  mpsEntryType *sourcePctr;
  Rctr.subContractionStart(sourcePctr,i+1);
  baseMeasurement::calcOuterContainerRight(i,sourcePctr,outercontainer);
}

//---------------------------------------------------------------------------------------------------//

int iterativeMeasurement::calcCtrFull(int const direction){
  //Full calculation of the contraction is only required once: before the first sweep
  //This is just some ordinary iterative computation of the partial contraction Pctr (P=R,L)
  int L=MPOperator->length();
  if(direction==1){
    Rctr.global_access(L-1,0,0,0)=1.0;
    for(int i=L-2;i>=0;--i){
      calcCtrIterRight(i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Lctr.global_access(0,0,0,0)=1.0;
      for(int i=1;i<L;++i){
	calcCtrIterLeft(i);
      }
      return 0;
    }
    else{
      return -1;
    }
  }
}
