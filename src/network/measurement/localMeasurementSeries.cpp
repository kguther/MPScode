#include "localMeasurementSeries.h"
#include "optHMatrix.h"
#include <memory>

localMeasurementSeries::localMeasurementSeries(localMpo<mpsEntryType > *const MPOperator, mps *const MPState):
  iterativeMeasurement(MPOperator,MPState),
  localMPOperator(MPOperator)
{}

//---------------------------------------------------------------------------------------------------//

void localMeasurementSeries::measureFull(std::vector<mpsEntryType > &lambda){
  int const L=MPOperator->length();
  mpsEntryType  result;
  int const operatorSize=localMPOperator->width();
  MPOperator->setUpSparse();
  //The input operator is stored since the measure sweep destroys its form
  localMpo<mpsEntryType > backup=*localMPOperator;
  calcCtrFull(1);
  lambda.clear();
  calcCtrIterRightBase(-1,&result);
  lambda.push_back(real(result));
  Lctr.global_access(0,0,0,0)=1.0;
  for(int i=1;i<=localMPOperator->currentSite();++i){
    calcCtrIterLeft(i);
  }
  //Beginning from the initial site, we sweep to the right and compute the expectation value on each site, using the unchanged partial contractions as intermediate results.
  for(int i=localMPOperator->currentSite();i<L-1;++i){
    localMPOperator->stepRight();
    for(int m=operatorSize-1;m>=0;--m){
      calcCtrIterLeft(i+1-m);
    }
    getCurrentValue(lambda,localMPOperator->currentSite());
    
  }
  //Here, the input operator is restored to its original form
  *localMPOperator=backup;
}

//---------------------------------------------------------------------------------------------------//

void localMeasurementSeries::getCurrentValue(std::vector<mpsEntryType > &lambda, int const i){
  mpsEntryType simpleContainer=0.0;
  mpsEntryType *LTerm, *RTerm, *currentM;
  getLocalDimensions(i);
  Rctr.subContractionStart(RTerm,i);
  Lctr.subContractionStart(LTerm,i);
  MPState->subMatrixStart(currentM,i);
  //Only the current site has to be contracted explicitly, the rest is stored in the partial contraction. Explicit contraction is carried out using the optHMatrix class multiplication.
  std::unique_ptr<mpsEntryType[]> siteMatrixContainerP(new mpsEntryType [ld*lDR*lDL]);
  mpsEntryType *siteMatrixContainer=siteMatrixContainerP.get();
  optHMatrix gather(RTerm,LTerm,MPOperator,MPState->getDimInfo(),MPOperator->maxDim(),i,0,0,0);
  gather.MultMv(currentM,siteMatrixContainer);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer+=conj(currentM[stateIndex(si,ai,aim)])*siteMatrixContainer[stateIndex(si,ai,aim)];
      }
    }
  }
  lambda.push_back((simpleContainer));
}
