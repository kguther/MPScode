#include "localMeasurementSeries.h"
#include "optHMatrix.h"
#include "globalMeasurement.h"
#include <iostream>

localMeasurementSeries::localMeasurementSeries(localMpo<lapack_complex_double> *const MPOperator, mps *const MPState):
  iterativeMeasurement(MPOperator,MPState),
  localMPOperator(MPOperator)
{}

//---------------------------------------------------------------------------------------------------//

void localMeasurementSeries::measureFull(std::vector<lapack_complex_double> &lambda){
  int const L=MPOperator->length();
  lapack_complex_double result;
  int const operatorSize=localMPOperator->width();
  MPOperator->setUpSparse();
  //The input operator is stored since the measure sweep destroys its form
  localMpo<lapack_complex_double> backup=*localMPOperator;
  calcCtrFull(1);
  lambda.clear();
  calcCtrIterRightBase(-1,&result);
  lambda.push_back(real(result));
  Lctr.global_access(0,0,0,0)=1.0;
  for(int i=1;i<=localMPOperator->currentSite();++i){
    calcCtrIterLeft(i);
  }
  globalMeasurement test;
  double buffer;
  //Beginning from the initial site, we sweep to the right and compute the expectation value on each site, using the unchanged partial contractions as intermediate results.
  for(int i=localMPOperator->currentSite();i<L-1;++i){
    localMPOperator->stepRight();
    /*
    test.setupMeasurement(localMPOperator,MPState);
    test.measureFull(buffer);
    lambda.push_back(buffer);
    */
    for(int m=operatorSize-1;m>=0;--m){
      calcCtrIterLeft(i+1-m);
    }
    getCurrentValue(lambda,localMPOperator->currentSite());
  }
  //Here, the input operator is restored to its original form
  *localMPOperator=backup;
}

//---------------------------------------------------------------------------------------------------//

void localMeasurementSeries::getCurrentValue(std::vector<lapack_complex_double> &lambda, int const i){
  lapack_complex_double simpleContainer=0.0;
  lapack_complex_double *LTerm, *RTerm, *currentM;
  getLocalDimensions(i);
  Rctr.subContractionStart(RTerm,i);
  Lctr.subContractionStart(LTerm,i);
  MPState->subMatrixStart(currentM,i);
  //Only the current site has to be contracted explicitly, the rest is stored in the partial contraction. Explicit contraction is carried out using the optHMatrix class multiplication.
  lapack_complex_double *siteMatrixContainer=new lapack_complex_double [ld*lDR*lDL];
  optHMatrix gather(RTerm,LTerm,MPOperator,MPState->dimInfo,MPOperator->maxDim(),i,0,0,0);
  gather.MultMv(currentM,siteMatrixContainer);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer+=conj(currentM[stateIndex(si,ai,aim)])*siteMatrixContainer[stateIndex(si,ai,aim)];
      }
    }
  }
  lambda.push_back((simpleContainer));
  delete[] siteMatrixContainer;
}
