#include "localMeasurementSeries.h"
#include "optHMatrix.h"

localMeasurementSeries::localMeasurementSeries(localMpo<lapack_complex_double> *const MPOperator, mps *const MPState):
  iterativeMeasurement(MPOperator,MPState),
  localMPOperator(MPOperator)
{}

void localMeasurementSeries::measureFull(std::vector<lapack_complex_double> &lambda){
  int const L=MPOperator->length();
  lapack_complex_double result;
  MPOperator->setUpSparse();
  localMpo<lapack_complex_double> backup=*localMPOperator;
  calcCtrFull(1);
  lambda.clear();
  calcCtrIterRightBase(-1,&result);
  lambda.push_back(real(result));
  Lctr.global_access(0,0,0,0)=1.0;
  for(int i=1;i<=localMPOperator->currentSite();++i){
    calcCtrIterLeft(i);
  }
  for(int i=localMPOperator->currentSite();i<L-1;++i){
    localMPOperator->stepRight();
    calcCtrIterLeft(i+1);
    getCurrentValue(lambda,localMPOperator->currentSite());
  }
  *localMPOperator=backup;
}

void localMeasurementSeries::getCurrentValue(std::vector<lapack_complex_double> &lambda, int const i){
  lapack_complex_double simpleContainer=0.0;
  lapack_complex_double *LTerm, *RTerm, *currentM;
  getLocalDimensions(i);
  Rctr.subContractionStart(RTerm,i);
  Lctr.subContractionStart(LTerm,i);
  MPState->subMatrixStart(currentM,i);
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
