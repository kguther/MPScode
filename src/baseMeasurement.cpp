#include "baseMeasurement.h"
#include "mpo.h"
#include "mps.h"
#include "tmpContainer.h"
#include "pContraction.h"

//---------------------------------------------------------------------------------------------------//
// Constructor and initialize() function for the baseMeasurement class.
//---------------------------------------------------------------------------------------------------//

baseMeasurement::baseMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

baseMeasurement::baseMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  MPOperator(MPOperatorIn),
  MPState(MPStateIn)
{
  initializeBase(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::initializeBase(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  Dw=(*MPOperator).maxDim();
  D=(*MPState).maxDim();
  Rctr.initialize((*MPOperator).length(),(*MPState).maxDim(),(*MPOperator).maxDim());
}

//---------------------------------------------------------------------------------------------------//  

void baseMeasurement::setupMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//  
// This does the same as the getLocalDimensions of the network
//---------------------------------------------------------------------------------------------------//

void baseMeasurement::getLocalDimensions(int const i){
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

void baseMeasurement::calcCtrIterRightBase(int const i, lapack_complex_double *targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr;
  lapack_complex_double *siteMatrixState, *siteMatrixH;
  Rctr.subContractionStart(sourcePctr,i+1);
  (*MPState).subMatrixStart(siteMatrixState,i+1);
  (*MPOperator).subMatrixStart(siteMatrixH,i+1);
  getLocalDimensions(i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
  tmpContainer<lapack_complex_double> outercontainer(ld,lDwL,lDL,lDR);
  for(int sip=0;sip<ld;++sip){                                                       
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=sourcePctr[pctrIndex(ai,bi,aip)]*siteMatrixState[stateIndex(sip,aip,aimp)];
	  }
	  innercontainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	for(int ai=0;ai<lDR;++ai){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=siteMatrixH[operatorIndex(si,sip,bi,bim)]*innercontainer.global_access(sip,bi,ai,aimp);
	    }
	  }
	  outercontainer.global_access(si,bim,aimp,ai)=simpleContainer;
	}
      }
    }
  }
  for(int aim=0;aim<lDL;++aim){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=conj(siteMatrixState[stateIndex(si,ai,aim)])*outercontainer.global_access(si,bim,aimp,ai);
	  }
	}
	targetPctr[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
}
