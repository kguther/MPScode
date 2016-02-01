#include "baseMeasurement.h"
#include "mpo.h"
#include "mps.h"
#include "tmpContainer.h"
#include "pContraction.h"
#include <iostream>
#include <time.h>

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
  Dw=MPOperator->maxDim();
  D=MPState->maxDim();
  Rctr.initialize((*MPOperator).length(),(*MPState).maxDim(),(*MPOperator).maxDim());
  MPOperator->setUpSparse();
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
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=MPOperator->numEls(i+1);
  int biS, bimS, siS, sipS;
  clock_t curtime;
  Rctr.subContractionStart(sourcePctr,i+1);
  MPState->subMatrixStart(siteMatrixState,i+1);
  MPOperator->sparseSubMatrixStart(siteMatrixH,i+1);
  MPOperator->biSubIndexArrayStart(biIndices,i+1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i+1);
  MPOperator->siSubIndexArrayStart(siIndices,i+1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i+1);
  getLocalDimensions(i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
  tmpContainer<lapack_complex_double> outercontainer(lDL,lDwL,ld,lDR);
  curtime=clock();
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
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
  for(int aimp=0;aimp<lDL;++aimp){
    for(int ai=0;ai<lDR;++ai){
      for(int si=0;si<ld;++si){
	for(int bim=0;bim<lDwL;++bim){
	  outercontainer.global_access(aimp,bim,si,ai)=0;
	}
      }
      for(int nSparse=0;nSparse<sparseSize;++nSparse){
	biS=biIndices[nSparse];
	bimS=bimIndices[nSparse];
	siS=siIndices[nSparse];
	sipS=sipIndices[nSparse];
	outercontainer.global_access(aimp,bimS,siS,ai)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,biS,ai,aimp);
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Outer contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
  for(int aim=0;aim<lDL;++aim){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=conj(siteMatrixState[stateIndex(si,ai,aim)])*outercontainer.global_access(aimp,bim,si,ai);
	  }
	}
	targetPctr[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Total contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
}
