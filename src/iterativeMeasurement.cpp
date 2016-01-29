#include <iostream>
#include <time.h>
#include "iterativeMeasurement.h"

iterativeMeasurement::iterativeMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

iterativeMeasurement::iterativeMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{
  Lctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::initialize(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase(MPOperatorIn,MPStateIn);
  Lctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state.
//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterLeft(int const i){
  lapack_complex_double *targetPctr;
  Lctr.subContractionStart(targetPctr,i);
  calcCtrIterLeft(i,targetPctr);
}

void iterativeMeasurement::calcCtrIterLeft(int const i, lapack_complex_double *targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr;
  lapack_complex_double *siteMatrixState, *siteMatrixH;
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=MPOperator->numEls(i-1);
  int biS, bimS, siS, sipS;
  clock_t curtime;
  MPState->subMatrixStart(siteMatrixState,i-1);
  MPOperator->sparseSubMatrixStart(siteMatrixH,i-1);
  MPOperator->biSubIndexArrayStart(biIndices,i-1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i-1);
  MPOperator->siSubIndexArrayStart(siIndices,i-1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i-1);
  Lctr.subContractionStart(sourcePctr,i-1);
  getLocalDimensions(i-1);
  //container arrays to significantly reduce computational effort by storing intermediate results
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwL,lDL,lDR);
  tmpContainer<lapack_complex_double> outercontainer(ld,lDwR,lDR,lDL);
  curtime=clock();
  //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
  for(int sip=0;sip<ld;++sip){
    for(int bim=0;bim<lDwL;++bim){
      for(int aim=0;aim<lDL;++aim){
	for(int aip=0;aip<lDR;++aip){
	  simpleContainer=0;
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=sourcePctr[pctrIndex(aim,bim,aimp)]*siteMatrixState[stateIndex(sip,aip,aimp)];
	  }
	  innercontainer.global_access(sip,bim,aim,aip)=simpleContainer;
	}
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
  for(int aim=0;aim<lDL;++aim){
    for(int aip=0;aip<lDR;++aip){
      for(int si=0;si<ld;++si){
	for(int bi=0;bi<lDwR;++bi){
	      outercontainer.global_access(si,bi,aip,aim)=0;
	    }
	  }
      for(int nSparse=0;nSparse<sparseSize;++nSparse){
	biS=biIndices[nSparse];
	bimS=bimIndices[nSparse];
	siS=siIndices[nSparse];
	sipS=sipIndices[nSparse];
	outercontainer.global_access(siS,biS,aip,aim)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,bimS,aim,aip);
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
  for(int ai=0;ai<lDR;++ai){
    for(int bi=0;bi<lDwR;++bi){
      for(int aip=0;aip<lDR;++aip){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int aim=0;aim<lDL;++aim){
	    simpleContainer+=conj(siteMatrixState[stateIndex(si,ai,aim)])*outercontainer.global_access(si,bi,aip,aim);
	  }
	}
	targetPctr[pctrIndex(ai,bi,aip)]=simpleContainer;
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Total left contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterRight(int const i){
  lapack_complex_double *target;
  Rctr.subContractionStart(target,i);
  calcCtrIterRightBase(i,target);
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
      std::cout<<"CRITICAL ERROR: Invalid sweep direction identifier in calculation of partial contractions\n";
      return -1;
    }
  }
}
