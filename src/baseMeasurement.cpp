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

baseMeasurement::baseMeasurement(mpo<lapack_complex_double> *const MPOperatorIn, mps *const MPStateIn):
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
}

//---------------------------------------------------------------------------------------------------//  

void baseMeasurement::setupMeasurement(mpo<lapack_complex_double> *const MPOperatorIn, mps *const MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase();
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

void baseMeasurement::calcCtrIterLeftBase(int const i, lapack_complex_double *const source, lapack_complex_double *const targetPctr){
  calcCtrIterLeftBaseQNOpt(i,source,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcOuterContainerLeft(int const i, lapack_complex_double *const source, tmpContainer<lapack_complex_double> &outercontainer){
  if(MPState->indexTable.nQNs()){
    calcOuterContainerLeftQNOpt(i,source,outercontainer);
  }
  else{
    //If one wanted to use the code without exploiting QNs, this had to be added
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterLeftBaseQNOpt(int const i, lapack_complex_double *const source, lapack_complex_double *const targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *siteMatrixState;
  int const numBlocks=MPState->indexTable.numBlocksLP(i-1);
  int aiB, aimB, siB;
  int lBlockSize, rBlockSize;
  clock_t curtime;
  MPState->subMatrixStart(siteMatrixState,i-1);
  getLocalDimensions(i-1);
  //container arrays to significantly reduce computational effort by storing intermediate results
  tmpContainer<lapack_complex_double> outercontainer(ld,lDwR,lDR,lDL);
  calcOuterContainerLeftQNOpt(i,source,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aiB,siB,aimB)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
      rBlockSize=MPState->indexTable.rBlockSizeLP(i-1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=MPState->indexTable.aiBlockIndexLP(i-1,iBlock,j);
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0;
	  for(int k=0;k<lBlockSize;++k){
	    siB=MPState->indexTable.siBlockIndexLP(i-1,iBlock,k);
	    aimB=MPState->indexTable.aimBlockIndexLP(i-1,iBlock,k);
	    simpleContainer+=conj(siteMatrixState[stateIndex(siB,aiB,aimB)])*outercontainer.global_access(siB,bi,aip,aimB);
	  }
	  targetPctr[pctrIndex(aiB,bi,aip)]=simpleContainer;
	}
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Total left contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
}

void baseMeasurement::calcOuterContainerLeftQNOpt(int const i, lapack_complex_double *const source, tmpContainer<lapack_complex_double> &outercontainer){
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  lapack_complex_double *siteMatrixH, *siteMatrixState;
  int biS, bimS, siS, sipS, aiB, aimB, siB;
  int lBlockSize, rBlockSize;
  int const sparseSize=MPOperator->numEls(i-1);
  int const numBlocks=MPState->indexTable.numBlocksLP(i-1);
  clock_t curtime;
  MPState->subMatrixStart(siteMatrixState,i-1);
  getLocalDimensions(i-1);
  MPOperator->sparseSubMatrixStart(siteMatrixH,i-1);
  MPOperator->biSubIndexArrayStart(biIndices,i-1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i-1);
  MPOperator->siSubIndexArrayStart(siIndices,i-1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i-1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDR,lDwL,lDL);
  curtime=clock();
  //horrible construct to efficiently compute the partial contraction
#pragma omp parallel for
  for(int sip=0;sip<ld;++sip){
    for(int bim=0;bim<lDwL;++bim){
      for(int aim=0;aim<lDL;++aim){
	for(int aip=0;aip<lDR;++aip){
	  innercontainer.global_access(sip,aip,bim,aim)=0;
	}
      }
    }
  }
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,aimB,siB)
  for(int bim=0;bim<lDwL;++bim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
      rBlockSize=MPState->indexTable.rBlockSizeLP(i-1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=MPState->indexTable.aimBlockIndexLP(i-1,iBlock,k);
	siB=MPState->indexTable.siBlockIndexLP(i-1,iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=MPState->indexTable.aiBlockIndexLP(i-1,iBlock,j);
	  for(int aim=0;aim<lDL;++aim){
	    innercontainer.global_access(siB,aiB,bim,aim)+=source[pctrIndex(aim,bim,aimB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
#pragma omp parallel for private(lBlockSize,rBlockSize,siB,aimB,siS,biS,bimS,sipS)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	siB=MPState->indexTable.siBlockIndexLP(i-1,iBlock,k);
	aimB=MPState->indexTable.aimBlockIndexLP(i-1,iBlock,k);
	for(int bi=0;bi<lDwR;++bi){
	  outercontainer.global_access(siB,bi,aip,aimB)=0;
	}
	for(int nSparse=0;nSparse<sparseSize;++nSparse){
	  siS=siIndices[nSparse];
	  if(siS==siB){
	    biS=biIndices[nSparse];
	    bimS=bimIndices[nSparse];
	    sipS=sipIndices[nSparse];
	    outercontainer.global_access(siB,biS,aip,aimB)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,aip,bimS,aimB);
	  }
	}
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterRightBaseQNOpt(int const i, lapack_complex_double *const sourcePctr, lapack_complex_double *const targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *siteMatrixState;
  int const numBlocks=MPState->indexTable.numBlocksRP(i+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  MPState->subMatrixStart(siteMatrixState,i+1);
  getLocalDimensions(i+1);
  tmpContainer<lapack_complex_double> outercontainer(lDL,lDwL,ld,lDR);
  //The calculation of the first two contractions has to be done in other functions, too. It therefore has an extra function. 
  calcOuterContainerRightQNOpt(i,sourcePctr,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aimB,aiB,siB)  
  for(int aimp=0;aimp<lDL;++aimp){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=MPState->indexTable.lBlockSizeRP(i+1,iBlock);
      rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=MPState->indexTable.aimBlockIndexRP(i+1,iBlock,k);
	for(int bim=0;bim<lDwL;++bim){
	  simpleContainer=0;
	  for(int j=0;j<rBlockSize;++j){
	    siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
	    aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
	    simpleContainer+=conj(siteMatrixState[stateIndex(siB,aiB,aimB)])*outercontainer.global_access(aimp,bim,siB,aiB);
	  }
	  targetPctr[pctrIndex(aimB,bim,aimp)]=simpleContainer;
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcOuterContainerRightQNOpt(int const i, lapack_complex_double *const sourcePctr, tmpContainer<lapack_complex_double> &outercontainer){
  lapack_complex_double *siteMatrixState;
  int const numBlocks=MPState->indexTable.numBlocksRP(i+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  MPState->subMatrixStart(siteMatrixState,i+1);
  getLocalDimensions(i+1);
  lapack_complex_double *siteMatrixH;
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=MPOperator->numEls(i+1);
  int biS, bimS, siS, sipS;
  MPOperator->sparseSubMatrixStart(siteMatrixH,i+1);
  MPOperator->biSubIndexArrayStart(biIndices,i+1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i+1);
  MPOperator->siSubIndexArrayStart(siIndices,i+1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
#pragma omp parallel for
  for(int sip=0;sip<ld;++sip){                                  
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  innercontainer.global_access(sip,bi,ai,aimp)=0;
	}
      }
    }
  }
  //Contraction works similar to the general case. Here, only entries fulfilling the QN constraint are used
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,siB,aimB)
  for(int ai=0;ai<lDR;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=MPState->indexTable.lBlockSizeRP(i+1,iBlock);
      rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
	siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimB=MPState->indexTable.aimBlockIndexRP(i+1,iBlock,k);
	  for(int bi=0;bi<lDwR;++bi){
	    innercontainer.global_access(siB,bi,ai,aimB)+=sourcePctr[pctrIndex(ai,bi,aiB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,siB,siS,sipS,biS,bimS)
  for(int aim=0;aim<lDL;++aim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
	aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
	for(int bim=0;bim<lDwL;++bim){
	  outercontainer.global_access(aim,bim,siB,aiB)=0;
	}
	for(int nSparse=0;nSparse<sparseSize;++nSparse){
	  siS=siIndices[nSparse];
	  if(siS==siB){
	    biS=biIndices[nSparse];
	    bimS=bimIndices[nSparse];
	    sipS=sipIndices[nSparse];
	    outercontainer.global_access(aim,bimS,siB,aiB)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,biS,aiB,aim);
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterRightBase(int const i, lapack_complex_double *const sourcePctr, lapack_complex_double *const targetPctr){
  if(MPState->indexTable.nQNs()){
    //Only this is usually relevant
    calcCtrIterRightBaseQNOpt(i,sourcePctr,targetPctr);
  }
  else{
    lapack_complex_double simpleContainer;
    lapack_complex_double *siteMatrixState, *siteMatrixH;
    int *biIndices, *bimIndices, *siIndices, *sipIndices;
    int const sparseSize=MPOperator->numEls(i+1);
    int biS, bimS, siS, sipS;
    clock_t curtime;
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
  }
}
