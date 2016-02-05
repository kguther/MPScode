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
  Lctr.initialize((*MPOperator).length(),(*MPState).maxDim(),(*MPOperator).maxDim());
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

void baseMeasurement::calcCtrIterLeftBase(int const i, lapack_complex_double *targetPctr){
  calcCtrIterLeftBaseQNOpt(i,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcOuterContainerLeft(int const i, tmpContainer<lapack_complex_double> &outercontainer){
  if(MPState->indexTable.nQNs()){
    calcOuterContainerLeftQNOpt(i,outercontainer);
  }
  else{
    //If one wanted to use the code without exploiting QNs, this had to be added
  }
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterLeftBaseQNOpt(int const i, lapack_complex_double *targetPctr){
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
  calcOuterContainerLeftQNOpt(i,outercontainer);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
    rBlockSize=MPState->indexTable.rBlockSizeLP(i-1,iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiB=MPState->indexTable.aiBlockIndexLP(i-1,iBlock,j);
      for(int bi=0;bi<lDwR;++bi){
	for(int aip=0;aip<lDR;++aip){
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

void baseMeasurement::calcOuterContainerLeftQNOpt(int const i, tmpContainer<lapack_complex_double> &outercontainer){
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  lapack_complex_double *sourcePctr, *siteMatrixH, *siteMatrixState;
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
  Lctr.subContractionStart(sourcePctr,i-1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDR,lDwL,lDL);
  curtime=clock();
  //horrible construct to efficiently compute the partial contraction
  for(int sip=0;sip<ld;++sip){
    for(int bim=0;bim<lDwL;++bim){
      for(int aim=0;aim<lDL;++aim){
	for(int aip=0;aip<lDR;++aip){
	  innercontainer.global_access(sip,aip,bim,aim)=0;
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
    rBlockSize=MPState->indexTable.rBlockSizeLP(i-1,iBlock);
    for(int k=0;k<lBlockSize;++k){
      aimB=MPState->indexTable.aimBlockIndexLP(i-1,iBlock,k);
      siB=MPState->indexTable.siBlockIndexLP(i-1,iBlock,k);
      for(int j=0;j<rBlockSize;++j){
	aiB=MPState->indexTable.aiBlockIndexLP(i-1,iBlock,j);
	for(int bim=0;bim<lDwL;++bim){
	  for(int aim=0;aim<lDL;++aim){
	    innercontainer.global_access(siB,aiB,bim,aim)+=sourcePctr[pctrIndex(aim,bim,aimB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
  curtime=clock()-curtime;
  //std::cout<<"Inner contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  curtime=clock();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeLP(i-1,iBlock);
    for(int k=0;k<lBlockSize;++k){
      siB=MPState->indexTable.siBlockIndexLP(i-1,iBlock,k);
      aimB=MPState->indexTable.aimBlockIndexLP(i-1,iBlock,k);
      for(int aip=0;aip<lDR;++aip){
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
