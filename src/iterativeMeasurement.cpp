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

void iterativeMeasurement::calcCtrIterLeftBase(int const i, lapack_complex_double *targetPctr){
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
  calcOuterContainerLeft(i,outercontainer);
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

void iterativeMeasurement::calcOuterContainerLeft(int const i, tmpContainer<lapack_complex_double> &outercontainer){
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
