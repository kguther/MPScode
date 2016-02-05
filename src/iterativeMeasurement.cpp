#include <iostream>
#include <time.h>
#include "iterativeMeasurement.h"

iterativeMeasurement::iterativeMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

iterativeMeasurement::iterativeMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  baseMeasurement(MPOperatorIn,MPStateIn)
{
  Rctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::initialize(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase(MPOperatorIn,MPStateIn);
  Rctr.initialize(MPOperator->length(),MPState->maxDim(),MPOperator->maxDim());
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state.
//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterLeft(int const i){
  lapack_complex_double *targetPctr;
  Lctr.subContractionStart(targetPctr,i);
  calcCtrIterLeftBase(i,targetPctr);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterRight(int const i){
  lapack_complex_double *target;
  Rctr.subContractionStart(target,i);
  calcCtrIterRightBase(i,target);
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterRightBase(int const i, lapack_complex_double *targetPctr){
  if(MPState->indexTable.nQNs()){
    calcCtrIterRightBaseQNOpt(i,targetPctr);
  }
  else{
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

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcOuterContainerRight(int const i, tmpContainer<lapack_complex_double> &outercontainer){
  if(MPState->indexTable.nQNs()){
    calcOuterContainerRightQNOpt(i,outercontainer);
  }
  else{
    //If one wanted to use the code without exploiting QNs, this had to be added
  }
}

void iterativeMeasurement::calcCtrIterRightBaseQNOpt(int const i, lapack_complex_double *targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *siteMatrixState;
  int const numBlocks=MPState->indexTable.numBlocksRP(i+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  MPState->subMatrixStart(siteMatrixState,i+1);
  getLocalDimensions(i+1);
  tmpContainer<lapack_complex_double> outercontainer(lDL,lDwL,ld,lDR);
  calcOuterContainerRightQNOpt(i,outercontainer);

  /*
  lapack_complex_double *sourcePctr;
  lapack_complex_double *siteMatrixH;
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=MPOperator->numEls(i+1);
  int biS, bimS, siS, sipS;
  Rctr.subContractionStart(sourcePctr,i+1);
  MPOperator->sparseSubMatrixStart(siteMatrixH,i+1);
  MPOperator->biSubIndexArrayStart(biIndices,i+1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i+1);
  MPOperator->siSubIndexArrayStart(siIndices,i+1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
  curtime=clock();
  for(int sip=0;sip<ld;++sip){                                  
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  innercontainer.global_access(sip,bi,ai,aimp)=0;
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeRP(i+1,iBlock);
    rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
      siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimB=MPState->indexTable.aimBlockIndexRP(i+1,iBlock,k);
	for(int bi=0;bi<lDwR;++bi){
	  for(int ai=0;ai<lDR;++ai){
	    innercontainer.global_access(siB,bi,ai,aimB)+=sourcePctr[pctrIndex(ai,bi,aiB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
    for(int j=0;j<rBlockSize;++j){
      siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
      aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
      for(int aim=0;aim<lDL;++aim){
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
  */
  
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeRP(i+1,iBlock);
    rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
    for(int k=0;k<lBlockSize;++k){
      aimB=MPState->indexTable.aimBlockIndexRP(i+1,iBlock,k);
      for(int bim=0;bim<lDwL;++bim){
	for(int aimp=0;aimp<lDL;++aimp){
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

void iterativeMeasurement::calcOuterContainerRightQNOpt(int const i, tmpContainer<lapack_complex_double> &outercontainer){
  lapack_complex_double *siteMatrixState;
  int const numBlocks=MPState->indexTable.numBlocksRP(i+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  MPState->subMatrixStart(siteMatrixState,i+1);
  getLocalDimensions(i+1);
 
  lapack_complex_double *sourcePctr;
  lapack_complex_double *siteMatrixH;
  int *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=MPOperator->numEls(i+1);
  int biS, bimS, siS, sipS;
  Rctr.subContractionStart(sourcePctr,i+1);
  MPOperator->sparseSubMatrixStart(siteMatrixH,i+1);
  MPOperator->biSubIndexArrayStart(biIndices,i+1);
  MPOperator->bimSubIndexArrayStart(bimIndices,i+1);
  MPOperator->siSubIndexArrayStart(siIndices,i+1);
  MPOperator->sipSubIndexArrayStart(sipIndices,i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
  for(int sip=0;sip<ld;++sip){                                  
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  innercontainer.global_access(sip,bi,ai,aimp)=0;
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=MPState->indexTable.lBlockSizeRP(i+1,iBlock);
    rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
      siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimB=MPState->indexTable.aimBlockIndexRP(i+1,iBlock,k);
	for(int bi=0;bi<lDwR;++bi){
	  for(int ai=0;ai<lDR;++ai){
	    innercontainer.global_access(siB,bi,ai,aimB)+=sourcePctr[pctrIndex(ai,bi,aiB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    rBlockSize=MPState->indexTable.rBlockSizeRP(i+1,iBlock);
    for(int j=0;j<rBlockSize;++j){
      siB=MPState->indexTable.siBlockIndexRP(i+1,iBlock,j);
      aiB=MPState->indexTable.aiBlockIndexRP(i+1,iBlock,j);
      for(int aim=0;aim<lDL;++aim){
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
