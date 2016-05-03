#ifndef CALCULATION_OF_TENSOR_CONTRACTIONS
#define CALCULATION_OF_TENSOR_CONTRACTIONS

#include "dimensionTable.h"
#include "tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// Core functionality of locally contracting some MPS given by siteMatrixState with some MPO given by siteMatrixH. Does not care how the MPS is build -> featured in both baseMeasurement and uncachedMeasurement.
//---------------------------------------------------------------------------------------------------//

class contractor{
 public:
  contractor(){}
 contractor(int DwIn, dimensionTable const &dimInfoIn, basisQNOrderMatrix const &indexTableIn):dimInfo(dimInfoIn), indexTable(&indexTableIn),Dw(DwIn){nQNs=indexTable->nQNs(); D=dimInfo.D();}

  template<typename T>
    void calcLeftContraction(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcRightContraction(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcLeftOuterContainer(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);
 template<typename T>
    void calcRightOuterContainer(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);

 private:
  dimensionTable dimInfo;
  int lDL, lDR, ld, lDwR, lDwL;
  basisQNOrderMatrix const *indexTable;
  int nQNs;
  int D, Dw;
  //Index functions
  int pctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*D*Dw;}
  int stateIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
  int operatorIndex(int const si, int const sip, int const bi, int const bim) {return bim+bi*Dw+sip*Dw*Dw+si*ld*Dw*Dw;}

 void getLocalDimensions(int i){
    lDL=dimInfo.locDimL(i);
    lDR=dimInfo.locDimR(i);
    ld=dimInfo.locd(i);
  }
  
  //Versions of contractions optimized for systems with conserved QNs (are called automatically when available)
  template<typename T>
    void calcLeftContractionQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcRightContractionQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcLeftOuterContainerQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);
  template<typename T>
    void calcRightOuterContainerQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);
};

//---------------------------------------------------------------------------------------------------//

//Currently, these are just forwarding to the QN-optimized versions. This is very efficient when symmetries are used but gives an overhead when they are disabled. For specialization for unsymmetric systems, reimplement the general versions.
template<typename T>
void contractor::calcLeftContraction(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  if(nQNs){
    calcLeftContractionQNOpt(i,internalSite,siteMatrixState,H,source,target);
  }
}

template<typename T>
void contractor::calcRightContraction(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  if(nQNs){
    calcRightContractionQNOpt(i,internalSite,siteMatrixState,H,source,target);
  }
}

template<typename T>
void contractor::calcLeftOuterContainer(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  if(nQNs){
    calcLeftOuterContainerQNOpt(i,internalSite,siteMatrixState,H,source,outerContainer);
  }
}

template<typename T>
void contractor::calcRightOuterContainer(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  if(nQNs){
    calcRightOuterContainerQNOpt(i,internalSite,siteMatrixState,H,source,outerContainer);
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcLeftContractionQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  int const numBlocks=indexTable->numBlocksLP(internalSite-1);
  int lBlockSize, rBlockSize;
  int aiB, siB, aimB;
  getLocalDimensions(i-1);
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  //container arrays to significantly reduce computational effort by storing intermediate results
  tmpContainer<T> outercontainer(ld,lDwR,lDR,lDL);
  calcLeftOuterContainerQNOpt(i,internalSite,siteMatrixState,H,source,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aiB,siB,aimB)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(internalSite-1,iBlock);
      rBlockSize=indexTable->rBlockSizeLP(internalSite-1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=indexTable->aiBlockIndexLP(i-1,iBlock,j);
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0;
	  for(int k=0;k<lBlockSize;++k){
	    siB=indexTable->siBlockIndexLP(internalSite-1,iBlock,k);
	    aimB=indexTable->aimBlockIndexLP(internalSite-1,iBlock,k);
	    simpleContainer+=conj(siteMatrixState[stateIndex(siB,aiB,aimB)])*outercontainer.global_access(siB,bi,aip,aimB);
	  }
	  target[pctrIndex(aiB,bi,aip)]=simpleContainer;
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcLeftOuterContainerQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  int const *biIndices, *bimIndices, *siIndices, *sipIndices;
  T const *siteMatrixH;
  int biS, bimS, siS, sipS, aiB, aimB, siB;
  int lBlockSize, rBlockSize;
  int const sparseSize=H.numEls();
  int const numBlocks=indexTable->numBlocksLP(internalSite-1);
  getLocalDimensions(i-1);
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  H.sparseSubMatrixStart(siteMatrixH);
  H.biSubIndexArrayStart(biIndices);
  H.bimSubIndexArrayStart(bimIndices);
  H.siSubIndexArrayStart(siIndices);
  H.sipSubIndexArrayStart(sipIndices);
  tmpContainer<T> innercontainer(ld,lDR,lDwL,lDL);
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
      lBlockSize=indexTable->lBlockSizeLP(internalSite-1,iBlock);
      rBlockSize=indexTable->rBlockSizeLP(internalSite-1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=indexTable->aimBlockIndexLP(internalSite-1,iBlock,k);
	siB=indexTable->siBlockIndexLP(internalSite-1,iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=indexTable->aiBlockIndexLP(internalSite-1,iBlock,j);
	  for(int aim=0;aim<lDL;++aim){
	    innercontainer.global_access(siB,aiB,bim,aim)+=source[pctrIndex(aim,bim,aimB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }

#pragma omp parallel for private(lBlockSize,rBlockSize,siB,aimB,siS,biS,bimS,sipS)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(internalSite-1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	siB=indexTable->siBlockIndexLP(internalSite-1,iBlock,k);
	aimB=indexTable->aimBlockIndexLP(internalSite-1,iBlock,k);
	for(int bi=0;bi<lDwR;++bi){
	  outerContainer.global_access(siB,bi,aip,aimB)=0;
	}
	for(int nSparse=0;nSparse<sparseSize;++nSparse){
	  siS=siIndices[nSparse];
	  if(siS==siB){
	    biS=biIndices[nSparse];
	    bimS=bimIndices[nSparse];
	    sipS=sipIndices[nSparse];
	    outerContainer.global_access(siB,biS,aip,aimB)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,aip,bimS,aimB);
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcRightContractionQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  int const numBlocks=indexTable->numBlocksRP(internalSite+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  getLocalDimensions(i+1);
  lDwL=H.locDimL();
  lDwR=H.locDimR();
  tmpContainer<T> outercontainer(lDL,lDwL,ld,lDR);
  //The calculation of the first two contractions has to be done in other functions, too. It therefore has an extra function. 
  calcRightOuterContainerQNOpt(i,internalSite,siteMatrixState,H,source,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aimB,aiB,siB)  
  for(int aimp=0;aimp<lDL;++aimp){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeRP(internalSite+1,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(internalSite+1,iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=indexTable->aimBlockIndexRP(internalSite+1,iBlock,k);
	for(int bim=0;bim<lDwL;++bim){
	  simpleContainer=0;
	  for(int j=0;j<rBlockSize;++j){
	    siB=indexTable->siBlockIndexRP(internalSite+1,iBlock,j);
	    aiB=indexTable->aiBlockIndexRP(internalSite+1,iBlock,j);
	    simpleContainer+=conj(siteMatrixState[stateIndex(siB,aiB,aimB)])*outercontainer.global_access(aimp,bim,siB,aiB);
	  }
	  target[pctrIndex(aimB,bim,aimp)]=simpleContainer;
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcRightOuterContainerQNOpt(int i, int internalSite, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  int const numBlocks=indexTable->numBlocksRP(internalSite+1);
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  getLocalDimensions(i+1);
  T const *siteMatrixH;
  int const *biIndices, *bimIndices, *siIndices, *sipIndices;
  int const sparseSize=H.numEls();
  int biS, bimS, siS, sipS;
  lDwL=H.locDimL();
  lDwR=H.locDimR();
  H.sparseSubMatrixStart(siteMatrixH);
  H.biSubIndexArrayStart(biIndices);
  H.bimSubIndexArrayStart(bimIndices);
  H.siSubIndexArrayStart(siIndices);
  H.sipSubIndexArrayStart(sipIndices);
  tmpContainer<T> innercontainer(ld,lDwR,lDR,lDL);
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
      lBlockSize=indexTable->lBlockSizeRP(internalSite+1,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(internalSite+1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=indexTable->aiBlockIndexRP(internalSite+1,iBlock,j);
	siB=indexTable->siBlockIndexRP(internalSite+1,iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimB=indexTable->aimBlockIndexRP(internalSite+1,iBlock,k);
	  for(int bi=0;bi<lDwR;++bi){
	    innercontainer.global_access(siB,bi,ai,aimB)+=source[pctrIndex(ai,bi,aiB)]*siteMatrixState[stateIndex(siB,aiB,aimB)];
	  }
	}
      }
    }
  }
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,siB,siS,sipS,biS,bimS)
  for(int aim=0;aim<lDL;++aim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      rBlockSize=indexTable->rBlockSizeRP(internalSite+1,iBlock);
      for(int j=0;j<rBlockSize;++j){
	siB=indexTable->siBlockIndexRP(internalSite+1,iBlock,j);
	aiB=indexTable->aiBlockIndexRP(internalSite+1,iBlock,j);
	for(int bim=0;bim<lDwL;++bim){
	  outerContainer.global_access(aim,bim,siB,aiB)=0;
	}
	for(int nSparse=0;nSparse<sparseSize;++nSparse){
	  siS=siIndices[nSparse];
	  if(siS==siB){
	    biS=biIndices[nSparse];
	    bimS=bimIndices[nSparse];
	    sipS=sipIndices[nSparse];
	    outerContainer.global_access(aim,bimS,siB,aiB)+=siteMatrixH[nSparse]*innercontainer.global_access(sipS,biS,aiB,aim);
	  }
	}
      }
    }
  }
}


#endif
