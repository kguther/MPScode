#ifndef CALCULATION_OF_TENSOR_CONTRACTIONS
#define CALCULATION_OF_TENSOR_CONTRACTIONS

#include "dimensionTable.h"
#include "templates/tmpContainer.h"
#include "templates/mpoSiteTensor.h"
#include "siteQNOrderMatrix.h"

//---------------------------------------------------------------------------------------------------//
// Core functionality of locally contracting some MPS given by siteMatrixState with some MPO given by siteMatrixH. Does not care how the MPS is build -> featured in both baseMeasurement and uncachedMeasurement.
//---------------------------------------------------------------------------------------------------//

class contractor{
 public:
  contractor(){}
 contractor(int DwIn, dimensionTable const &dimInfoIn):dimInfo(dimInfoIn), Dw(DwIn){D=dimInfo.D();}

  template<typename T>
    void calcLeftContraction(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcRightContraction(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  
  template<typename T>
    void calcLeftOuterContainer(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);
 template<typename T>
    void calcRightOuterContainer(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);

  
  //Versions of contractions optimized for systems with conserved QNs (are called automatically when available)
  template<typename T>
    void calcLeftContraction(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcRightContraction(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target);
  template<typename T>
    void calcLeftOuterContainer(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);
  template<typename T>
    void calcRightOuterContainer(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer);

 private:
  dimensionTable dimInfo;
  int lDL, lDR, ld, lDwR, lDwL;
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
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcLeftContraction(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  getLocalDimensions(i-1);
  tmpContainer<T> outerContainer(ld,lDwR,lDR,lDL);
  calcLeftOuterContainer(i,siteMatrixState,H,source,outerContainer);
  for(int bi=0;bi<lDwR;++bi){
    for(int aip=0;aip<lDR;++aip){
      for(int ai=0;ai<lDR;++ai){
	simpleContainer=0.0;
	for(int si=0;si<ld;++si){
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=outerContainer.global_access(si,bi,ai,aimp)*conj(siteMatrixState[stateIndex(si,aip,aimp)]);
	  }
	}
	target[pctrIndex(aip,bi,ai)]=simpleContainer;
      }
    }
  }
}

template<typename T>
void contractor::calcRightContraction(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  getLocalDimensions(i+1);
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  tmpContainer<T> outerContainer(lDL,lDwL,ld,lDR);
  calcRightOuterContainer(i,siteMatrixState,H,source,outerContainer);
  for(int aim=0;aim<lDL;++aim){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=conj(siteMatrixState[stateIndex(si,ai,aim)])*outerContainer.global_access(aimp,bim,si,ai);
	  }
	}
	target[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
}


template<typename T>
void contractor::calcLeftOuterContainer(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  T simpleContainer;
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  getLocalDimensions(i-1);
  tmpContainer<T> innerContainer(ld,lDwL,lDL,lDR);
  for(int si=0;si<ld;++si){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	for(int ai=0;ai<lDR;++ai){
	  simpleContainer=0.0;
	  for(int aim=0;aim<lDL;++aim){
	    simpleContainer+=siteMatrixState[stateIndex(si,ai,aim)]*source[pctrIndex(aimp,bim,aim)];
	  }
	  innerContainer.global_access(si,bim,aimp,ai)=simpleContainer;
	}
      }
    }
  }
  T const *siteMatrixH;
  H.subMatrixStart(siteMatrixH);
  for(int sip=0;sip<ld;++sip){
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0.0;
	  for(int si=0;si<ld;++si){
	    for(int bim=0;bim<lDwL;++bim){
	      simpleContainer+=innerContainer.global_access(si,bim,aimp,ai)*siteMatrixH[operatorIndex(sip,si,bi,bim)];
	    }
	  }
	  outerContainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
}

template<typename T>
void contractor::calcRightOuterContainer(int i, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  T simpleContainer;
  getLocalDimensions(i+1);
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  tmpContainer<T > innerContainer(ld,lDwR,lDR,lDL);
  for(int sip=0;sip<ld;++sip){                                                       
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0.0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=source[pctrIndex(ai,bi,aip)]*siteMatrixState[stateIndex(sip,aip,aimp)];
	  }
	  innerContainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  T const *siteMatrixH;
  H.subMatrixStart(siteMatrixH);
  for(int aimp=0;aimp<lDL;++aimp){
    for(int bim=0;bim<lDwL;++bim){
      for(int ai=0;ai<lDR;++ai){
	for(int si=0;si<ld;++si){
	  simpleContainer=0.0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=innerContainer.global_access(sip,bi,ai,aimp)*siteMatrixH[operatorIndex(si,sip,bi,bim)];
	    }
	  }
	  outerContainer.global_access(aimp,bim,si,ai)=simpleContainer;
	}
      }
    }
  }
}


//---------------------------------------------------------------------------------------------------//

template<typename T>
void contractor::calcLeftContraction(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  int const numBlocks=localIndexTable.numBlocksLP();
  int lBlockSize, rBlockSize;
  int aiB, siB, aimB;
  getLocalDimensions(i-1);
  lDwR=H.locDimR();
  lDwL=H.locDimL();
  //container arrays to significantly reduce computational effort by storing intermediate results
  tmpContainer<T> outercontainer(ld,lDwR,lDR,lDL);
  calcLeftOuterContainer(i,localIndexTable,siteMatrixState,H,source,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aiB,siB,aimB)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0.0;
	  for(int k=0;k<lBlockSize;++k){
	    siB=localIndexTable.siBlockIndexLP(iBlock,k);
	    aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
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
void contractor::calcLeftOuterContainer(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  int const *biIndices, *bimIndices, *siIndices, *sipIndices;
  T const *siteMatrixH;
  int biS, bimS, siS, sipS, aiB, aimB, siB;
  int lBlockSize, rBlockSize;
  int const sparseSize=H.numEls();
  int const numBlocks=localIndexTable.numBlocksLP();
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
	  innercontainer.global_access(sip,aip,bim,aim)=0.0;
	}
      }
    }
  }
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,aimB,siB)
  for(int bim=0;bim<lDwL;++bim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
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
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	for(int bi=0;bi<lDwR;++bi){
	  outerContainer.global_access(siB,bi,aip,aimB)=0.0;
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
void contractor::calcRightContraction(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, T *const target){
  T simpleContainer;
  int const numBlocks=localIndexTable.numBlocksRP();
  int aiB, siB, aimB;
  int lBlockSize, rBlockSize;
  getLocalDimensions(i+1);
  lDwL=H.locDimL();
  lDwR=H.locDimR();
  tmpContainer<T> outercontainer(lDL,lDwL,ld,lDR);
  //The calculation of the first two contractions has to be done in other functions, too. It therefore has an extra function. 
  calcRightOuterContainer(i,localIndexTable,siteMatrixState,H,source,outercontainer);
#pragma omp parallel for private(simpleContainer,lBlockSize,rBlockSize,aimB,aiB,siB)  
  for(int aimp=0;aimp<lDL;++aimp){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeRP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeRP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexRP(iBlock,k);
	for(int bim=0;bim<lDwL;++bim){
	  simpleContainer=0.0;
	  for(int j=0;j<rBlockSize;++j){
	    siB=localIndexTable.siBlockIndexRP(iBlock,j);
	    aiB=localIndexTable.aiBlockIndexRP(iBlock,j);
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
void contractor::calcRightOuterContainer(int i, siteQNOrderMatrix const &localIndexTable, T *const siteMatrixState, mpoSiteTensor<T> const &H, T *const source, tmpContainer<T> &outerContainer){
  int const numBlocks=localIndexTable.numBlocksRP();
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
	  innercontainer.global_access(sip,bi,ai,aimp)=0.0;
	}
      }
    }
  }
  //Contraction works similar to the general case. Here, only entries fulfilling the QN constraint are used
#pragma omp parallel for private(lBlockSize,rBlockSize,aiB,siB,aimB)
  for(int ai=0;ai<lDR;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeRP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeRP(iBlock);
      for(int j=0;j<rBlockSize;++j){
	aiB=localIndexTable.aiBlockIndexRP(iBlock,j);
	siB=localIndexTable.siBlockIndexRP(iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimB=localIndexTable.aimBlockIndexRP(iBlock,k);
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
      rBlockSize=localIndexTable.rBlockSizeRP(iBlock);
      for(int j=0;j<rBlockSize;++j){
	siB=localIndexTable.siBlockIndexRP(iBlock,j);
	aiB=localIndexTable.aiBlockIndexRP(iBlock,j);
	for(int bim=0;bim<lDwL;++bim){
	  outerContainer.global_access(aim,bim,siB,aiB)=0.0;
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
