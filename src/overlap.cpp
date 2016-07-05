#include "overlap.h"
#include "tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// Tested computation of partial contractions, they work as intended. Checked computation of F (and its
// usage in network).
//---------------------------------------------------------------------------------------------------//

overlap::overlap(){
  Lctr=0;
  Rctr=0;
  psi=0;
  phi=0;
}

//---------------------------------------------------------------------------------------------------//

//Comes useful for testing normalization
overlap::overlap(mps const*const psi, mps const*const phi):
  Lctr(0),
  Rctr(0)
{
  loadMPS(psi,phi);
}

//---------------------------------------------------------------------------------------------------//

overlap::overlap(overlap const &source){
  ovCpy(source);
}

//---------------------------------------------------------------------------------------------------//

overlap::~overlap(){
  delete[] Lctr;
  delete[] Rctr;
}

//---------------------------------------------------------------------------------------------------//

overlap& overlap::operator=(overlap const &source){
  ovCpy(source);
  return *this;
}

//---------------------------------------------------------------------------------------------------//

void overlap::loadMPS(mps const*const psiIn, mps const*const phiIn){
  //Get two states of which the overlap shall be computed
  phi=phiIn;
  psi=psiIn;
  D=psiIn->maxDim();
  L=psiIn->length();
  d=psiIn->siteDim();
  delete[] Lctr;
  delete[] Rctr;
  //Adjust buffers for caching
  Lctr=new lapack_complex_double[L*D*D];
  Rctr=new lapack_complex_double[L*D*D];
  //Initialize F-Matrix
  F.initialize(psiIn->getDimInfo());
  getF();
}

//---------------------------------------------------------------------------------------------------//

void overlap::ovCpy(overlap const &source){
  if(source.psi && source.phi){
    loadMPS(source.psi,source.phi);
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::subContractionStartLeft(lapack_complex_double *&pStart, int i){
  pStart=Lctr+i*D*D;
}

//---------------------------------------------------------------------------------------------------//

void overlap::subContractionStartRight(lapack_complex_double *&pStart, int i){
  pStart=Rctr+i*D*D;
}

//---------------------------------------------------------------------------------------------------//
// This is basically a more efficient way of measuring identity. It works just the same way as in 
// the measurement classes, just without the intermediate step of inserting some MPO.
//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterLeft(int i){
  int lDR, lDL, ld;
  lDL=phi->locDimL(i);
  lDR=phi->locDimR(i);
  ld=phi->locd(i);
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  lapack_complex_double simpleContainer;
  lapack_complex_double *source, *target;
  if(i>0){
    subContractionStartLeft(source,i-1);
  }
  else{
    lapack_complex_double zone=1.0;
    source=&zone;
  }
  subContractionStartLeft(target,i);
  if(phi->indexTable().nQNs() && psi->indexTable().nQNs()){
    calcCtrIterLeftQNOpt(i,source,target);
  }
  else{
    for(int si=0;si<ld;++si){
      for(int aip=0;aip<lDR;++aip){
	for(int aim=0;aim<lDL;++aim){
	  simpleContainer=0;
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=source[pCtrLocalIndex(aimp,aim)]*phi->global_access(i,si,aip,aimp);
	  }
	  innerContainer.global_access(0,si,aip,aim)=simpleContainer;
	}
      }
    }
    for(int aip=0;aip<lDR;++aip){
      for(int ai=0;ai<lDR;++ai){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int aim=0;aim<lDL;++aim){
	    simpleContainer+=innerContainer.global_access(0,si,aip,aim)*conj(psi->global_access(i,si,ai,aim));
	  }
	}
	target[pCtrLocalIndex(aip,ai)]=simpleContainer;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterRight(int i){
  int lDR, lDL, ld;
  lDL=phi->locDimL(i);
  lDR=phi->locDimR(i);
  ld=phi->locd(i);
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  lapack_complex_double simpleContainer;
  lapack_complex_double *source, *target;
  if(i<(L-1)){
    subContractionStartRight(source,i+1);
  }
  else{
    lapack_complex_double zone=1.0;
    source=&zone;
  }
  subContractionStartRight(target,i);
  //Using QN-fixed MPS and QN-nonfixed MPS in the same run is not supported
  if(phi->indexTable().nQNs() && psi->indexTable().nQNs()){
    //Version optimized for conserved QNs. 
    calcCtrIterRightQNOpt(i,source,target);
  }
  else{
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=source[pCtrLocalIndex(aip,ai)]*phi->global_access(i,si,aip,aimp);
	  }
	  innerContainer.global_access(0,si,ai,aimp)=simpleContainer;
	}
      }
    }
    for(int aimp=0;aimp<lDL;++aimp){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=innerContainer.global_access(0,si,ai,aimp)*conj(psi->global_access(i,si,ai,aim));
	  }
	}
	target[pCtrLocalIndex(aimp,aim)]=simpleContainer;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// Optimized versions of calcCtrIterLeft/Right for the case both, psi and phi are QN-Blocked.
//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterLeftQNOpt(int i, lapack_complex_double const*const source, lapack_complex_double *const target){
  int lDL=phi->locDimL(i);
  int ld=phi->locd(i);
  int lDR=phi->locDimR(i);
  siteQNOrderMatrix localIndexTable=phi->indexTable().getLocalIndexTable(i);
  lapack_complex_double const *localMatrix;
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	innerContainer.global_access(0,si,ai,aim)=0;
      }
    }
  }
  phi->subMatrixStart(localMatrix,i);
  int numBlocks, lBlockSize, rBlockSize;
  int aiB, aimB, siB;
  numBlocks=localIndexTable.numBlocksLP();
  for(int aim=0;aim<lDL;++aim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
	  innerContainer.global_access(0,siB,aiB,aim)+=source[pCtrLocalIndex(aimB,aim)]*localMatrix[aimB+aiB*lDL+siB*lDL*lDR];
	}
      }
    }
  }
  for(int ai=0;ai<lDR;++ai){
    for(int aip=0;aip<lDR;++aip){
      target[pCtrLocalIndex(ai,aip)]=0;
    }
  }
  localIndexTable=psi->indexTable().getLocalIndexTable(i);
  psi->subMatrixStart(localMatrix,i);
  numBlocks=localIndexTable.numBlocksLP();
  for(int ai=0;ai<lDR;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
	  target[pCtrLocalIndex(ai,aiB)]+=innerContainer.global_access(0,siB,ai,aimB)*conj(localMatrix[aimB+aiB*lDL+siB*lDL*lDR]);
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterRightQNOpt(int i, lapack_complex_double const*const source, lapack_complex_double *const target){
  //works also if phi and psi have different L/R-basis at site i
  int lDL=phi->locDimL(i);
  int ld=phi->locd(i);
  int lDR=phi->locDimR(i);
  siteQNOrderMatrix localIndexTable=phi->indexTable().getLocalIndexTable(i);
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  lapack_complex_double const *localMatrix;
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	innerContainer.global_access(0,si,ai,aim)=0;
      }
    }
  }
  int numBlocks, lBlockSize, rBlockSize;
  int aiB, aimB, siB;
  numBlocks=localIndexTable.numBlocksLP();
  phi->subMatrixStart(localMatrix,i);
  for(int ai=0;ai<lDR;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
	  innerContainer.global_access(0,siB,ai,aimB)+=source[pCtrLocalIndex(aiB,ai)]*localMatrix[aimB+aiB*lDL+siB*lDL*lDR];
	}
      }
    }
  }
  for(int aim=0;aim<lDL;++aim){
    for(int aimp=0;aimp<lDL;++aimp){
      target[pCtrLocalIndex(aim,aimp)]=0;
    }
  }
  localIndexTable=psi->indexTable().getLocalIndexTable(i);
  numBlocks=localIndexTable.numBlocksLP();
  psi->subMatrixStart(localMatrix,i);
  for(int aim=0;aim<lDL;++aim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimB=localIndexTable.aimBlockIndexLP(iBlock,k);
	siB=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiB=localIndexTable.aiBlockIndexLP(iBlock,j);
	  target[pCtrLocalIndex(aim,aimB)]+=innerContainer.global_access(0,siB,aiB,aim)*conj(localMatrix[aimB+aiB*lDL+siB*lDL*lDR]);
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// The F matrix is given by F=d/dconj(M)_i <psi|phi>, that is, the complete contraction of the two 
// states except for site i of psi. 
// Structurally, F is just an mps but without normalization methods.
// updateF(..) is used to update the F matrix of some site after 
// one of the partial contractions Lctr or Rctr has changed.
// getF() computes the F matrix for all sites.
//---------------------------------------------------------------------------------------------------//

void overlap::updateF(int i){
  int lDR, lDL, ld;
  lDL=phi->locDimL(i);
  lDR=phi->locDimR(i);
  ld=phi->locd(i);
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  lapack_complex_double simpleContainer;
  lapack_complex_double zone=1.0;
  lapack_complex_double *leftPart, *rightPart;
  lapack_complex_double const *localMatrix;
  if(i>0){
    subContractionStartLeft(leftPart,i-1);
  }
  else{
    leftPart=&zone;
  }
  if(i<L-1){
    subContractionStartRight(rightPart,i+1);
  }
  else{
    rightPart=&zone;
  }
  phi->subMatrixStart(localMatrix,i);
  for(int si=0;si<ld;++si){
    for(int aip=0;aip<lDR;++aip){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer+=leftPart[pCtrLocalIndex(aimp,aim)]*localMatrix[aimp+aip*lDL+si*lDL*lDR];
	}
	innerContainer.global_access(0,si,aip,aim)=simpleContainer;
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int aip=0;aip<lDR;++aip){
	  simpleContainer+=innerContainer.global_access(0,si,aip,aim)*rightPart[pCtrLocalIndex(aip,ai)];
	}
	F.global_access(i,si,ai,aim)=simpleContainer;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::getF(){
  for(int i=0;i<L-1;++i){
    calcCtrIterLeft(i);
  }
  updateF(L-1);
  for(int i=L-1;i>0;--i){
    calcCtrIterRight(i);
    updateF(i-1);
  }
}

//---------------------------------------------------------------------------------------------------//
// These functions update the complete overlap after the on-site matrices of the state psi (first argument)
// have changed on site i.
//---------------------------------------------------------------------------------------------------//

void overlap::stepLeft(int i){
  //i is the source site of the step, i.e. site i-1 is updated
  calcCtrIterRight(i);
  if(i>0){
    updateF(i-1);
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::stepRight(int i){
  calcCtrIterLeft(i);
  updateF(i+1);
}

//---------------------------------------------------------------------------------------------------//
// These two functions return the scalar product of psi and phi. fullOverlap() only returns the value, 
// but it has to be computed manually previously whereas getFullOverlap() re-evaluates the complete
// contraction.
//---------------------------------------------------------------------------------------------------//

lapack_complex_double overlap::fullOverlap(){
  return Rctr[0];
}

//---------------------------------------------------------------------------------------------------//

lapack_complex_double overlap::getFullOverlap(){
  for(int i=L-1;i>=0;--i){
    calcCtrIterRight(i);
  }
  return Rctr[0];
}

//---------------------------------------------------------------------------------------------------//
// Not in use currently. 
//---------------------------------------------------------------------------------------------------//

lapack_complex_double overlap::applyF(lapack_complex_double *vec, int i){
  lapack_complex_double simpleContainer=0.0;
  int lDR, lDL, ld;
  lDL=(*phi).locDimL(i);
  lDR=(*phi).locDimR(i);
  ld=(*phi).locd(i);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer+=F.global_access(i,si,ai,aim)*conj(vec[aim+lDL*ai+lDR*lDL*si]);
      }
    }
  }
  return simpleContainer;
}


