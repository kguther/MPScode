#include <iostream>
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
  phi=phiIn;
  psi=psiIn;
  D=psiIn->maxDim();
  L=psiIn->length();
  d=psiIn->siteDim();
  delete[] Lctr;
  delete[] Rctr;
  Lctr=new lapack_complex_double[L*D*D];
  Rctr=new lapack_complex_double[L*D*D];
  F.generate(psiIn->dimInfo);
  getF();
}

//---------------------------------------------------------------------------------------------------//

void overlap::ovCpy(overlap const &source){
  if(source.psi && source.phi){
    loadMPS(source.psi,source.phi);
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::subContractionStartLeft(lapack_complex_double *&pStart, int const i){
  pStart=Lctr+i*D*D;
}

//---------------------------------------------------------------------------------------------------//

void overlap::subContractionStartRight(lapack_complex_double *&pStart, int const i){
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

//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterRight(int i){
  int lDR, lDL, ld;
  lDL=(*phi).locDimL(i);
  lDR=(*phi).locDimR(i);
  ld=(*phi).locd(i);
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

//---------------------------------------------------------------------------------------------------//
// Optimized versions of calcCtrIterLeft/Right for the case both, psi and phi are QN-Blocked.
//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterLeftQNOpt(int i){
  
}

//---------------------------------------------------------------------------------------------------//
// The F matrix is given by F=d/dconj(M)_i <psi|phi>, that is, the complete contraction of the two 
// states except for site i of psi. 
// Structurally, F is just an mps but without normalization methods.
// updateF(..) is used to update the F matrix of some site after 
// one of the partial contractions Lctr or Rctr has changed.
// getF() computes the F matrix for all sites.
//---------------------------------------------------------------------------------------------------//

void overlap::updateF(int const i){
  int lDR, lDL, ld;
  lDL=(*phi).locDimL(i);
  lDR=(*phi).locDimR(i);
  ld=(*phi).locd(i);
  tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
  lapack_complex_double simpleContainer;
  lapack_complex_double zone=1.0;
  lapack_complex_double *leftPart, *rightPart;
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
  for(int si=0;si<ld;++si){
    for(int aip=0;aip<lDR;++aip){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer+=leftPart[pCtrLocalIndex(aimp,aim)]*phi->global_access(i,si,aip,aimp);
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

void overlap::stepLeft(int const i){
  //i is the source site of the step, i.e. site i-1 is updated
  calcCtrIterRight(i);
  if(i>0){
    updateF(i-1);
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::stepRight(int const i){
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

lapack_complex_double overlap::applyF(lapack_complex_double *vec, int const i){
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


