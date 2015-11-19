#include <lapacke.h>
#include <iostream>
#include "overlap.h"
#include "tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// Tested computation of partial contractions, they work as intended. Check computation of F (and its
// usage in network).
//---------------------------------------------------------------------------------------------------//

overlap::overlap(){
  Lctr=0;
  Rctr=0;
}

//---------------------------------------------------------------------------------------------------//

overlap::~overlap(){
  delete[] Lctr;
  delete[] Rctr;
}

//---------------------------------------------------------------------------------------------------//

void overlap::loadMPS(mps *psiIn, mps *phiIn){
  phi=phiIn;
  psi=psiIn;
  D=(*psiIn).maxDim();
  L=(*psiIn).length();
  d=(*psiIn).siteDim();
  Lctr=new lapack_complex_double[L*D*D];
  Rctr=new lapack_complex_double[L*D*D];
  F.generate(d,D,L);
  getF();
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

void overlap::calcCtrIterLeft(int const i){
  int lDR, lDL, ld;
  lDL=(*phi).locDimL(i);
  lDR=(*phi).locDimR(i);
  ld=(*phi).locd(i);
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
	  simpleContainer+=source[pCtrLocalIndex(aimp,aim)]*(*phi).global_access(i,si,aip,aimp);
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
	  simpleContainer+=innerContainer.global_access(0,si,aip,aim)*(*psi).global_access(i,si,ai,aim);
	}
      }
      target[pCtrLocalIndex(aip,ai)]=simpleContainer;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void overlap::calcCtrIterRight(int const i){
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
	  simpleContainer+=source[pCtrLocalIndex(aip,ai)]*(*phi).global_access(i,si,aip,aimp);
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
	  simpleContainer+=innerContainer.global_access(0,si,ai,aimp)*(*psi).global_access(i,si,ai,aim);
	}
      }
      target[pCtrLocalIndex(aimp,aim)]=simpleContainer;
    }
  }
}

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
	  simpleContainer+=leftPart[pCtrLocalIndex(aimp,aim)]*(*phi).global_access(i,si,aip,aimp);
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


