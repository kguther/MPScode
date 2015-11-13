#include <lapacke.h>
#include "overlap.h"

overlap::overlap(){
  Lctr=0;
  Rctr=0;
}

overlap::~overlap(){
  delete[] Lctr;
  delete[] Rctr;
}

void overlap::loadMPS(mps *psiIn, mps *phiIn){
  phi=phiIn;
  psi=psiIn;
  D=(*psiIn).maxDim();
  L=(*psiIn).length();
  d=(*psiIn).siteDim();
  Lctr=new lapack_complex_double[L*D*D];
  Rctr=new lapack_complex_double[L*D*D];
  Lctr[0]=1;
  Rctr[(L-1)*D*D]=1;
}

void overlap::subContractionStartLeft(lapack_complex_double *&pStart, int const i){
  pStart=Lctr+i*D*D;
}

void overlap::subContractionStartRight(lapack_complex_double *&pStart, int const i){
  pStart=Rctr+i*D*D;
}

void overlap::calcCtrIterLeft(int const i){
  if(i>0){
    int lDR, lDL, ld;
    lDL=(*phi).locDimL(i);
    lDR=(*phi).locDimR(i);
    ld=(*phi).locd(i);
    tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
    lapack_complex_double simpleContainer;
    lapack_complex_double *source, *target;
    subContractionStartLeft(source,i-1);
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
	    simpleContainer+=innerContainer.global_access[0,si,aip,aim]*(*psi).global_access(i,si,ai,aim);
	  }
	}
	target[pCtrLocalIndex(aip,ai)]=simpleContainer;
      }
    }
  }
  else{
    Lctr[0]=1;
  }
}

void overlap::calcCtrIterRight(int const i){
  if(i<(L-1)){
    int lDR, lDL, ld;
    lDL=(*phi).locDimL(i);
    lDR=(*phi).locDimR(i);
    ld=(*phi).locd(i);
    tmpContainer<lapack_complex_double> innerContainer(1,ld,lDR,lDL);
    lapack_complex_double simpleContainer;
    lapack_complex_double *source, *target;
    subContractionStartLeft(source,i+1);
    subContractionStartLeft(target,i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=source[pCtrLocalIndex(aip,ai)]*(*phi).global_access(i,si,aip,aimp);
	  }
	  innerContainer.global_access(0,si,ai,aim)=simpleContainer;
	}
      }
    }
    for(int aimp=0;aimp<lDL;++aimp){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=innerContainer.global_access[0,si,ai,aim]*(*psi).global_access(i,si,ai,aim);
	  }
	}
	target[pCtrLocalIndex(aimp,aim)]=simpleContainer;
      }
    }
  }
  else{
    Rctr[(L-1)*D*D]=1;
  }
}


