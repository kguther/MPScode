#include <arcomp.h>
#include <iostream>
#include "parameters.h"
#include "optHMatrix.h"
#include "tmpContainer.h"

optHMatrix::optHMatrix(arcomplex<double> *Rin, arcomplex<double> *Lin, arcomplex<double> *Hin, parameters pars, int i):
  Rctr(Rin),
  Lctr(Lin),
  H(Hin),
  d(pars.d),
  D(pars.D),
  L(pars.L),
  Dw(pars.Dw)
{
  icrit=L/2;
  for(int j=0;j<L/2;j++){
    if(pow(d,j+1)>D){
      icrit=j;
      break;
    }
  }
  //Lengthy initialization of local Matrix dimension
  lDwL=Dw;
  lDwR=Dw;
  if(i==0){
    lDwL=1;
  }
  if(i==(L-1)){
    lDwR=1;
  }
  if(i<=icrit){
    lDL=pow(d,i);
  }
  else{
    if(i<=L-icrit-1){
    lDL=D;
  }
    else{
      lDL=pow(d,L-i);
    }
  }
  if(i<icrit){
    lDR=pow(d,i+1);
  }
  else{
    if(i<=L-icrit-2){
    lDR=D;
  }
    else{
      lDR=pow(d,L-i-1);
    }
  }
  //Dimension of H is obviously a necessary information for ARPACK++
  dimension=d*lDL*lDR;
  std::cout<<"Current eigenvalue problem Dimension: "<<dimension<<std::endl;
}

//---------------------------------------------------------------------------------------------------//
// This is the multiplication of some vector v with an optHMatrix, the output is returned in w
// The functiion arguments are required this way by ARPACK++ 
//---------------------------------------------------------------------------------------------------//

void optHMatrix::MultMv(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innercontainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outercontainer(d,lDwL,lDR,lDL);
  int nNzero;
  arcomplex<double> simpleContainer;
  //Similar to the calculation of partial contractions, we use optimal bracketing to reuse any intermediate results. This greatly reduces the computational effort and is much faster than storing H in a sparse format and using the internal ARPACK++ matrix classes
  for(int sip=0;sip<d;++sip){
    for(int aimp=0;aimp<lDL;++aimp){
      for(int ai=0;ai<lDR;++ai){
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=Rctr[ctrIndex(ai,bi,aip)]*v[vecIndex(sip,aip,aimp)];
	  }
	  innercontainer.global_access(sip,aimp,ai,bi)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<d;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int bim=0;bim<lDwL;++bim){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0;
	  for(int sip=0;sip<d;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=innercontainer.global_access(sip,aimp,ai,bi)*H[hIndex(si,sip,bi,bim)];
	    }
	  }
	  outercontainer.global_access(si,bim,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<d;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	simpleContainer=0;
	for(int bim=0;bim<lDwL;++bim){
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=outercontainer.global_access(si,bim,ai,aimp)*Lctr[ctrIndex(aim,bim,aimp)];
	  }
	}
	w[vecIndex(si,ai,aim)]=simpleContainer;
      }
    }
  }
}
