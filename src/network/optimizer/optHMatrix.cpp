#include <arcomp.h>
#include <iostream>
#include <time.h>
#include "parameters.h"
#include "optHMatrix.h"
#include "tmpContainer.h"

optHMatrix::optHMatrix(arcomplex<double> *Rin, arcomplex<double> *Lin, mpo<arcomplex<double> > const *Hin, dimensionTable const &dimInfo, int Dwin, int iIn, projector *excitedStateP, double shiftin, std::vector<quantumNumber> *conservedQNsin):
  Rctr(Rin),
  Lctr(Lin),
  Dw(Dwin),
  i(iIn),
  shift(shiftin),
  P(excitedStateP),
  HMPO(Hin),
  conservedQNs(conservedQNsin)
{
  D=dimInfo.D();
  lDwR=Hin->locDimR(i);
  lDwL=Hin->locDimL(i);
  lDR=dimInfo.locDimR(i);
  lDL=dimInfo.locDimL(i);
  d=dimInfo.locd(i);
  //Dimension of H is obviously a necessary information for ARPACK++
  dimension=d*lDL*lDR;
  //std::cout<<"Current eigenvalue problem dimension: "<<dimension<<std::endl;
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
  if(P){
    P->project(v,i);
  }
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
  int const *biIndices, *siIndices, *bimIndices, *sipIndices;
  int const sparseSize=HMPO->numEls(i);
  arcomplex<double> const *H;
  HMPO->biSubIndexArrayStart(biIndices,i);
  HMPO->siSubIndexArrayStart(siIndices,i);
  HMPO->bimSubIndexArrayStart(bimIndices,i);
  HMPO->sipSubIndexArrayStart(sipIndices,i);
  HMPO->sparseSubMatrixStart(H,i);

  #pragma omp parallel for schedule(static,1)
  for(int si=0;si<d;++si){
    for(int bim=0;bim<lDwL;++bim){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  outercontainer.global_access(si,bim,ai,aim)=0.0;
	}
      }
    }
  }

  for(int ai=0;ai<lDR;++ai){
    for(int aimp=0;aimp<lDL;++aimp){
      for(int nSparse=0;nSparse<sparseSize;++nSparse){
	outercontainer.global_access(siIndices[nSparse],bimIndices[nSparse],ai,aimp)+=innercontainer.global_access(sipIndices[nSparse],aimp,ai,biIndices[nSparse])*H[nSparse];
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
	w[vecIndex(si,ai,aim)]=simpleContainer+shift*v[vecIndex(si,ai,aim)];
      }
    }
  }
  if(P){
    P->project(w,i);
  }
}

//---------------------------------------------------------------------------------------------------//

void optHMatrix::MultMvQNConserving(arcomplex<double> *v, arcomplex<double> *w){
  //FOR MORE SOPHISTICATED QN CONSERVING MULTIPLICATION USE THE blockHMatrix CLASS
  clock_t curtime;
  curtime=clock();
  projectQN(v);
  MultMv(v,w);
  projectQN(w);
  if(0){
  curtime=clock()-curtime;
  std::cout<<"Matrix multiplication took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  exit(1);
  }
}

//---------------------------------------------------------------------------------------------------//

void optHMatrix::projectQN(arcomplex<double> *v){
  //This just sets all entries forbidden by QN constraint to zero. Important: the QN constraint is sufficient, but not necessary to have a state with some fixed QN.
  for(int si=0;si<d;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	for(int iQN=0;iQN<conservedQNs->size();++iQN){
	  if((*conservedQNs)[iQN].qnConstraint(i,si,ai,aim)){
	    v[vecIndex(si,ai,aim)]=0;
	  }
	}
      }
    }
  }
}
