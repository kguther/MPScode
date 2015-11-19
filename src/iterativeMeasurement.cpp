#include <lapacke.h>
#include <iostream>
#include "mpo.h"
#include "mps.h"
#include "iterativeMeasurement.h"
#include "tmpContainer.h"
#include "pContraction.h"

iterativeMeasurement::iterativeMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::initialize(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  MPOperator=MPOperatorIn;
  MPState=MPStateIn;
  initializeBase(MPOperatorIn,MPStateIn);
  Lctr.initialize((*MPOperator).length(),(*MPState).maxDim(),(*MPOperator).maxDim());
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state.
//---------------------------------------------------------------------------------------------------//

void iterativeMeasurement::calcCtrIterLeft(int const i){
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr, *targetPctr;
  Lctr.subContractionStart(sourcePctr,i-1);
  Lctr.subContractionStart(targetPctr,i);
  getLocalDimensions(i-1);
  //container arrays to significantly reduce computational effort by storing intermediate results
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwL,lDL,lDR);
  tmpContainer<lapack_complex_double> outercontainer(ld,lDwR,lDR,lDL);
  //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
  for(int sip=0;sip<ld;++sip){
    for(int bim=0;bim<lDwL;++bim){
      for(int aim=0;aim<lDL;++aim){
	for(int aip=0;aip<lDR;++aip){
	  simpleContainer=0;
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=sourcePctr[pctrIndex(aim,bim,aimp)]*(*MPState).global_access(i-1,sip,aip,aimp);
	  }
	  innercontainer.global_access(sip,bim,aim,aip)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int bi=0;bi<lDwR;++bi){
      for(int aip=0;aip<lDR;++aip){
	for(int aim=0;aim<lDL;++aim){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bim=0;bim<lDwL;++bim){
	      simpleContainer+=(*MPOperator).global_access(i-1,si,sip,bi,bim)*innercontainer.global_access(sip,bim,aim,aip);
	    }
	  }
	  outercontainer.global_access(si,bi,aip,aim)=simpleContainer;
	}
      }
    }
  }
  for(int ai=0;ai<lDR;++ai){
    for(int bi=0;bi<lDwR;++bi){
      for(int aip=0;aip<lDR;++aip){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int aim=0;aim<lDL;++aim){
	    simpleContainer+=conj((*MPState).global_access(i-1,si,ai,aim))*outercontainer.global_access(si,bi,aip,aim);
	  }
	}
	targetPctr[pctrIndex(ai,bi,aip)]=simpleContainer;
      }
    }
  }
  std::cout<<"Completed calculation of partial contraction at site "<<i<<std::endl;
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
  int L=(*MPOperator).length();
  if(direction==1){
    Rctr.global_access(L-1,0,0,0)=lapack_make_complex_double(1.0,0.0);
    for(int i=L-2;i>=0;--i){
      calcCtrIterRight(i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Lctr.global_access(0,0,0,0)=lapack_make_complex_double(1.0,0.0);
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
