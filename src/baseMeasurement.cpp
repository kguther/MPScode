#include "baseMeasurement.h"
#include <lapacke.h>
#include <iostream>
#include "mpo.h"
#include "mps.h"
#include "tmpContainer.h"
#include "pContraction.h"

baseMeasurement::baseMeasurement(){
}

//---------------------------------------------------------------------------------------------------//

baseMeasurement::baseMeasurement(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn):
  MPOperator(MPOperatorIn),
  MPState(MPStateIn)
{
  initializeBase(MPOperatorIn,MPStateIn);
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::initializeBase(mpo<lapack_complex_double> *MPOperatorIn, mps *MPStateIn){
  Dw=(*MPOperator).maxDim();
  D=(*MPState).maxDim();
  Rctr.initialize((*MPOperator).length(),(*MPState).maxDim(),(*MPOperator).maxDim());
}

//---------------------------------------------------------------------------------------------------//  

void baseMeasurement::getLocalDimensions(int const i){
  lDL=(*MPState).locDimL(i);
  lDR=(*MPState).locDimR(i);
  ld=(*MPState).locd(i);
  lDwR=(*MPOperator).locDimR(i);
  lDwL=(*MPOperator).locDimL(i);
}

//---------------------------------------------------------------------------------------------------//

void baseMeasurement::calcCtrIterRightBase(int const i, lapack_complex_double *targetPctr){
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr;
  Rctr.subContractionStart(sourcePctr,i+1);
  getLocalDimensions(i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,lDwR,lDR,lDL);
  tmpContainer<lapack_complex_double> outercontainer(ld,lDwL,lDL,lDR);
  for(int sip=0;sip<ld;++sip){                                                       
    for(int bi=0;bi<lDwR;++bi){
      for(int ai=0;ai<lDR;++ai){
	for(int aimp=0;aimp<lDL;++aimp){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=sourcePctr[pctrIndex(ai,bi,aip)]*(*MPState).global_access(i+1,sip,aip,aimp);
	  }
	  innercontainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	for(int ai=0;ai<lDR;++ai){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=(*MPOperator).global_access(i+1,si,sip,bi,bim)*innercontainer.global_access(sip,bi,ai,aimp);
	    }
	  }
	  outercontainer.global_access(si,bim,aimp,ai)=simpleContainer;
	}
      }
    }
  }
  for(int aim=0;aim<lDL;++aim){
    for(int bim=0;bim<lDwL;++bim){
      for(int aimp=0;aimp<lDL;++aimp){
	simpleContainer=0;
	for(int si=0;si<ld;++si){
	  for(int ai=0;ai<lDR;++ai){
	    simpleContainer+=conj((*MPState).global_access(i+1,si,ai,aim))*outercontainer.global_access(si,bim,aimp,ai);
	  }
	}
	targetPctr[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
  std::cout<<"Completed calculation of partial contraction at site "<<i<<std::endl;
}
