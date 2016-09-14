#include "verifyQN.h"
#include "quantumNumber.h"
#include "mpstype.h"
#include <iostream>

//Verifying the QN constraint for testing purposes

int checkQNConstraint(mps &test, int i){
  int const ld=test.getDimInfo().locd(i);
  int const lDL=test.getDimInfo().locDimL(i);
  int const lDR=test.getDimInfo().locDimR(i);
  int failed=0;
  pseudoQuantumNumber *gqn=&test.getConservedQNs()[0];
  mpsEntryType *subMatrix;
  mpsEntryType mEl;
  for(int si=0;si<ld;++si){
    test.subMatrixStart(subMatrix,i,si);
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	mEl=subMatrix[aim+lDL*ai];	
	if(gqn->groupOperation(gqn->QNLabel(i-1,aim),gqn->QNLabel(si))!=gqn->QNLabel(i,ai) && std::abs(mEl)>1e-10){
	  std::cout<<"Matrix element at left labels: "<<gqn->QNLabel(i-1,aim)<<" "<<gqn->QNLabel(si)<<" right label: "<<gqn->QNLabel(i,ai)<<" : "<<mEl<<" with indices "<<"("<<aim<<","<<si<<")"<<" "<<ai<<std::endl;	
	  failed=1;
	}
      }
    }
  }
  return failed;
}

/*
int checkQNConstraint(impBase &test, int i){
  int const ld=test.getDimInfo().locd(i);
  int const lDL=test.getDimInfo().locDimL(i);
  int const lDR=test.getDimInfo().locDimR(i);
  int failed=0;
  pseudoQuantumNumber *gqn=test.getConservedQNs(0);
  mpsEntryType *subMatrix;
  mpsEntryType mEl;
  for(int si=0;si<ld;++si){
    test.subMatrixStart(subMatrix,i,si);
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	mEl=subMatrix[aim+lDL*ai];	
	if(gqn->groupOperation(gqn->QNLabel(i-1,aim),gqn->QNLabel(si))!=gqn->QNLabel(i,ai) && std::abs(mEl)>1e-10){
	  std::cout<<"Matrix element at left labels: "<<gqn->QNLabel(i-1,aim)<<" "<<gqn->QNLabel(si)<<" right label: "<<gqn->QNLabel(i,ai)<<" : "<<mEl<<std::endl;	
	  failed=1;
	}
      }
    }
  }
  return failed;
}

int checkQNConstraint(impBase &test){
  int info=0;
  for(int i=0;i<test.length();++i){
    info=checkQNConstraint(test,i);
    if(info){
      return i+1;
    }
  }
  return 0;
}

int twositeCheck(impBase &test, mpsEntryType *M){
  int i=test.currentSite();
  int const ld=test.getDimInfo().locd(i);
  int const lDL=test.getDimInfo().locDimL(i);
  pseudoQuantumNumber *gqn=test.getConservedQNs(0);
  mpsEntryType mEl;
  for(int si=0;si<ld;++si){
    for(int sip=0;sip<ld;++sip){
      for(int aim=0;aim<lDL;++aim){
	for(int air=0;air<lDL;++air){
	  mEl=M[aim+lDL*air+sip*lDL*lDL+si*ld*lDL*lDL];
	  if(gqn->groupOperation(gqn->QNLabel(i-1,aim),gqn->QNLabel(sip))!=gqn->groupOperation(gqn->QNLabel(i+1,air),gqn->QNLabel(si),-1) && std::abs(mEl)>1e-10){   
	    std::cout<<"Violation at matrix element with labels: "<<gqn->QNLabel(i-1,aim)<<" "<<gqn->QNLabel(sip)<<" right labels: "<<gqn->QNLabel(i+1,air)<<gqn->QNLabel(si)<<" : "<<mEl<<std::endl;
	    return 1;
	  }
	}
      }
    }
  }
  return 0;
}
*/
