#include <iostream>
#include <cstdlib>
#include "quantumNumber.h"
#include "math.h"

quantumNumber::quantumNumber(){
  indexLabel=0;
  leftLabel=0;
  rightLabel=0;
}

//---------------------------------------------------------------------------------------------------//

quantumNumber::~quantumNumber(){
  delete[] indexLabel;
  delete[] leftLabel;
  delete[] rightLabel;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initialize(dimensionTable &dimInfoin, int const Nin, int *QNlocin, int mult){
  int violation;
  N=Nin;
  dimInfo=dimInfoin;
  QNloc=QNlocin;
  QNlocMax=QNloc[0];
  QNlocMin=QNloc[0];
  for(int si=1;si<dimInfo.d();++si){
    if(QNloc[si]>QNlocMax){
      QNlocMax=QNloc[si];
    }
    if(QNloc[si]<QNlocMin){
      QNlocMin=QNloc[si];
    }
  }
  iLRSwap=dimInfo.L()/2;
  parityNumber=mult;
  initializeLabelList();
  for(int i=dimInfo.L()-1;i>0;--i){
    violation=0;
    for(int ai=0;ai<dimInfo.locDimR(i);++ai){
      if(rightLabel[ai+i*dimInfo.D()]==0){
	++violation;
      }
    }
    if(violation>1){
      iLRSwap=i+1;
      break;
    }
  }      
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::QNLabel(int const si){
  return QNloc[si];
}

//---------------------------------------------------------------------------------------------------//



//---------------------------------------------------------------------------------------------------//

int quantumNumber::groupOperation(int const label, int const labelp){
  if(parityNumber){
    return label*labelp;
  }
  return label+labelp;
}

int quantumNumber::groupInverse(int const label){
  if(parityNumber){
    return label;
  }
  return -label;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::QNLowerCheck(int i, int aim){
  if(QNLabel(i,aim)-(i+1)*QNlocMin>=0 && QNLabel(i,aim)-(i+1)*QNlocMax<=0){
    return 0;
  }
  return 1;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::QNUpperCheck(int i, int ai){
  if(QNLabel(i,ai)+(dimInfo.L()-(i))*QNlocMax>=N && QNLabel(i,ai)+(dimInfo.L()-(i))*QNlocMin<=N){
    return 0;
  }
  return 1;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::qnCriterium(int const i, int const si, int const ai, int const aim){
  int qnCriteriumViolation=qnConstraint(i,si,ai,aim);
  if(qnCriteriumViolation || QNUpperCheck(i,ai) || QNLowerCheck(i-1,aim)){
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::qnConstraint(int const i, int const si, int const ai, int const aim){
  if(QNLabel(i,ai)!=groupOperation(QNLabel(i-1,aim),QNLabel(si)) && QNLabel(i,ai)>-2){
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::exactLabel(int const i, int const ai){
  int aiReduced=ai;
  int QNSum=0;
  int sigma;
  if(i==-1){
    if(parityNumber){
      return 1;
    }
    return 0;
  }
  for(int j=i;j>=0;--j){
    sigma=aiReduced/pow(dimInfo.d(),j);
    QNSum=groupOperation(QNSum,QNLabel(sigma));
    aiReduced-=sigma*pow(dimInfo.d(),j);
  }
  return QNSum;

}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::setParameterD(int const Dnew){
  dimInfo.setParameterD(Dnew);
  initializeLabelList();
}

/*
void quantumNumber::initializeLabelList(){
  //This is for spin quantum numbers
  int const Lred=dimInfo.L();
  leftLabel=new int[Lred*Lred+1];
  rightLabel=new int[Lred*Lred+1];
  for(int i=0;i<=dimInfo.L();++i){
    if(i<iLRSwap){
      for(int ai=0;ai<=i;++ai){
	leftLabel[ai+i*Lred]=i-2*ai;
      }
    }
    else{
      for(int ai=0;ai<=dimInfo.L()-i;++ai){
	rightLabel[ai+i*Lred]=N-(Lred-i)+2*ai;
      }
    }
  }
}

int quantumNumber::QNLabel(int const i, int const ai){
  if((i+1)<iLRSwap && ai<=(i+1)){
    return leftLabel[ai+(i+1)*dimInfo.L()];
  }
  if(ai<=(dimInfo.L()-(i+1)) && (i+1)>=iLRSwap){
    return rightLabel[ai+(i+1)*dimInfo.L()];
  }
  return dimInfo.L()*dimInfo.L();
}
*/

/*
void quantumNumber::initializeLabelList(){
  //This is for particle number quantum numbers
  int const Lred=dimInfo.L();
  int const Dred=dimInfo.D();
  int lDR;
  indexLabel=new int[Dred*(Lred+1)];
  for(int j=0;j<=Lred;++j){
    lDR=(j<Lred)?dimInfo.locDimR(j-1):1;
    for(int ai=0;ai<lDR;++ai){
      if(j<=Lred/2){
	indexLabel[ai+j*dimInfo.D()]=exactLabel(j-1,ai);
      }
      else{
	indexLabel[ai+j*dimInfo.D()]=N-exactLabel(Lred-j-1,ai);
      }
      std::cout<<"Label for ai= "<<ai<<" at site "<<j<<": "<<indexLabel[ai+j*dimInfo.D()]<<std::endl;
    }
  }
}

int quantumNumber::QNLabel(int const i, int const ai){
  return indexLabel[ai+(i+1)*dimInfo.D()];
}
*/

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelList(){
  delete[] indexLabel;
  delete[] rightLabel;
  delete[] leftLabel;
  rightLabel=new int[dimInfo.D()*(dimInfo.L()+1)];
  leftLabel=new int[dimInfo.D()*(dimInfo.L()+1)];
  indexLabel=new int[dimInfo.D()*(dimInfo.L()+1)];
  int lDR;
  for(int i=0;i<=dimInfo.iCrit();++i){
    for(int ai=0;ai<dimInfo.locDimL(i);++ai){
      leftLabel[ai+i*dimInfo.D()]=truncLabel(i-1,ai);
      rightLabel[ai+i*dimInfo.D()]=groupOperation(N,groupInverse(abs(truncLabel(dimInfo.L()-1-i,ai))));
    }
  }
  for(int i=dimInfo.L()-1;i>=dimInfo.L()-dimInfo.iCrit()-1;--i){
    for(int ai=0;ai<dimInfo.locDimR(i);++ai){
      leftLabel[ai+(i+1)*dimInfo.D()]=truncLabel(i,ai);
      rightLabel[ai+(i+1)*dimInfo.D()]=groupOperation(N,groupInverse(abs(truncLabel(dimInfo.L()-i-2,ai))));
    }
  }
  for(int i=dimInfo.iCrit()+1;i<=dimInfo.L()-dimInfo.iCrit()-1;++i){
    for(int ai=0;ai<dimInfo.locDimL(i);++ai){
      leftLabel[ai+i*dimInfo.D()]=truncLabel(i-1,ai,0,N);
      rightLabel[ai+i*dimInfo.D()]=groupOperation(N,groupInverse(abs(truncLabel(dimInfo.L()-1-i,ai,N,0))));
    }
  }
}

int quantumNumber::truncLabel(int const i, int const ai){
  return truncLabel(i,ai,0,N);
}

int quantumNumber::truncLabel(int const i, int const ai, int const leftVacuum ,int const rightVacuum){
  int minimalLabel, maximalLabel, labelRange;
  int aux, treshold, offset;
  // Use the minimal index twice if minimalLabel!=0 and the maximal twice if maximalLabel==N since these can be reached in more than one way. The other ones are unique and only one index can exist (else the block structure is corrupted. Be careful to consider the apt parity if an index only appears once.
  minimalLabel=(N-2*(dimInfo.L()-i)>0)?N-2*(dimInfo.L()-i-1):0;
  maximalLabel=(N>2*(i+1))?2*(i+1):N;
  labelRange=maximalLabel-minimalLabel;
  offset=(minimalLabel==0)?1:0;
  if(maximalLabel!=2*(i+1)){
    treshold=2*labelRange+2-offset;
  }
  else{
    treshold=2*labelRange+1-offset;
  }
  if(parityNumber){
    if(ai%2 || ai+1==2*labelRange){
      return 1;
    }
    return -1;
  }
  else{
    if(i==-1){
      return 0;
    }
    aux=ai;
    if(aux<treshold){
      return (aux+offset)/2+minimalLabel;
    }
    return -100;
  }
}

int quantumNumber::QNLabel(int const i, int const ai){
  if((i+1)<iLRSwap){
    return leftLabel[ai+(i+1)*dimInfo.D()];
  }
  return rightLabel[ai+(i+1)*dimInfo.D()];
}

