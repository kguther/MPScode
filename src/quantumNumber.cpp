#include <iostream>
#include <cstdlib>
#include "quantumNumber.h"
#include "math.h"

quantumNumber::quantumNumber(){
  indexLabel=0;
  leftLabel=0;
  rightLabel=0;
}

quantumNumber::~quantumNumber(){
  delete[] indexLabel;
  delete[] leftLabel;
  delete[] rightLabel;
}

void quantumNumber::initialize(dimensionTable &dimInfoin, int const Nin, int *QNlocin){
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
  initializeLabelList();
}

int quantumNumber::QNLabel(int const si){
  return QNloc[si];
}

int quantumNumber::QNLowerCheck(int i, int aim){
  if(QNLabel(i,aim)-(i+1)*QNlocMin>=0 && QNLabel(i,aim)-(i+1)*QNlocMax<=0){
    return 0;
  }
  return 1;
}

int quantumNumber::QNUpperCheck(int i, int ai){
  if(QNLabel(i,ai)+(dimInfo.L()-(i))*QNlocMax>=N && QNLabel(i,ai)+(dimInfo.L()-(i))*QNlocMin<=N){
    return 0;
  }
  return 1;
}

int quantumNumber::qnCriterium(int const i, int const si, int const ai, int const aim){
  int qnCriteriumViolation=qnConstraint(i,si,ai,aim);
  if(qnCriteriumViolation || QNUpperCheck(i,ai) || QNLowerCheck(i-1,aim)){
    return 1;
  }
  return 0;
}

int quantumNumber::qnConstraint(int const i, int const si, int const ai, int const aim){
  if(QNLabel(i,ai)-QNLabel(i-1,aim)-QNLabel(si)){
    return 1;
  }
  return 0;
}

int quantumNumber::exactLabel(int const i, int const ai){
  int aiReduced=ai;
  int QNSum=0;
  int sigma;
  if(i==-1){
    return 0;
  }
  for(int j=i;j>=0;--j){
    sigma=aiReduced/pow(dimInfo.d(),j);
    QNSum+=QNLabel(sigma);
    aiReduced-=sigma*pow(dimInfo.d(),j);
  }
  return QNSum;
}

void quantumNumber::setParameterD(int const Dnew){
  dimInfo.setParameterD(Dnew);
  initializeLabelList();
}


void quantumNumber::initializeLabelList(){
  //This is for spin quantum numbers
  int const Lred=dimInfo.L();
  leftLabel=new int[Lred*Lred+1];
  rightLabel=new int[Lred*Lred+1];
  for(int i=0;i<=dimInfo.L();++i){
    if(i<dimInfo.L()/2){
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
  if((i+1)<dimInfo.L()/2 && ai<=(i+1)){
    return leftLabel[ai+(i+1)*dimInfo.L()];
  }
  if(ai<=(dimInfo.L()-(i+1)) && (i+1)>=dimInfo.L()/2){
    return rightLabel[ai+(i+1)*dimInfo.L()];
  }
  return dimInfo.L()*dimInfo.L();
}

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
