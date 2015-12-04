#include "quantumNumber.h"

quantumNumber::quantumNumber(){
  leftLabel=0;
  rightLabel=0;
}

quantumNumber::~quantumNumber(){
  delete[] leftLabel;
  delete[] rightLabel;
}

void quantumNumber::initialize(int const din, int const Lin, int const Nin, int *QNlocin){
  d=din;
  L=Lin;
  N=Nin;
  QNloc=QNlocin;
  QNlocMax=QNloc[0];
  QNlocMin=QNloc[0];
  for(int si=1;si<d;++si){
    if(QNloc[si]>QNlocMax){
      QNlocMax=QNloc[si];
    }
    if(QNloc[si]<QNlocMin){
      QNlocMin=QNloc[si];
    }
  }
  initializeLabelList();
}


void quantumNumber::initializeLabelList(){
  //This is for particle number quantum numbers
  int const Lred=L;
  leftLabel=new int[Lred*Lred];
  rightLabel=new int[Lred*Lred];
  for(int i=0;i<=L;++i){
    if(i<L/2){
      for(int ai=0;ai<i;++ai){
	leftLabel[ai+i*Lred]=i-2*ai;
      }
    }
    else{
      for(int ai=0;ai<=L-i;++ai){
	rightLabel[ai+i*Lred]=N-(L-i)+2*ai;
      }
    }
  }
}

int quantumNumber::QNLabel(int const i, int const ai){
  if((i+1)<L/2 && ai<(i+1)){
    return leftLabel[ai+(i+1)*L];
  }
  if(ai<=(L-(i+1))){
    return rightLabel[ai+(i+1)*L];
  }
  return N*L;
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
  if(QNLabel(i,ai)+(L-(i))*QNlocMax>=N && QNLabel(i,ai)+(L-(i))*QNlocMin<=N){
    return 0;
  }
  return 1;
}
