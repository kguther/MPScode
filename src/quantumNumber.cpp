#include <iostream>
#include <cstdlib>
#include "quantumNumber.h"
#include "math.h"

quantumNumber::quantumNumber(){
  leftLabel=0;
}

//---------------------------------------------------------------------------------------------------//

quantumNumber::~quantumNumber(){
  delete[] leftLabel;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initialize(dimensionTable &dimInfoin, int const Nin, int *QNlocin, int Pin, int mult){
  int violation;
  N=Nin;
  dimInfo=dimInfoin;
  QNloc=QNlocin;
  auxiliaryParityNumber=Pin;
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
  parityNumber=mult;
  initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::QNLabel(int const si){
  return QNloc[si];
}


//---------------------------------------------------------------------------------------------------//

int quantumNumber::QNLabel(int const i, int const ai){
  return leftLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::groupOperation(int const label, int const labelp){
  if(parityNumber!=0){
    return label*labelp;
  }
  return label+labelp;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::groupInverse(int const label){
  if(parityNumber){
    return label;
  }
  return -label;
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
  delete[] leftLabel;
  leftLabel=new int[dimInfo.D()*(dimInfo.L()+1)];
  int lDR;
  for(int i=0;i<dimInfo.L();++i){
    for(int ai=0;ai<dimInfo.locDimL(i);++ai){
      leftLabel[ai+i*dimInfo.D()]=truncLabel(i,ai);
    }
  }
  for(int ai=0;ai<dimInfo.locDimR(dimInfo.L()-1);++ai){
    leftLabel[ai+dimInfo.L()*dimInfo.D()]=truncLabel(dimInfo.L(),ai);
  }
}

//---------------------------------------------------------------------------------------------------//
// The next function defines the QN Labels in the following way: Every particle number that allows for
// reaching the left and rigth vacuum QNs appears twice, once with each parity. Only the maximal and
// minimal particle numbers at a site that are possible appear once, if they have the right parity
// and not at all if they have the wrong parity.
//---------------------------------------------------------------------------------------------------//

int quantumNumber::truncLabel(int const i, int const ai){
  int minimalLabel, maximalLabel, labelRange;
  int aux, treshold, offset;
  minimalLabel=(0>N-2*(dimInfo.L()-i))?0:N-2*(dimInfo.L()-i);
  maximalLabel=(2*i>N)?N:2*i;
  labelRange=maximalLabel-minimalLabel;
  offset=1;
  treshold=(2*labelRange>0)?2*labelRange:1;
  if(parityNumber!=0){
    if(i==dimInfo.L()){
      return parityNumber;
    }
    if(ai==0){
      if(minimalLabel==0){
	return 1;
      }
      else{
	return integerParity(dimInfo.L()-i)*parityNumber;
      }
    }
    if(ai==2*labelRange-1){
      if(maximalLabel==2*i){
	return integerParity(i);
      }
      else{
	return parityNumber;
      }
    }
    return integerParity(ai);
  }
  else{
    if(i==dimInfo.L()){
      return N;
    }
    if(i==0){
      return 0;
    }
    aux=ai;
    if(aux<treshold){
      if(ai==0 && 0==N-2*(dimInfo.L()-i) && integerParity(dimInfo.L()-i)!=auxiliaryParityNumber){
	return -100;
      }
      if(ai==treshold-1 && 2*i==N && integerParity(i)!=auxiliaryParityNumber){
	return -100;
      }
      return (aux+offset)/2+minimalLabel;
    }
    return -100;
  }
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::integerParity(int const n){
  if(n%2){
    return -1;
  }
  return 1;
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
