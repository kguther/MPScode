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

void quantumNumber::initialize(dimensionTable &dimInfoin, std::complex<int> const Nin, std::complex<int> *QNlocin){
  int violation;
  N=Nin;
  dimInfo=dimInfoin;
  QNloc=QNlocin;
  delete[] leftLabel;
  leftLabel=new std::complex<int>[dimInfo.D()*(dimInfo.L()+1)];
  initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabel(int const si){
  return QNloc[si];
}


//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabel(int const i, int const ai){
  return leftLabel[ai+(i+1)*dimInfo.D()];
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
// For our specific system, we use a dedicated minimal labeling, which is hardcoded. Might be improved
// later on, it might even be neccessary to do so, but for now, it is enough. Plus, since it is minimal,
// the number of nonzero matrix elements of the MPO matrices is bound by 4*D, which is really great.
//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelList(){
  for(int i=0;i<dimInfo.L();++i){
    initializeLabelList(i);
  }
  for(int ai=0;ai<dimInfo.locDimR(dimInfo.L()-1);++ai){
    leftLabel[ai+dimInfo.L()*dimInfo.D()]=truncLabel(dimInfo.L(),ai);
  }
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelList(int const i){
  for(int ai=0;ai<dimInfo.locDimL(i);++ai){
    leftLabel[ai+i*dimInfo.D()]=truncLabel(i,ai);
  }
}

//---------------------------------------------------------------------------------------------------//
// The next function defines the QN Labels in the following way: Every particle number that allows for
// reaching the left and rigth vacuum QNs appears twice, once with each parity. Only the maximal and
// minimal particle numbers at a site that are possible appear once, if they have the right parity
// and not at all if they have the wrong parity.
//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::truncLabel(int const i, int const ai){
  int minimalLabel, maximalLabel, labelRange;
  int aux, treshold, deltaTreshold;
  std::complex<int> label;
  minimalLabel=(0>real(N)-2*(dimInfo.L()-i))?0:real(N)-2*(dimInfo.L()-i);
  maximalLabel=(2*i>real(N))?real(N):2*i;
  labelRange=maximalLabel-minimalLabel;
  treshold=(2*labelRange>0)?2*labelRange:1;
  imag(label)=integerParity(ai);
  if(i==dimInfo.L()){
    imag(label)=imag(N);
  }
  else{
    if(ai==0){
      if(minimalLabel==0){
	imag(label)=1;
      }
      else{
	imag(label)=integerParity(dimInfo.L()-i)*imag(N);
      }
    }
    else{
      if(ai==2*labelRange-1){
	if(maximalLabel==2*i){
	  imag(label)=integerParity(i);
	}
	else{
	  imag(label)=imag(N);
	}
      }
    }
  }
  if(ai==0){
    std::cout<<treshold*2-2<<std::endl;
  }
  aux=ai;
  real(label)=-100;
  if(aux<treshold){
    if(ai==0 && 0==real(N)-2*(dimInfo.L()-i) && integerParity(dimInfo.L()-i)!=imag(N)){
      real(label)=-100;
    }
    else{
      if(ai==treshold-1 && 2*i==real(N) && integerParity(i)!=imag(N)){
	real(label)=-100;
      }
      else{
	real(label)=(aux+1)/2+minimalLabel;
      }
    }
  }
  else{
  deltaTreshold=2;
  while(deltaTreshold==2){
    aux-=treshold;
    treshold-=deltaTreshold;
    ++minimalLabel;
    if(treshold>0 && aux<treshold){
      real(label)=aux/2+minimalLabel;
    }
    if(deltaTreshold==2){
      deltaTreshold=4;
    }
  }
  }
  if(i==dimInfo.L()){
    real(label)=real(N);
  }
  if(i==0){
    real(label)=0;
  }
  return label;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::gatherBlocks(int const i, std::vector<std::vector<int> > &aimIndices, std::vector<std::vector<int> > &siIndices){
  int isNew;
  std::vector<std::complex<int> > qnLabels;
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<dimInfo.locDimL(i);++aim){
      isNew=1;
      for(int iBlock=0;iBlock<aimIndices.size();++iBlock){
	if(QNLabel(i-1,aim)==qnLabels[iBlock]){
	  isNew=0;
	}
      }
      if(isNew){
	qnLabels.push_back(QNLabel(i,aim));
      }
    }
  }
}


//---------------------------------------------------------------------------------------------------//

int quantumNumber::primaryIndex(int const i, int const ai){
  int const j=i+1;
  if(j==dimInfo.L() || j==0){
    return 1;
  }
  int minimalLabel, maximalLabel, labelRange;
  int treshold;
  minimalLabel=(0>real(N)-2*(dimInfo.L()-j))?0:real(N)-2*(dimInfo.L()-j);
  maximalLabel=(2*j>real(N))?real(N):2*j;
  labelRange=maximalLabel-minimalLabel;
  treshold=(2*labelRange>0)?2*labelRange:1;
  if(ai<treshold){
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::integerParity(int const n){
  if(n%2){
    return -1;
  }
  return 1;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::qnConstraint(int const i, int const si, int const ai, int const aim){
  if(QNLabel(i,ai)!=groupOperation(QNLabel(i-1,aim),QNLabel(si)) && real(QNLabel(i,ai))>-2){
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::exactLabel(int const i, int const ai){
  int aiReduced=ai;
  std::complex<int> QNSum=0;
  int sigma;
  for(int j=i;j>=0;--j){
    sigma=aiReduced/pow(dimInfo.d(),j);
    QNSum=groupOperation(QNSum,QNLabel(sigma));
    aiReduced-=sigma*pow(dimInfo.d(),j);
  }
  return QNSum;
}

std::complex<int> quantumNumber::groupOperation(std::complex<int> a, std::complex<int> b){
  std::complex<int> result;
  real(result)=real(a)+real(b);
  imag(result)=imag(a)*imag(b);
  return result;
}
