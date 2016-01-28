#include <iostream>
#include <cstdlib>
#include "quantumNumber.h"
#include "math.h"

quantumNumber::quantumNumber():
  leftLabel(0),
  rightLabel(0),
  indexLabel(0)
{}

//---------------------------------------------------------------------------------------------------//

quantumNumber::~quantumNumber(){
  delete[] leftLabel;
  delete[] rightLabel;
  delete[] indexLabel;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initialize(dimensionTable &dimInfoin, std::complex<int> const Nin, std::complex<int> *QNlocin){
  int violation;
  N=Nin;
  dimInfo=dimInfoin;
  QNloc=QNlocin;
  delete[] leftLabel;
  delete[] rightLabel;
  delete[] indexLabel;
  leftLabel=new std::complex<int>[dimInfo.D()*(dimInfo.L()+1)];
  rightLabel=new std::complex<int>[dimInfo.D()*(dimInfo.L()+1)];
  indexLabel=new std::complex<int>[dimInfo.D()*(dimInfo.L()+1)];
  initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabel(int const si){
  return QNloc[si];
}


//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabel(int const i, int const ai){
  return indexLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabelLP(int const i, int const ai){
  return leftLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabelRP(int const i, int const ai){
  return rightLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::setParameterD(int const Dnew){
  dimInfo.setParameterD(Dnew);
  initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//
// For our specific system, we use a dedicated minimal labeling, which is hardcoded. Might be improved
// later on, it might even be neccessary to do so, but for now, it is enough. Plus, since it is minimal,
// the number of nonzero matrix elements of the MPO matrices is bound by 4*D, which is really great.
//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelList(){
  initializeLabelListLP();
  initializeLabelListRP();
  primaryIndices.resize(dimInfo.L()+1);
  for(int i=0;i<=dimInfo.L();++i){
    initializeLabelList(i);
  }
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelListLP(){
  for(int i=0;i<=dimInfo.L();++i){
    initializeLabelList(i,0);
  }
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::initializeLabelListRP(){
  for(int i=dimInfo.L();i>=0;--i){
    initializeLabelList(i,1);
  }
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::initializeLabelList(int const i, int const direction){
  //direction==1 is RP, direction==0 is LP and -1 is final index
  int minimalLabel, maximalLabel;
  int blockCounter=0;
  int validBlock, cBlock, allowedBlockSize;
  std::complex<int> *cLabel;
  std::complex<int> label;
  std::vector<std::vector<int> > aimIndices;
  std::vector<std::vector<int> > aiIndices;
  std::vector<std::vector<int> > aiIndicesReordered;
  std::vector<std::vector<int> > aimIndicesReordered;
  std::vector<std::complex<int> > qnLabels;
  std::vector<std::complex<int> > qnLabelsRP;
  std::vector<std::complex<int> > qnLabelsLP;
  std::vector<std::complex<int> > validQNLabels;
  std::vector<int> blockOccupations;
  std::vector<int> maxBlockSizes;
  minimalLabel=(0>real(N)-2*(dimInfo.L()-i))?0:real(N)-2*(dimInfo.L()-i);
  maximalLabel=(2*i>real(N))?real(N):2*i;
  if(direction==1){
    cLabel=rightLabel;
  }
  if(direction==0){
    cLabel=leftLabel;
  }
  if(direction==-1){
    cLabel=indexLabel;
    primaryIndices[i].clear();
  }
  if(i!=0 && i!=dimInfo.L()){
    if(direction!=-1){
      gatherBlocks(i-1+direction,aimIndices,qnLabels,direction);
    }
    else{
      gatherBlocks(i-1,aimIndices,qnLabelsLP,0);
      gatherBlocks(i,aiIndices,qnLabelsRP,1);
      qnLabels.clear();
      aimIndicesReordered.clear();
      aiIndicesReordered.clear();
      for(int iLP=0;iLP<qnLabelsLP.size();++iLP){
	for(int iRP=0;iRP<qnLabelsRP.size();++iRP){
	  if(qnLabelsRP[iRP]==qnLabelsLP[iLP]){
	    aimIndicesReordered.push_back(aimIndices[iLP]);
	    aiIndicesReordered.push_back(aiIndices[iRP]);
	    qnLabels.push_back(qnLabelsLP[iLP]);
	  }
	}
      }
    }
    for(int iBlock=0;iBlock<qnLabels.size();++iBlock){
      validBlock=1;
      if(real(qnLabels[iBlock])>maximalLabel || real(qnLabels[iBlock])<minimalLabel){
	validBlock=0;
      }
      if(real(qnLabels[iBlock])==0 && real(qnLabels[iBlock])==real(N)-2*(dimInfo.L()-i) && integerParity(dimInfo.L()-i)!=imag(N)){
	validBlock=0;
      }
      if(real(qnLabels[iBlock])==real(N) && real(qnLabels[iBlock])==2*i && integerParity(i)!=imag(N)){
	validBlock=0;
      }
      if((real(qnLabels[iBlock])==0 && imag(qnLabels[iBlock])!=1) || (real(qnLabels[iBlock])==real(N)-2*(dimInfo.L()-i) && imag(qnLabels[iBlock])!=integerParity(dimInfo.L()-i)*imag(N))){
	validBlock=0;
      }
      if((real(qnLabels[iBlock])==2*i && imag(qnLabels[iBlock])!=integerParity(i)) || (real(qnLabels[iBlock])==real(N) && imag(qnLabels[iBlock])!=imag(N))){
	validBlock=0;
      }
      if(validBlock){
	if(direction!=-1){
	  allowedBlockSize=aimIndices[iBlock].size();
	}
	if(direction==-1){
	  allowedBlockSize=(aimIndicesReordered[iBlock].size()>aiIndicesReordered[iBlock].size())?aiIndicesReordered[iBlock].size():aimIndicesReordered[iBlock].size();
	}
	validQNLabels.push_back(qnLabels[iBlock]);
	maxBlockSizes.push_back(allowedBlockSize);
      }
    }
    blockOccupations.resize(validQNLabels.size());
    for(int j=0;j<validQNLabels.size();++j){
      blockOccupations[j]=0;
    }
    for(int aim=0;aim<dimInfo.locDimL(i);++aim){
      for(int iBlock=0;iBlock<validQNLabels.size();++iBlock){
	cBlock=(blockCounter+iBlock)%(validQNLabels.size());
	if(blockOccupations[cBlock]<maxBlockSizes[cBlock]){
	  cLabel[aim+i*dimInfo.D()]=validQNLabels[cBlock];
	  if(blockOccupations[cBlock]==0 && direction==-1){
	    primaryIndices[i].push_back(aim);
	  }
	  ++(blockOccupations[cBlock]);
	  break;
	}
	if(iBlock==validQNLabels.size()-1){
	  cLabel[aim+i*dimInfo.D()]=std::complex<int>(-100,1);
	}
      }
      blockCounter=(blockCounter+1)%(validQNLabels.size());
    }
  }
  if(i==0){
    imag(label)=1;
    real(label)=0;
    cLabel[i*dimInfo.D()]=label;
    if(direction==-1){
      primaryIndices[i].push_back(0);
    }
  }
  if(i==dimInfo.L()){
    cLabel[i*dimInfo.D()]=N;
    if(direction==-1){
      primaryIndices[i].push_back(0);
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// The next function defines the QN Labels in the following way: Every particle number that allows for
// reaching the left and rigth vacuum QNs appears twice, once with each parity. Only the maximal and
// minimal particle numbers at a site that are possible appear once, if they have the right parity
// and not at all if they have the wrong parity.
//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::truncLabel(int const i, int const ai){
  /*
  int minimalLabel, maximalLabel, labelRange;
  int aux, treshold, deltaTreshold;
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
  */
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::gatherBlocks(int const i, std::vector<std::vector<int> > &aimIndices, std::vector<std::complex<int> > &qnLabels, int const direction){
  //direction==1 is RP, direction==0 is LP
  int isNew, matchBlock;
  int pre, lD;
  qnLabels.clear();
  aimIndices.clear();
  std::complex<int> (quantumNumber::*labelFunction)(int const i, int const ai);
  if(direction==1){
    labelFunction=&quantumNumber::QNLabelRP;
    pre=-1;
    lD=dimInfo.locDimR(i);
  }
  if(direction==0){
    labelFunction=&quantumNumber::QNLabelLP;
    pre=1;
    lD=dimInfo.locDimL(i);
  }
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<lD;++aim){
      isNew=1;
      for(int iBlock=0;iBlock<qnLabels.size();++iBlock){
	if(groupOperation((this->*labelFunction)(i-1+direction,aim),QNLabel(si),pre)==qnLabels[iBlock]){
	  isNew=0;
	}
      }
      if(isNew){
	qnLabels.push_back(groupOperation((this->*labelFunction)(i-1+direction,aim),QNLabel(si),pre));
      }
    }
  }
  aimIndices.resize(qnLabels.size());
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<lD;++aim){
      for(int iBlock=0;iBlock<qnLabels.size();++iBlock){
	if((groupOperation((this->*labelFunction)(i-1+direction,aim),QNLabel(si),pre)==qnLabels[iBlock] && real((this->*labelFunction)(i-1+direction,aim))>-2)){
	  aimIndices[iBlock].push_back(aim);
	}
      }
    }
  }
}


//---------------------------------------------------------------------------------------------------//

int quantumNumber::primaryIndex(int const i, int const ai){
  for(int aim=0;aim<primaryIndices[i+1].size();++aim){
    if(ai==primaryIndices[i+1][aim]){
      return 1;
    }
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

std::complex<int> quantumNumber::groupOperation(std::complex<int> a, std::complex<int> b, int const pre){
  std::complex<int> result;
  real(result)=real(a)+pre*real(b);
  imag(result)=imag(a)*imag(b);
  return result;
}

//---------------------------------------------------------------------------------------------------//
