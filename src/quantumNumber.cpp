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
// This is the new, dynamic labeling scheme, where first, starting with the vacuum labels for the 
// leftmost index, the allowed blocks of the next site are constructed - these are those blocks
// that allow for fullfilling the QN constraints - for all sites from the left. 
// This yields an index list, called leftLabel. The same is now done from the right, starting with
// the final QN as initial labels, such that for each site, a set of indices that are reachable from
// the right side is generated, called rightLabel. 
// Now, for each site, those labels are picked that appear in both the left and the right labeling. 
// These are distributed to the bond indices such that none appears more often than in the either left
// or right labeling scheme, guaranteeing the normalizability of blocks. Also, each one is used at least
// once if the bond dimensions is at least equal to the number of labels.
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
// This is the most sophisticated QN labeling scheme I came up with. It guarantees for maximal
// usage of indices and does not show unnormalizable blocks. Also, all required labels appear (given
// a bond dimension of at least the maximal number of required labels) and no label appears more often
// than in the exact labeling. 
//---------------------------------------------------------------------------------------------------//

int quantumNumber::initializeLabelList(int const i, int const direction){
  //direction==1 is RP, direction==0 is LP and -1 is final index
  int minimalLabel, maximalLabel;
  int blockCounter=0;
  int validBlock, cBlock, allowedBlockSize;
  std::complex<int> *cLabel;
  std::complex<int> label;
  std::vector<int> aimIndices;
  std::vector<int> aiIndices;
  std::vector<int> aiIndicesReordered;
  std::vector<int> aimIndicesReordered;
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
      //First, get the possible indices reachable from those of the last site
      gatherBlocks(i-1+direction,aimIndices,qnLabels,direction);
    }
    else{
      //After the left and right label lists have been generated, get both, the list of indices reachable from right and from left and take each one that appears in both
      gatherBlocks(i-1,aimIndices,qnLabelsLP,0);
      gatherBlocks(i,aiIndices,qnLabelsRP,1);
      qnLabels.clear();
      aimIndicesReordered.clear();
      aiIndicesReordered.clear();
      for(int iLP=0;iLP<qnLabelsLP.size();++iLP){
	for(int iRP=0;iRP<qnLabelsRP.size();++iRP){
	  if(qnLabelsRP[iRP]==qnLabelsLP[iLP]){
	    //In this scheme, the labels might be reordered, therefore, the blocksizes have to be reordered, too
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
	  allowedBlockSize=aimIndices[iBlock];
	}
	if(direction==-1){
	  allowedBlockSize=(aimIndicesReordered[iBlock]>aiIndicesReordered[iBlock])?aiIndicesReordered[iBlock]:aimIndicesReordered[iBlock];
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

void quantumNumber::gatherBlocks(int const i, std::vector<int> &aimIndices, std::vector<std::complex<int> > &qnLabels, int const direction){
  //direction==1 is RP, direction==0 is LP
  int isNew, matchBlock;
  int pre, lD;
  std::vector<std::vector<int> > aimIndexTable;
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
  aimIndexTable.resize(qnLabels.size());
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<lD;++aim){
      for(int iBlock=0;iBlock<qnLabels.size();++iBlock){
	if((groupOperation((this->*labelFunction)(i-1+direction,aim),QNLabel(si),pre)==qnLabels[iBlock] && real((this->*labelFunction)(i-1+direction,aim))>-2)){
	  aimIndexTable[iBlock].push_back(aim);
	}
      }
    }
  }
  for(int iBlock=0;iBlock<qnLabels.size();++iBlock){
    aimIndices[iBlock]=aimIndexTable[iBlock].size();
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
  //This was the labeling scheme (up to permuations) if one would not truncate the MPS site matrices at all
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

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::groupOperation(std::complex<int> a, std::complex<int> b, int const pre){
  //Defines the real part as the U(1) part and the imaginary as the Z_2 part of a quantum number
  std::complex<int> result;
  real(result)=real(a)+pre*real(b);
  imag(result)=imag(a)*imag(b);
  return result;
}

//---------------------------------------------------------------------------------------------------//
