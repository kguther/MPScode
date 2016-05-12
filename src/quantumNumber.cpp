#include <iostream>
#include <cstdlib>
#include "quantumNumber.h"
#include "math.h"

quantumNumber::quantumNumber(){
}

//---------------------------------------------------------------------------------------------------//

quantumNumber::quantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin):
  pseudoQuantumNumber(dimInfoin,Nin,QNlocin),
  failed(0)
{
  if(imag(N)==0){
    N.imag(0);
    for(int m=0;m<QNloc.size();++m){
      (QNloc[m]).imag(0);
    }
  }
  int info;
  info=initializeLabelList();
  if(info){
    std::cout<<"Critical error: Target quantum number cannot be reached.\n";
    std::cout<<"Target quantum number: "<<N<<" with system size "<<dimInfo.L()<<std::endl;
    failed=1;
  }
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabelLP(int i, int ai){
  return leftLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::QNLabelRP(int i, int ai){
  return rightLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::setParameterD(int Dnew){
  dimInfo.setParameterD(Dnew);
  return initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::setParameterL(int Lnew){
  dimInfo.setParameterL(Lnew);
  return initializeLabelList();
}

//---------------------------------------------------------------------------------------------------//
// Functions for easy manipulation of the index labels during runtime. 

// grow() adds sites to reach length L at site i. The new bonds have undefined indices, except for the last one which takes source as indices. Then, the indices of the right part are shifted by
// targetQN-N such that the new lattice has total charge targetQN.

// refine() defines the indices at site i by using the entries of source. It returns -1 if source is of insufficient size. Use to define labels for bond indices introduces with grow()
//---------------------------------------------------------------------------------------------------//

int quantumNumber::grow(int L, int i, std::complex<int> const &targetQN, std::vector<std::complex<int> > const &source){
  int const D=dimInfo.D();
  /*
  std::cout<<"SITE: "<<i<<std::endl;

  for(int i=0;i<dimInfo.L()+1;++i){
    for(int ai=0;ai<D;++ai){
      std::cout<<indexLabel[ai+D*i]<<"\t";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  */
  int const dL=L-dimInfo.L();
  indexLabel.insert(indexLabel.begin()+(i+1)*D,dL*D,invalidQN);
  int const lDL=dimInfo.locDimL(i);
  for(int aim=0;aim<lDL;++aim){
    indexLabel[aim+(i+dL)*D]=source[aim];
  }
  int const deltaN=targetQN.real()-N.real();

  std::cout<<"Old target QN: "<<N.real()<<" New target QN: "<<targetQN.real()<<std::endl;

  for(int j=(i+dL)*D;j<D*(L+1);++j){
    indexLabel[j].real(indexLabel[j].real()+deltaN);
  }
  
  dimInfo.setParameterL(L);
  N=targetQN;
  
  /*
  std::cout<<std::endl;
  for(int i=0;i<dimInfo.L()+1;++i){
    for(int ai=0;ai<D;++ai){
      std::cout<<indexLabel[ai+D*i]<<"\t";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  */
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::refine(int i, std::vector<std::complex<int> > const &source){
  int const lDL=dimInfo.locDimL(i);
  int const D=dimInfo.D();
  if(source.size()<lDL){
    std::cout<<"Invalid input in refinement\n";
    return -1;
  }

  
  for(int aim=0;aim<lDL;++aim){
    std::cout<<"Old label: "<<QNLabel(i-1,aim)<<"\t";
    indexLabel[aim+i*D]=source[aim];
    std::cout<<"New label: "<<QNLabel(i-1,aim)<<"\t"<<"Index: "<<aim<<"\n";
  }
  
  std::cout<<std::endl;
  return 0;
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

int quantumNumber::initializeLabelList(){
  int info;
  leftLabel.resize(dimInfo.D()*(dimInfo.L()+1));
  rightLabel.resize(dimInfo.D()*(dimInfo.L()+1));
  indexLabel=std::vector<std::complex<int> >(dimInfo.D()*(dimInfo.L()+1),invalidQN);
  info=initializeLabelListLP();
  if(info)
    return info;
  initializeLabelListRP();
  if(info)
    return info;
  primaryIndices.resize(dimInfo.L()+1);
  for(int i=0;i<=dimInfo.L();++i){
    info=initializeLabelList(i);
    if(info)
      return info;
  }
  leftLabel.clear();
  rightLabel.clear();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::initializeLabelListLP(){
  int info;
  for(int i=0;i<=dimInfo.L();++i){
    info=initializeLabelList(i,0);
    if(info)
      return info;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::initializeLabelListRP(){
  int info;
  for(int i=dimInfo.L();i>=0;--i){
    info=initializeLabelList(i,1);
    if(info)
      return info;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// This is the most sophisticated QN labeling scheme I came up with. It guarantees for maximal
// usage of indices and does not show unnormalizable blocks. Also, all required labels appear (given
// a bond dimension of at least the maximal number of required labels) and no label appears more often
// than in the exact labeling. 
//---------------------------------------------------------------------------------------------------//

int quantumNumber::initializeLabelList(int i, int direction){
  //direction==1 is RP, direction==0 is LP and -1 is final index
  int validBlock, cBlock, allowedBlockSize;
  std::vector<std::complex<int> > *cLabel;
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
  if(direction==1){
    cLabel=&rightLabel;
  }
  if(direction==0){
    cLabel=&leftLabel;
  }
  if(direction==-1){
    cLabel=&indexLabel;
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
      validBlock=validQN(i,qnLabels[iBlock]);
      //Check if the Block can be meaningful (i.e. can be reached from both sides). This should always be true in the final run, but better check it.
      
      //I did not believe this, but this checks in the warmup are actually required. The scheme does not work without

      if(!validBlock && direction==-1){
	std::cout<<"Invalid label at site "<<i<<std::endl;
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
    //Evenly distribute the bond indices to the available labels - this might not be optimal, but it always results in a valid scheme
    int blockCounter=0;
    for(int aim=0;aim<dimInfo.locDimL(i);++aim){
      for(int iBlock=0;iBlock<validQNLabels.size();++iBlock){
	cBlock=(blockCounter+iBlock)%(validQNLabels.size());
	if(blockOccupations[cBlock]<maxBlockSizes[cBlock]){
	  (*cLabel)[aim+i*dimInfo.D()]=validQNLabels[cBlock];
	  if(blockOccupations[cBlock]==0 && direction==-1){
	    primaryIndices[i].push_back(aim);
	  }
	  ++(blockOccupations[cBlock]);
	  break;
	}
	if(iBlock==validQNLabels.size()-1){
	  (*cLabel)[aim+i*dimInfo.D()]=invalidQN;
	}
      }
      if(validQNLabels.size()==0){
	return -1;
      }
      blockCounter=(blockCounter+1)%(validQNLabels.size());
    }
  }
  if(i==0){
    label.imag(1);
    label.real(0);
    (*cLabel)[i*dimInfo.D()]=label;
    if(direction==-1){
      primaryIndices[i].push_back(0);
    }
  }
  if(i==dimInfo.L()){
    (*cLabel)[i*dimInfo.D()]=N;
    if(direction==-1){
      primaryIndices[i].push_back(0);
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void quantumNumber::gatherBlocks(int i, std::vector<int> &aimIndices, std::vector<std::complex<int> > &qnLabels, int direction){
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
	if((groupOperation((this->*labelFunction)(i-1+direction,aim),QNLabel(si),pre)==qnLabels[iBlock] && !isInvalid((this->*labelFunction)(i-1+direction,aim)))){
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

int quantumNumber::primaryIndex(int i, int ai){
  for(int aim=0;aim<primaryIndices[i+1].size();++aim){
    if(ai==primaryIndices[i+1][aim]){
      return 1;
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int quantumNumber::integerParity(int n) const{
  if(n%2){
    return -1;
  }
  return 1;
}
  

//---------------------------------------------------------------------------------------------------//

int quantumNumber::qnConstraint(int i, int si, int ai, int aim){
  if(QNLabel(i,ai)!=groupOperation(QNLabel(i-1,aim),QNLabel(si)) || isInvalid(QNLabel(i,ai))){
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> quantumNumber::exactLabel(int i, int ai){
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

int quantumNumber::validQN(int i, std::complex<int> const &label) const{
  int validBlock=1;
  int const leftVacuum=0;
  int const maxCharge=2*i;
  int const minimalLabel=(leftVacuum>real(N)-2*(dimInfo.L()-i))?leftVacuum:real(N)-2*(dimInfo.L()-i);
  int const maximalLabel=(maxCharge>real(N))?real(N):maxCharge;
  if(real(label)>maximalLabel || real(label)<minimalLabel){
    validBlock=0;
  }
  if(real(label)==leftVacuum && real(label)==real(N)-2*(dimInfo.L()-i) && integerParity(dimInfo.L()-i)!=imag(N)){
    if(imag(N)){
      validBlock=0;
    }
  }
  if(real(label)==real(N) && real(label)==maxCharge && integerParity(i)!=imag(N)){
    if(imag(N)){
      validBlock=0;
    }
  }
  if((real(label)==leftVacuum && imag(label)!=1) || (real(label)==real(N)-2*(dimInfo.L()-i) && imag(label)!=integerParity(dimInfo.L()-i)*imag(N))){
    if(imag(N)){
      validBlock=0;
    }
  }
  if((real(label)==maxCharge && imag(label)!=integerParity(i)) || (real(label)==real(N) && imag(label)!=imag(N))){
    if(imag(N)){
      validBlock=0;
    }
  }
  return validBlock;
}

