#include "basisQNOrderMatrix.h"
#include <iostream>

basisQNOrderMatrix::basisQNOrderMatrix():
  aiBlockIndicesLP(0),
  siaimBlockIndicesLP(0),
  aimBlockIndicesRP(0),
  siaiBlockIndicesRP(0),
  aimBlockIndicesSplit(0),
  siBlockIndicesSplit(0),
  siBlockIndicesSplitFixedaim(0),
  aiBlockIndicesSplit(0)
{}

//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin):
  dimInfo(dimin),
  conservedQNs(conservedQNsin),
  aiBlockIndicesLP(0),
  siaimBlockIndicesLP(0),
  aimBlockIndicesRP(0),
  siaiBlockIndicesRP(0),
  aimBlockIndicesSplit(0),
  siBlockIndicesSplit(0),
  siBlockIndicesSplitFixedaim(0),
  aiBlockIndicesSplit(0)
{}

//---------------------------------------------------------------------------------------------------//

basisQNOrderMatrix::~basisQNOrderMatrix(){
  deleteTables();
}

//---------------------------------------------------------------------------------------------------//

void basisQNOrderMatrix::deleteTables(){
  delete[] aiBlockIndicesLP;
  delete[] siaimBlockIndicesLP;
  delete[] siaiBlockIndicesRP;
  delete[] aimBlockIndicesRP;
  delete[] aimBlockIndicesSplit;
  delete[] siBlockIndicesSplit;
  delete[] siBlockIndicesSplitFixedaim;
  delete[] aiBlockIndicesSplit;
}

//---------------------------------------------------------------------------------------------------//

void basisQNOrderMatrix::initialize(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin){
  dimInfo=dimin;
  conservedQNs=conservedQNsin;
}

//---------------------------------------------------------------------------------------------------//

void basisQNOrderMatrix::generateQNIndexTables(){
  int L=dimInfo.L();
  deleteTables();
  aiBlockIndicesLP=new std::vector<std::vector<int> >[L];
  siaimBlockIndicesLP=new std::vector<std::vector<multInt> >[L];
  aimBlockIndicesRP=new std::vector<std::vector<int> >[L];
  siaiBlockIndicesRP=new std::vector<std::vector<multInt> >[L];
  siBlockIndicesSplit=new std::vector<int>[L];
  aimBlockIndicesSplit=new std::vector<int>[L];
  aiBlockIndicesSplit=new std::vector<int>[L];
  siBlockIndicesSplitFixedaim=new std::vector<std::vector<int> >[L];
  for(int i=0;i<dimInfo.L();++i){
    blockStructure(i,0,aiBlockIndicesLP[i],siaimBlockIndicesLP[i]);
    blockStructure(i,1,aimBlockIndicesRP[i],siaiBlockIndicesRP[i]);
    splitIndexTables(i);
    if(i==17 && 1){
      /*
      std::cout<<"Left labels:\n";
      for(int aim=0;aim<dimInfo.locDimL(i+1);++aim){
	std::cout<<aim<<" with label "<<(*conservedQNs)[0].QNLabel(i,aim)<<std::endl;
      }
      */
      for(int iBlock=0;iBlock<numBlocksRP(i);++iBlock){
	std::cout<<"Right indices: "<<std::endl;
	for(int j=0;j<rBlockSizeRP(i,iBlock);++j){
	  std::cout<<aiBlockIndexRP(i,iBlock,j)<<" with label ("<<(*conservedQNs)[0].QNLabel(i,aiBlockIndexRP(i,iBlock,j))<<","<<(*conservedQNs)[1].QNLabel(i,aiBlockIndexRP(i,iBlock,j))<<")\t"<<siBlockIndexRP(i,iBlock,j)<<" with label ("<<(*conservedQNs)[0].QNLabel(siBlockIndexRP(i,iBlock,j))<<","<<(*conservedQNs)[1].QNLabel(siBlockIndexRP(i,iBlock,j))<<")"<<std::endl;
	}
	std::cout<<"Left indices: \n";
	for(int j=0;j<lBlockSizeRP(i,iBlock);++j){
	  std::cout<<aimBlockIndexRP(i,iBlock,j)<<" with label ("<<(*conservedQNs)[0].QNLabel(i-1,aimBlockIndexRP(i,iBlock,j))<<","<<(*conservedQNs)[1].QNLabel(i-1,aimBlockIndexRP(i,iBlock,j))<<")"<<std::endl;
	}
      }
      /*
      for(int iBlock=0;iBlock<numBlocksLP(i);++iBlock){
	std::cout<<"Left indices: "<<std::endl;
	for(int j=0;j<lBlockSizeLP(i,iBlock);++j){
	  std::cout<<aimBlockIndexLP(i,iBlock,j)<<" with label "<<(*conservedQNs)[0].QNLabel(i-1,aimBlockIndexLP(i,iBlock,j))<<"\t"<<siBlockIndexLP(i,iBlock,j)<<" with label "<<(*conservedQNs)[0].QNLabel(siBlockIndexLP(i,iBlock,j))<<std::endl;
	}
	std::cout<<"Right indices: \n";
	for(int j=0;j<rBlockSizeLP(i,iBlock);++j){
	  std::cout<<aiBlockIndexLP(i,iBlock,j)<<" with label "<<(*conservedQNs)[0].QNLabel(i,aiBlockIndexLP(i,iBlock,j))<<std::endl;
	}
      }
      */
      /*
      std::cout<<"aim indices split: \n";
      for(int j=0;j<aimBlockSizeSplit(i,0);++j){
	std::cout<<aimBlockIndexSplit(i,0,j);
	for(int k=0;k<siBlockSizeSplitFixedaim(i,0,j);++k){
	  std::cout<<"\t"<<siBlockIndexSplitFixedaim(i,0,j,k);
	}
	std::cout<<std::endl;
      }*/
      exit(1);
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// The blockStructure function goes through the indices of a given matrix and determines those index
// triplets which fullfill the QN constraint, i.e. which belong to one of the blocks. These are stored
// together with information on their block and can then be accessed by the index functions.
//---------------------------------------------------------------------------------------------------//

int basisQNOrderMatrix::blockStructure(int const i, int const direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices){
  int const nQNs=(*conservedQNs).size();
  if(nQNs==0){
    return 1;
  }
  std::vector<int> *qnLabels;
  int isNew=1;
  int lDsingle;
  int lDpaired;
  int pre;
  int qnConstraint;
  int matchBlock;
  int numBlocks;
  multInt cMultInd;
  if(direction==0){
    lDsingle=dimInfo.locDimR(i);
    lDpaired=dimInfo.locDimL(i);
    pre=1;
  }
  if(direction==1){
    lDsingle=dimInfo.locDimL(i);
    lDpaired=dimInfo.locDimR(i);
    pre=-1;
  }
  qnLabels=new std::vector<int>[nQNs];
  for(int ai=0;ai<lDsingle;++ai){
    isNew=1;
    for(int iBlock=0;iBlock<qnLabels[0].size();++iBlock){
      for(int iQN=0;iQN<nQNs;++iQN){
	if((*conservedQNs)[iQN].QNLabel(i-direction,ai)!=qnLabels[iQN][iBlock]){
	  break;
	}
	if(iQN==nQNs-1){
	  isNew=0;
	}
      }
    }
    if(isNew){
      for(int iQN=0;iQN<nQNs;++iQN){
	qnLabels[iQN].push_back((*conservedQNs)[iQN].QNLabel(i-direction,ai));
      }
    }
  }
  numBlocks=qnLabels[0].size();
  siaimIndices.resize(numBlocks);
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<lDpaired;++aim){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	matchBlock=1;
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(qnCriterium(iQN,i,aim,si,direction,pre)!=qnLabels[iQN][iBlock] || qnLabels[iQN][iBlock]<-2){
	    matchBlock=0;
	  }
	}
	if(matchBlock){
	  cMultInd.si=si;
	  cMultInd.aim=aim;
	  siaimIndices[iBlock].push_back(cMultInd);
	}
      }
    }
  }
  aiIndices.resize(numBlocks);
  for(int ai=0;ai<lDsingle;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      matchBlock=1;
      for(int iQN=0;iQN<nQNs;++iQN){
	if((*conservedQNs)[iQN].QNLabel(i-direction,ai)!=qnLabels[iQN][iBlock]){
	  matchBlock=0;
	}
      }
      if(matchBlock){
	aiIndices[iBlock].push_back(ai);
      }
    }
  }
  delete[] qnLabels;
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void basisQNOrderMatrix::splitIndexTables(int const i){
  int newIndex;
  for(int iBlock=0;iBlock<aiBlockIndicesLP[i].size();++iBlock){
    for(int k=0;k<lBlockSizeLP(i,iBlock);++k){
      newIndex=1;
      for(int splitIndex=0;splitIndex<aimBlockIndicesSplit[i].size();++splitIndex){
	if(aimBlockIndexLP(i,iBlock,k)==aimBlockIndicesSplit[i][splitIndex]){
	  newIndex=0;
	}
      }
      if(newIndex){
	aimBlockIndicesSplit[i].push_back(aimBlockIndexLP(i,iBlock,k));
      }
    }
  }
  siBlockIndicesSplitFixedaim[i].resize(aimBlockSizeSplit(i,0));
  for(int iBlock=0;iBlock<aiBlockIndicesLP[i].size();++iBlock){
    for(int k=0;k<aimBlockSizeSplit(i,iBlock);++k){
      for(int kp=0;kp<lBlockSizeLP(i,iBlock);++kp){
	newIndex=1;
	if(aimBlockIndexLP(i,iBlock,kp)!=aimBlockIndexSplit(i,iBlock,k)){
	  newIndex=0;
	}
	if(newIndex){
	  for(int splitIndex=0;splitIndex<siBlockIndicesSplitFixedaim[i][k].size();++splitIndex){
	    if(siBlockIndexLP(i,iBlock,kp)==siBlockIndicesSplitFixedaim[i][k][splitIndex]){
	      newIndex=0;
	    }
	  }
	}
	if(newIndex){
	  siBlockIndicesSplitFixedaim[i][k].push_back(siBlockIndexLP(i,iBlock,kp));
	}
      }
    }
  }
  for(int iBlock=0;iBlock<numBlocksRP(i);++iBlock){
    for(int j=0;j<rBlockSizeRP(i,iBlock);++j){
      newIndex=1;
      for(int splitIndex=0;splitIndex<aiBlockIndicesSplit[i].size();++splitIndex){
	if(aiBlockIndexRP(i,iBlock,j)==aiBlockIndicesSplit[i][splitIndex]){
	  newIndex=0;
	}
      }
      if(newIndex){
	aiBlockIndicesSplit[i].push_back(aiBlockIndexRP(i,iBlock,j));
      }
      newIndex=1;
      for(int splitIndex=0;splitIndex<siBlockIndicesSplit[i].size();++splitIndex){
	if(siBlockIndexRP(i,iBlock,j)==siBlockIndicesSplit[i][splitIndex]){
	  newIndex=0;
	}
      }
      if(newIndex){
	siBlockIndicesSplit[i].push_back(siBlockIndexRP(i,iBlock,j));
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

int basisQNOrderMatrix::qnCriterium(int const iQN, int const i, int const aim, int const si, int const direction, int const pre){
  if((*conservedQNs)[iQN].parityType()){
    return (*conservedQNs)[iQN].QNLabel(i-1+direction,aim)*(*conservedQNs)[iQN].QNLabel(si);
  }
  return (*conservedQNs)[iQN].QNLabel(i-1+direction,aim)+pre*(*conservedQNs)[iQN].QNLabel(si);
}
