#include "basisQNOrderMatrix.h"
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// There is not really much to see here, a bunch of vectors is supplied which store the indices of 
// the blocks of the MPS matrices.
//---------------------------------------------------------------------------------------------------//

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
  int cumulativeBlockSize;
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
    if(i==0 && 0){
      // This part is used to test QN labeling schemes for their useability. It prints out the block indices and their QN labels.
      std::cout<<"Right labels:\n";
      for(int aim=0;aim<dimInfo.locDimL(i+1);++aim){
	std::cout<<aim<<" with label ("<<(*conservedQNs)[0].QNLabel(i,aim)<<","<<(*conservedQNs)[1].QNLabel(i,aim)<<")"<<std::endl;
      }
      /*
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
      */
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
    cumulativeBlockSize=0;
    for(int iBlock=0;iBlock<numBlocksLP(i);++iBlock){
      cumulativeBlockSize+=lBlockSizeLP(i,iBlock);
    }
    if(cumulativeBlockSize==0){
      // If the cumulativeBlockSize is zero, then there are no indices fullfilling the QN constraint on this site. That means the right vacuum QN is invalid, for example N=80 for a chain of length 20.
      std::cout<<"Critical error: Invalid quantum number. Terminating process.\n";
      exit(2);
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
  //There are two modes for constructing the Blocks: left pairing (LP) and right pairing (RP). The mode indicates whether si is paired with aim (LP) or with ai (RP) to form a new block index. The other index remains unpaired and determines the QN of the block
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
  //First, determine the number of blocks and their QNs by counting over the different QNs of the unpaired index.
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
  //Now, for each block, collect those pairs of indices from the paired indices that yield the QN of the block
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
  //And, for each block, collect all unpaired indices with the QN of that block. For a minimal labeling, this is not necessary, but we want to keep it general.
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
// Outdated function that splits the paired indices in some weird ways required to build a matrix
// multiplication based on the block structure. It turned out that this is very inefficient, so it is
// not used anymore.
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
