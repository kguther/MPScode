#include "basisQNOrderMatrix.h"
#include <iostream>

basisQNOrderMatrix::basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin):
  dimInfo(dimin),
  conservedQNs(conservedQNsin)
{

}

int basisQNOrderMatrix::blockStructure(int const i, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices){
  int const nQNs=(*conservedQNs).size();
  if(nQNs==0){
    return 1;
  }
  std::vector<int> *qnLabels;
  int isNew=1;
  int matchBlock;
  int numBlocks;
  multInt cMultInd;
  qnLabels=new std::vector<int>[nQNs];
  for(int ai=0;ai<dimInfo.locDimR(i);++ai){
    isNew=1;
    for(int j=0;j<qnLabels[0].size();++j){
      for(int iQN=0;iQN<nQNs;++iQN){
	if((*conservedQNs)[iQN].QNLabel(i,ai)!=qnLabels[iQN][j]){
	  break;
	}
	if(iQN==nQNs-1){
	  isNew=0;
	}
      }
    }
    if(isNew){
      for(int iQN=0;iQN<nQNs;++iQN){
	qnLabels[iQN].push_back((*conservedQNs)[iQN].QNLabel(i,ai));
      }
    }
  }
  numBlocks=qnLabels[0].size();
  aiIndices.resize(numBlocks);
  for(int ai=0;ai<dimInfo.locDimR(i);++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      matchBlock=1;
      for(int iQN=0;iQN<nQNs;++iQN){
	if((*conservedQNs)[iQN].QNLabel(i,ai)!=qnLabels[iQN][iBlock]){
	  matchBlock=0;
	}
      }
      if(matchBlock){
	aiIndices[iBlock].push_back(ai);
      }
    }
  }
  siaimIndices.resize(numBlocks);
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int aim=0;aim<dimInfo.locDimL(i);++aim){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	matchBlock=1;
	for(int iQN=0;iQN<nQNs;++iQN){
	  if((*conservedQNs)[iQN].QNLabel(i-1,aim)+(*conservedQNs)[iQN].QNLabel(si)!=qnLabels[iQN][iBlock]){
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
  delete[] qnLabels;
  return 0;
}
