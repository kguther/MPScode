#include "basisQNOrderMatrix.h"
#include <iostream>

basisQNOrderMatrix::basisQNOrderMatrix():
  aiBlockIndicesLP(0),
  siaimBlockIndicesLP(0),
  aimBlockIndicesRP(0),
  siaiBlockIndicesRP(0)
{}

basisQNOrderMatrix::basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin):
  dimInfo(dimin),
  conservedQNs(conservedQNsin),
  aiBlockIndicesLP(0),
  siaimBlockIndicesLP(0),
  aimBlockIndicesRP(0),
  siaiBlockIndicesRP(0)
{}

basisQNOrderMatrix::~basisQNOrderMatrix(){
  delete[] aiBlockIndicesLP;
  delete[] siaimBlockIndicesLP;
  delete[] siaiBlockIndicesRP;
  delete[] aimBlockIndicesRP;
}

void basisQNOrderMatrix::initialize(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin){
  dimInfo=dimin;
  conservedQNs=conservedQNsin;
}

void basisQNOrderMatrix::generateQNIndexTables(){
  int L=dimInfo.L();
  aiBlockIndicesLP=new std::vector<std::vector<int> >[L];
  siaimBlockIndicesLP=new std::vector<std::vector<multInt> >[L];
  aimBlockIndicesRP=new std::vector<std::vector<int> >[L];
  siaiBlockIndicesRP=new std::vector<std::vector<multInt> >[L];
  for(int i=0;i<dimInfo.L();++i){
    blockStructure(i,0,aiBlockIndicesLP[i],siaimBlockIndicesLP[i]);
    blockStructure(i,1,aimBlockIndicesRP[i],siaiBlockIndicesRP[i]);
  }
}

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
    /*for(int si=0;si<dimInfo.locd(i);++si){
      for(int aim=0;aim<lDpaired;++aim){
	qnConstraint=nQNs;
	for(int iQN=0;iQN<nQNs;++iQN){
	  if((*conservedQNs)[iQN].qnConstraint(i,si,ai,aim)==0){
	    qnConstraint-=1;
	  }
	}
      }
    }
    */
    qnConstraint=0;
    if(isNew && !qnConstraint){
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
	  if((*conservedQNs)[iQN].QNLabel(i-1+direction,aim)+pre*(*conservedQNs)[iQN].QNLabel(si)!=qnLabels[iQN][iBlock]){
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
