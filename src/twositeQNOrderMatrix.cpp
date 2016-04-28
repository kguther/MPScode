#include "twositeQNOrderMatrix.h"

twositeQNOrderMatrix::twositeQNOrderMatrix()
{}

//---------------------------------------------------------------------------------------------------//

twositeQNOrderMatrix::twositeQNOrderMatrix(int i, dimensionTable const &dimIn, pseudoQuantumNumber *conservedQNsin):
  conservedQNs(conservedQNsin),
  dimInfo(dimIn),
  site(i)
{
  generateQNIndexTable();
}

//---------------------------------------------------------------------------------------------------//

int twositeQNOrderMatrix::generateQNIndexTable(){
  if(nQNs()==0){
    return 1;
  }
  qnLabels.resize(nQNs());
  
  int isNew=1;
  int const ld=dimInfo.locd(site);
  int const ldp=dimInfo.locd(site+1);
  int const lDL=dimInfo.locDimL(site);
  int const lDRR=dimInfo.locDimR(site+1);
  int const ldDL=ld*lDL;
  int const ldDRR=ldp*lDRR;
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDL;++aim){
      isNew=1;
      for(int iBlock=0;iBlock<qnLabels[0].size();++iBlock){
	for(int iQN=0;iQN<nQNs();++iQN){
	  if(qnCriterium(iQN,site-1,aim,si,1)!=qnLabels[iQN][iBlock] && qnCriterium(iQN,site-1,aim,si,1).real()>-1){
	    break;
	  }
	  if(iQN==nQNs()-1){
	    isNew=0;
	  }
	}
      }
      if(isNew){
	for(int iQN=0;iQN<nQNs();++iQN){
	  qnLabels[iQN].push_back(qnCriterium(iQN,site-1,aim,si,1));
	}
      }
    }
  }
  writeIndexTables(site-1,ld,lDL,lBlockIndices);
  writeIndexTables(site+1,ldp,lDRR,rBlockIndices);
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void twositeQNOrderMatrix::writeIndexTables(int i, int ld, int lD, std::vector<std::vector<multInt> > &target){
  int const numBlocks=qnLabels[0].size();
  int matchBlock;
  int const pre=(i==site-1)?1:-1;
  multInt cMultInd;
  target.resize(numBlocks);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lD;++aim){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	matchBlock=1;
	for(int iQN=0;iQN<nQNs();++iQN){
	  if(qnCriterium(iQN,i,aim,si,pre)!=qnLabels[iQN][iBlock]){
	    matchBlock=0;
	  }
	}
	if(matchBlock){
	  cMultInd.si=si;
	  cMultInd.aim=aim;
	  target[iBlock].push_back(cMultInd);
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> twositeQNOrderMatrix::qnCriterium(int iQN, int i, int ai, int si, int pre){
  return conservedQNs->groupOperation(conservedQNs->QNLabel(i,ai),conservedQNs->QNLabel(si),pre);
}
