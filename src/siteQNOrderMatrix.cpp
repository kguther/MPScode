#include "siteQNOrderMatrix.h"
#include "exceptionClasses.h"
#include <iostream>

siteQNOrderMatrix::siteQNOrderMatrix(int site, int lDLin, int lDRin, int ldIn, std::vector<pseudoQuantumNumber*> const &conservedQNsin):
  lDL(lDLin),
  lDR(lDRin),
  ld(ldIn),
  i(site),
  conservedQNs(conservedQNsin)
{
  setUpTable();
}

//---------------------------------------------------------------------------------------------------//

siteQNOrderMatrix::siteQNOrderMatrix(int site, int lDLin, int lDRin, int ldIn, std::vector<quantumNumber> *conservedQNsin):
  lDL(lDLin),
  lDR(lDRin),
  ld(ldIn),
  i(site)
{
  loadConservedQNs(conservedQNsin);
  setUpTable();
}

//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::generateFull(siteQNOrderMatrix const &source){
  i=source.i;
  lDL=source.lDL;
  lDR=source.lDR;
  ld=source.ld;
  conservedQNs=source.conservedQNs;
  setUpTableFull();
}

//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::loadConservedQNs(std::vector<quantumNumber> *conservedQNsin){
  conservedQNs.resize(conservedQNsin->size());
  for(int m=0;m<conservedQNs.size();++m){
    conservedQNs[m]=&((*conservedQNsin)[m]);
  }
}

//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::setUpTable(){
  blockStructure(0,aiBlockIndicesLP,siaimBlockIndicesLP);
  blockStructure(1,aimBlockIndicesRP,siaiBlockIndicesRP);
  //Check if there are any non-empty blocks (required). If this is true for LP, it is also for RP
  int cumulativeBlockSize=0;
  for(int iBlock=0;iBlock<numBlocksLP();++iBlock){
    cumulativeBlockSize+=lBlockSizeLP(iBlock)*rBlockSizeLP(iBlock);
  }
  if(cumulativeBlockSize==0){
    //no-throw guarantee if conservedQNs is unrefined
    throw empty_table(i);
  }
}

//---------------------------------------------------------------------------------------------------//
// setUpTableFull() uses the full blockstructure including empty blocks (see below)
//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::setUpTableFull(){
  blockStructureFull(0,aiBlockIndicesLP,siaimBlockIndicesLP);
  blockStructureFull(1,aimBlockIndicesRP,siaiBlockIndicesRP);
}

//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::blockStructure(int direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices){
  int const nQNs=conservedQNs.size();
  if(nQNs==0){
    //ADD EXCEPTION
  }
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
    lDsingle=lDR;
    lDpaired=lDL;
    pre=1;
  }
  if(direction==1){
    lDsingle=lDL;
    lDpaired=lDR;
    pre=-1;
  }
  //First, determine the number of blocks and their QNs by counting over the different QNs of the unpaired index.
  std::vector<std::vector<std::complex<int> > > qnLabels;
  qnLabels.resize(nQNs);
  for(int ai=0;ai<lDsingle;++ai){
    isNew=1;
    for(int iBlock=0;iBlock<qnLabels[0].size();++iBlock){
      for(int iQN=0;iQN<nQNs;++iQN){
	if((conservedQNs[iQN])->QNLabel(i-direction,ai)!=qnLabels[iQN][iBlock]){
	  break;
	}
	if(iQN==nQNs-1){
	  isNew=0;
	}
      }
    }
    if(isNew){
      for(int iQN=0;iQN<nQNs;++iQN){
	qnLabels[iQN].push_back((conservedQNs[iQN])->QNLabel(i-direction,ai));
      }
    }
  }
  //Now, for each block, collect those pairs of indices from the paired indices that yield the QN of the block
  numBlocks=qnLabels[0].size();
  siaimIndices.resize(numBlocks);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDpaired;++aim){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	matchBlock=1;
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(qnCriterium(iQN,aim,si,direction,pre)!=qnLabels[iQN][iBlock] || conservedQNs[iQN]->isInvalid(qnLabels[iQN][iBlock])){
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
	if((conservedQNs[iQN])->QNLabel(i-direction,ai)!=qnLabels[iQN][iBlock]){
	  matchBlock=0;
	}
      }
      if(matchBlock){
	aiIndices[iBlock].push_back(ai);
      }
    }
  }
  if(direction==0){
    qnLabelsLP=qnLabels[0];
  }
  else{
    qnLabelsRP=qnLabels[0];
  }
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> siteQNOrderMatrix::qnCriterium(int iQN, int aim, int si, int direction, int pre){
  std::complex<int> criterium;
  //criterium.real(real((conservedQNs[iQN])->QNLabel(i-1+direction,aim))+pre*real((conservedQNs[iQN])->QNLabel(si)));
  //criterium.imag(imag((conservedQNs[iQN])->QNLabel(i-1+direction,aim))*imag((conservedQNs[iQN])->QNLabel(si)));
  criterium=(conservedQNs[iQN])->groupOperation((conservedQNs[iQN])->QNLabel(i-1+direction,aim),(conservedQNs[iQN])->QNLabel(si),pre);
  return criterium;
}

//---------------------------------------------------------------------------------------------------//
// The validate() function returns -1 if there are non-normalizable blocks and 0 else.
//---------------------------------------------------------------------------------------------------//

int siteQNOrderMatrix::validate()const {
  int lBlockSize, rBlockSize;
  for(int iBlock=0;iBlock<numBlocksLP();++iBlock){
    lBlockSize=lBlockSizeLP(iBlock);
    rBlockSize=rBlockSizeLP(iBlock);
    //If the right basis of some site matrix contains more indices of some label than can be reached from the left basis, the labeling scheme is not valid and leads to non-normalizable blocks.
    if(lBlockSize<rBlockSize && lBlockSize>0){
      return -(i+1);
    }
  }
  for(int iBlock=0;iBlock<numBlocksRP();++iBlock){
    lBlockSize=lBlockSizeRP(iBlock);
    rBlockSize=rBlockSizeRP(iBlock);
    //The same is true for reaching the left basis from the right.
    if(rBlockSize<lBlockSize && rBlockSize>0){
      return i+1;
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// The index tables generated via blockStructureFull() contain all qnLabels reachable from the paired
// side. Thus, all possible labels are contained, but there will most likely be some empty blocks
// for those labels which are reachable, but do not appear on the other side. 
// Therefore, the full block structure is slower and not used for reading access. 
// It is only required when labels are adapted in the enrichment step.
//---------------------------------------------------------------------------------------------------//

void siteQNOrderMatrix::blockStructureFull(int direction, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices){
  int const nQNs=conservedQNs.size();
  if(nQNs==0){
    //ADD EXCEPTION
  }
  int isNew=1;
  int lDsingle;
  int lDpaired;
  int pre;
  int qnConstraint;
  int matchBlock;
  int numBlocks;
  int valid;
  multInt cMultInd;
  //There are two modes for constructing the Blocks: left pairing (LP) and right pairing (RP). The mode indicates whether si is paired with aim (LP) or with ai (RP) to form a new block index. The other index remains unpaired and determines the QN of the block
  if(direction==0){
    lDsingle=lDR;
    lDpaired=lDL;
    pre=1;
  }
  if(direction==1){
    lDsingle=lDL;
    lDpaired=lDR;
    pre=-1;
  }
  //First, determine the number of blocks and their QNs by counting over the different QNs of the unpaired index.
  std::vector<std::vector<std::complex<int> > > qnLabels;
  qnLabels.resize(nQNs);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDpaired;++aim){
      isNew=1;
      valid=1;
      for(int iBlock=0;iBlock<qnLabels[0].size();++iBlock){
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(qnCriterium(iQN,aim,si,direction,pre)!=qnLabels[iQN][iBlock]){
	    break;
	  }
	  if(iQN==nQNs-1){
	    isNew=0;
	  }
	}
      }
      if(isNew){
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(!(conservedQNs[iQN])->validQN(i+1-direction,qnCriterium(iQN,aim,si,direction,pre))){
	    valid=0;
	  }
	  for(int iQN=0;iQN<nQNs;++iQN){
	    if(valid){
	      qnLabels[iQN].push_back(qnCriterium(iQN,aim,si,direction,pre));
	    }
	  }
	}
      }
    }
  }
  //Now, for each block, collect those pairs of indices from the paired indices that yield the QN of the block
  numBlocks=qnLabels[0].size();
  siaimIndices.resize(numBlocks);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDpaired;++aim){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	matchBlock=1;
	for(int iQN=0;iQN<nQNs;++iQN){
	  if(qnCriterium(iQN,aim,si,direction,pre)!=qnLabels[iQN][iBlock] || conservedQNs[iQN]->isInvalid(qnLabels[iQN][iBlock])){
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
	if((conservedQNs[iQN])->QNLabel(i-direction,ai)!=qnLabels[iQN][iBlock]){
	  matchBlock=0;
	}
      }
      if(matchBlock){
	aiIndices[iBlock].push_back(ai);
      }
    }
  }
  if(direction==0){
    qnLabelsLP=qnLabels[0];
  }
  else{
    qnLabelsRP=qnLabels[0];
  }
}

//---------------------------------------------------------------------------------------------------//

void printIndexTable(siteQNOrderMatrix const &a){
  int const numBlocks=a.numBlocksRP();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    std::cout<<"Block with QN "<<a.qnLabelRP(iBlock)<<" of size "<<a.lBlockSizeRP(iBlock)<<"x"<<a.rBlockSizeRP(iBlock)<<std::endl;
  }
}
