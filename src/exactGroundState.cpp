#include "exactGroundState.h"
#include <vector>

exactGroundState::exactGroundState(std::complex<int> N):QNValue(N)
{
  if(imag(QNValue)==0){
    QNValue.imag(1);
  }
}

//---------------------------------------------------------------------------------------------------//

void exactGroundState::writeExactGroundState(mps &target){
  std::complex<int> localQNs[4]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  std::vector<std::complex<int> > localQNsVec;
  for(int m=0;m<4;++m){
    localQNsVec.push_back(localQNs[m]);
  }
  QNsVec.push_back(quantumNumber(target.dimInfo,QNValue,localQNsVec));
  mps proxy(target.dimInfo,QNsVec);
  generateExactState(proxy);
  target.stateArray::mpsCpy(proxy);
}

//---------------------------------------------------------------------------------------------------//

void exactGroundState::generateExactState(mps &target){
    //This is the exact ground state at the critical point for fixed particle number and subchain parity. It turns out that this is a nice guess for the ground state of the perturbed system (for small perturbations).
  int ld, lDL, lDR;
  for(int i=0;i<target.dimInfo.L();++i){
    lDL=target.locDimL(i);
    lDR=target.locDimR(i);
    ld=target.locd(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  target.global_access(i,si,ai,aim)=0;
	}
      }
    }
    int numBlocks, lBlockSize, rBlockSize;
    numBlocks=target.indexTable.numBlocksLP(i);
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      rBlockSize=target.indexTable.rBlockSizeLP(i,iBlock);
      lBlockSize=target.indexTable.lBlockSizeLP(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  if(QNsVec[0].primaryIndex(i,target.indexTable.aiBlockIndexLP(i,iBlock,j)) && QNsVec[0].primaryIndex(i-1,target.indexTable.aimBlockIndexLP(i,iBlock,k))){
	    target.global_access(i,target.indexTable.siBlockIndexLP(i,iBlock,k),target.indexTable.aiBlockIndexLP(i,iBlock,j),target.indexTable.aimBlockIndexLP(i,iBlock,k))=exactGroundStateEntry(i,target.indexTable.siBlockIndexLP(i,iBlock,k),target.indexTable.aiBlockIndexLP(i,iBlock,j),target.indexTable.aimBlockIndexLP(i,iBlock,k));
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

lapack_complex_double exactGroundState::exactGroundStateEntry(int i, int si, int ai, int aim){
  if(si==0 || si==2){
    return 1.0;
  }
  if(imag(QNsVec[0].QNLabel(i-1,aim))==-1){
    return -1.0;
  }
  return 1.0;
}
