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
  QNsVec.push_back(quantumNumber(target.getDimInfo(),QNValue,localQNsVec));
  mps proxy(target.getDimInfo(),QNsVec);
  generateExactState(proxy);
  target.setStateArray(proxy);
}

//---------------------------------------------------------------------------------------------------//

void exactGroundState::generateExactState(mps &target){
    //This is the exact ground state at the critical point for fixed particle number and subchain parity. It turns out that this is a nice guess for the ground state of the perturbed system (for small perturbations).
  int ld, lDL, lDR;
  siteQNOrderMatrix localIndexTable;
  for(int i=0;i<target.getDimInfo().L();++i){
    localIndexTable=target.indexTable().getLocalIndexTable(i);
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
    numBlocks=localIndexTable.numBlocksLP();
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  if(QNsVec[0].primaryIndex(i,localIndexTable.aiBlockIndexLP(iBlock,j)) && QNsVec[0].primaryIndex(i-1,localIndexTable.aimBlockIndexLP(iBlock,k))){
	    target.global_access(i,localIndexTable.siBlockIndexLP(iBlock,k),localIndexTable.aiBlockIndexLP(iBlock,j),localIndexTable.aimBlockIndexLP(iBlock,k))=exactGroundStateEntry(i,localIndexTable.siBlockIndexLP(iBlock,k),localIndexTable.aiBlockIndexLP(iBlock,j),localIndexTable.aimBlockIndexLP(iBlock,k));
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

mpsEntryType exactGroundState::exactGroundStateEntry(int i, int si, int ai, int aim){
  if(si==0 || si==2){
    return 1.0;
  }
  if(imag(QNsVec[0].QNLabel(i-1,aim))==-1){
    return std::complex<double>(0,-1.0);
  }
  return std::complex<double>(0,1.0);
}
