#include "imps.h"

imps::imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin):
  mps(dimInfo,conservedQNsin),
  centralIndexTable(twositeQNOrderMatrix(dimInfo.L()/2,dimInfo,&conservedQNsin))
{}

//---------------------------------------------------------------------------------------------------//

void imps::addSite(std::complex<int> *targetQN){
  stateArray::setParameterL(L+2);
  for(int iQN=0;iQN<nQNs;++iQN){
    //Not sure if this correctly updates quantum numbers since this relies on the layout being unchanged when QN is shifted
    conservedQNs[iQN]=quantumNumber(dimInfo,targetQN[iQN],conservedQNs[iQN].localQNValue());
  }
  setUpQNs(conservedQNs);
  centralIndexTable=twositeQNOrderMatrix(L/2,dimInfo,&conservedQNs);
}

//---------------------------------------------------------------------------------------------------//

void imps::exportState(mps &target){
  target=mps(dimInfo,conservedQNs);
  int lDL, lDR, ld;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    ld=locd(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  target.global_access(i,si,ai,aim)=global_access(i,si,ai,aim);
	}
      }
    }
  }
}
