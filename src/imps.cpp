#include "imps.h"

imps::imps():mps()
{}

//---------------------------------------------------------------------------------------------------//

imps::imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin):
  mps(dimInfo,conservedQNsin),
  centralIndexTable(twositeQNOrderMatrix((dimInfo.L()-1)/2,dimInfo,&conservedQNsin))
{}

//---------------------------------------------------------------------------------------------------//

void imps::addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN){
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].grow(Lnew,i,targetQN[iQN]);
  }
  stateArray::setParameterL(Lnew);
  setUpQNs(conservedQNs);
  centralIndexTable=twositeQNOrderMatrix(i,dimInfo,&conservedQNs);
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
  target.setUpQNs(conservedQNs);
}

//---------------------------------------------------------------------------------------------------//

int imps::refineQN(int i, std::vector<std::complex<int> > const &source){
  return conservedQNs[0].refine(i,source);
}
