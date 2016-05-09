#include "imps.h"

imps::imps():mps(),impBase()
{}

//---------------------------------------------------------------------------------------------------//

imps::imps(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin):
  impBase(dimInfoIn),
  mps(dimInfoIn,conservedQNsin)
{
  centralIndexTableVar=twositeQNOrderMatrix((dimInfo.L()-1)/2,dimInfo,conservedQNs);
}

imps::imps(mps const &source):
  impBase(source.getDimInfo()),
  mps(source)
{
  centralIndexTableVar=twositeQNOrderMatrix((dimInfo.L()-1)/2,dimInfo,conservedQNs);
}
  

//---------------------------------------------------------------------------------------------------//

int imps::addSite(int Lnew, int i){
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].grow(Lnew,i,internalTargetQNBuffer[iQN],internalSourceBuffer);
  }
  stateArray::setParameterL(Lnew);
  loadIndexTables();

  /*
  std::cout<<std::endl;
  int const D=dimInfo.D();
  for(int i=-1;i<Lnew;++i){
    for(int ai=0;ai<D;++ai){
      std::cout<<conservedQNs[0].QNLabel(i,ai)<<"\t";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  */

  centralIndexTableVar=twositeQNOrderMatrix(i,dimInfo,conservedQNs);
  return 0;
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
// Beware, refineQN only stores the rightSideLabels right now, the global adaption can only be done 
// when the system is grown. 
//---------------------------------------------------------------------------------------------------//

int imps::refineQN(int i, std::vector<std::complex<int> > const &leftSideLabels, std::vector<std::complex<int> > const &rightSideLabels, std::vector<std::complex<int> > const &targetQN){
  internalTargetQNBuffer=targetQN;
  internalSourceBuffer=rightSideLabels;
  return conservedQNs[0].refine(i,leftSideLabels);
}
