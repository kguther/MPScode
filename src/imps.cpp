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

int imps::addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN, std::vector<std::complex<int> > const &source){
  int info;
  for(int iQN=0;iQN<nQNs;++iQN){
    conservedQNs[iQN].grow(Lnew,i,targetQN[iQN],source);
  }
  stateArray::setParameterL(Lnew);
  info=loadIndexTables();

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
  return info;
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
