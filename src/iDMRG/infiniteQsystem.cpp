#include "infiniteQsystem.h"
#include "infiniteNetwork.h"
#include "quantumNumber.h"

infiniteQsystem::infiniteQsystem(problemParameters parsIn, simulationParameters simParsIn, mpo<std::complex<double> > const &HIn):
  pars(parsIn),
  simPars(simParsIn),
  unitCellSize(2),
  networkH(HIn)
{
  //Generate initial state with the desired quantum numbers
  dimensionTable dimInfo(simParsIn.D,unitCellSize,pars.d);
  std::complex<int> qnValue;
  std::vector<quantumNumber> conservedQNs;
  for(int iQN=0;iQN<parsIn.QNconserved.size();++iQN){
    qnValue=std::complex<int>(pars.filling[iQN]*unitCellSize,pars.QNconserved[iQN].imag());
    conservedQNs.push_back(quantumNumber(dimInfo,qnValue,pars.QNLocalList[iQN]));
  }
  networkState=timps(unitCellSize,dimInfo,conservedQNs);
}

void infiniteQsystem::solve(){
  infiniteNetwork sys(pars,simPars,&networkState);
  int const maxSteps=10;
  for(int nStep=0;nStep<maxSteps;++nStep){
    sys.iDMRGStep();
  }
}
