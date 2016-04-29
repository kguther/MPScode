#include "initStateGrow.h"
#include "infiniteNetwork.h"
#include "network.h"

initStateGrow::initStateGrow(problemParameters const &parsIn, simulationParameters const &simParsIn, mpo<std::complex<double > > const &H):
  pars(parsIn),
  simPars(simParsIn),
  networkH(H),
  initialLength(4)
{
  //Generate initial state with the desired quantum numbers
  dimensionTable dimInfo(simParsIn.D,initialLength,pars.d);
  std::complex<int> qnValue;
  std::vector<quantumNumber> conservedQNs;
  for(int iQN=0;iQN<parsIn.QNconserved.size();++iQN){
    qnValue=std::complex<int>(pars.filling[iQN]*2*dimInfo.L(),pars.QNconserved[iQN].imag());
    conservedQNs.push_back(quantumNumber(dimInfo,qnValue,pars.QNLocalList[iQN]));
  }
  networkState=imps(dimInfo,conservedQNs);
}

//---------------------------------------------------------------------------------------------------//

void initStateGrow::growSystem(){
  //Starting from a system of size initialLength, one of the target size L is generated
  int const L=pars.L;
  infiniteNetwork growingNetwork(pars,simPars,&networkState);
  growingNetwork.networkH=networkH;

  for(int m=0;m<(L-initialLength)/2+1;++m){
    //Since the iDMRG step starts with optimizing the two latest (i.e. middle) matrices, these are now optimized again, but this is just a small overhead. 
    //Inserting new matrices before would not work with the QN update scheme as this relies on the results of the SVD
    growingNetwork.iDMRGStep();
  }
  growingNetwork.addDiags();
}

//---------------------------------------------------------------------------------------------------//

void initStateGrow::prepareInitialState(mps &target){
  problemParameters parsPre=pars;
  parsPre.L=initialLength;
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    parsPre.QNconserved[iQN].real(pars.filling[iQN]*initialLength);
  }
  //The reduced system of size initialLength has to be solved (using the finite-size DMRG)
  network preRun(parsPre,simPars);
  //Copy a piece of networkH of size initialLength into the preparing network. Minimal size of H is still 3 sites.
  preRun.networkH.loadTI(networkH);
  std::vector<double> energies, dEnergies;
  preRun.solve(energies,dEnergies);
  //The network class works with mps while the infiniteNetwork uses the imps child class. A dummy mps is used for conversion.
  mps dummy;
  preRun.exportNetworkState(dummy);
  networkState=imps(dummy);

  growSystem();
  //Write the resulting initial state into the target
  exportState(target);
}

//---------------------------------------------------------------------------------------------------//

void initStateGrow::exportState(mps &target){
  networkState.exportState(target);
}
