#include "infiniteNetwork.h"
#include "twositeHMatrix.h"
#include "uncachedMeasurement.h"
#include <arcomp.h>
#include <arscomp.h>
#include <memory>

infiniteNetwork::infiniteNetwork(problemParameters parsIn, simulationParameters simParsIn):
  pars(parsIn),
  simPars(simParsIn),
  dimInfo(dimensionTable(simParsIn.D,2,parsIn.d))
{
  networkH=mpo<lapack_complex_double>(pars.d.maxd(),pars.Dw,pars.L);
  for(int iQN=0;iQN<parsIn.QNconserved.size();++iQN){
    conservedQNs.push_back(quantumNumber(dimInfo,pars.QNconserved[iQN],pars.QNLocalList[iQN]));
  }
  networkState=imps(dimInfo,conservedQNs);
  pCtr=uncachedMeasurement(&networkH,&networkState);
  diags.push_back(1.0);
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::iDMRGStep(){
  arcomplex<double> *bufferMatrix;
  i=(dimInfo.L()-1)/2;
  std::unique_ptr<arcomplex<double> > bufferMatrixP(new arcomplex<double>[dimInfo.locd(i)*dimInfo.locd(i+1)*dimInfo.locDimL(i)*dimInfo.locDimR(i+1)]);
  bufferMatrix=bufferMatrixP.get();
  
  //State prediction needs at least two optimized matrices, not available in the first step
  if(i!=0){
    statePrediction(bufferMatrix);
  }
  else{
    for(int m=0;m<dimInfo.locd(i)*dimInfo.locd(i+1);++m){
      bufferMatrix[m]=1;
    }
  }

  optimize(bufferMatrix);
  updateMPS(bufferMatrix);
  addSite();
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::addSite(){
  pCtr.update();
  std::vector<std::complex<int> > newQNs;
  for(int iQN=0;iQN<pars.QNconserved.size();++iQN){
    newQNs[iQN].real(static_cast<int>(pars.filling[iQN]*pars.L));
    newQNs[iQN].imag(pars.QNconserved[iQN].imag());
  }
  networkState.addSite(newQNs);
  dimInfo.setParameterL(L+2);
}
