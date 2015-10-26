#include "Qsystem.h"
#include "network.h"
#include "parameters.h"

Qsystem::Qsystem(parameters inputpars){
  pars=inputpars;
  parameters initialpars(pars.d,stageD(0),pars.L,pars.Dw,stageNSweeps(0),pars.nEigs,pars.nStages);
  TensorNetwork.generate(initialpars);
}

int Qsystem::getGroundState(){
  for(int iStage=0;iStage<pars.nStages;iStage++){
    TensorNetwork.setParameterNSweeps(stageNSweeps(iStage));
    TensorNetwork.setParameterD(pars.D);//stageD(iStage));
    E0=TensorNetwork.solve();
  }
  return 0;
}

int Qsystem::stageD(int nStage){
  int D0=(pars.D>50)?pars.D/10:5;
  return D0+nStage*(pars.D-D0)/(pars.nStages-1);
}

int Qsystem::stageNSweeps(int nStage){
  int nSweeps0=2;
  if(pars.nSweeps>2){
    return nSweeps0+nStage*(pars.nSweeps-nSweeps0)/(pars.nStages-1);
  }
  else{
    return 2;
  }
}
