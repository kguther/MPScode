#include <time.h>
#include <iostream>
#include "Qsystem.h"
#include "network.h"
#include "parameters.h"

Qsystem::Qsystem(parameters inputpars){
  pars=inputpars;
  parameters initialpars(pars.d,stageD(0),pars.L,pars.Dw,stageNSweeps(0),pars.nEigs,pars.nStages,pars.acc);
  TensorNetwork.initialize(initialpars);
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::getGroundState(){
  clock_t simTime;
  double alphaInitial=1e-2;
  int converged;
  simTime=clock();
  for(int iStage=0;iStage<pars.nStages;++iStage){
    //Start with low D to find a initial state, then increase D over the course of simulation
    TensorNetwork.setParameterNSweeps(stageNSweeps(iStage));
    TensorNetwork.setParameterD(stageD(iStage));
    TensorNetwork.setParameterAlpha(alphaInitial);
    converged=TensorNetwork.solve(E0);
    if(converged==0){
      std::cout<<"SIMULATION CONVERGED\n";
      break;
    }
  }
  simTime=clock()-simTime;
  std::cout<<"Total simulation took "<<(float)simTime/CLOCKS_PER_SEC<<" seconds\n";
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::stageD(int nStage){
  if(pars.nStages>1){
    int D0=(pars.D>10)?pars.D/4:5;
    return D0+nStage*(pars.D-D0)/(pars.nStages-1);
  }
  else{
    return pars.D;
  }
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::stageNSweeps(int nStage){
  //Also, use only a few sweeps in the warmup phase
  if(nStage==(pars.nStages-1)){
    return pars.nSweeps;
  }
  else{
    return 2;
  }
}
