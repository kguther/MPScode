#include <time.h>
#include <iostream>
#include <math.h>
#include "Qsystem.h"
#include "network.h"
#include "parameters.h"

Qsystem::Qsystem(problemParameters inputpars, simulationParameters inputsimPars){
  pars=inputpars;
  simPars=inputsimPars;
  DMax=simPars.D;
  nSweepsMax=simPars.nSweeps;
  tolInitialMax=simPars.tolInitial;
  simPars.D=stageD(0);
  simPars.nSweeps=stageNSweeps(0);
  simPars.tolInitial=stageTolInitial(0);
  TensorNetwork.initialize(pars,simPars);
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::getGroundState(){
  clock_t simTime;
  int converged;
  simTime=clock();
  for(int iStage=0;iStage<simPars.nStages;++iStage){
    //Start with low D to find a initial state, then increase D over the course of simulation
    if(iStage!=0){
      simPars.D=stageD(iStage);
      simPars.nSweeps=stageNSweeps(iStage);
      simPars.tolInitial=stageTolInitial(iStage);
      TensorNetwork.setSimParameters(simPars);
    }
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

int Qsystem::stageD(int const nStage){
  if(simPars.nStages>1){
    int D0=(DMax>10)?DMax/4:5;
    return D0+nStage*(DMax-D0)/(simPars.nStages-1);
  }
  return DMax;
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::stageNSweeps(int const nStage){
  //Also, use only a few sweeps in the warmup phase
  if(nStage==(simPars.nStages-1)){
    return nSweepsMax;
  }
  return 1;
}

double Qsystem::stageTolInitial(int const nStage){
  if(nStage==(simPars.nStages-1)){
    return tolInitialMax;
  }
  double warmupTol=1e-2;
  return warmupTol*pow(tolInitialMax/warmupTol,static_cast<float>(nStage)/(simPars.nStages-1));
}
