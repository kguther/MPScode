#include <iostream>
#include <math.h>
#include "Qsystem.h"
#include "network.h"
#include "parameters.h"

//---------------------------------------------------------------------------------------------------//
/*
void testIDMRG(problemParameters pars, simulationParameters simPars, mpo<mpsEntryType > const &H){
  infiniteQsystem test(pars,simPars,H);
  test.solve();
  exit(1);
}
*/
//---------------------------------------------------------------------------------------------------//

Qsystem::Qsystem(problemParameters &inputpars, simulationParameters &inputsimPars):
  pars(inputpars),
  simPars(inputsimPars),
  DMax(inputsimPars.D),
  nSweepsMax(inputsimPars.nSweeps),
  tolInitialMax(inputsimPars.tolInitial),
  TensorNetwork(network(inputpars,inputsimPars))
{
  simPars.nSweeps=stageNSweeps(0);
  simPars.tolInitial=stageTolInitial(0);
}

//---------------------------------------------------------------------------------------------------//

void Qsystem::getGroundState(){
  int converged;
  for(int iStage=0;iStage<simPars.nStages;++iStage){
    //Start with low D to find a initial state, then increase D over the course of simulation
    if(iStage!=0){
      simPars.D=stageD(iStage);
      simPars.nSweeps=stageNSweeps(iStage);
      simPars.tolInitial=stageTolInitial(iStage);
    }
    TensorNetwork.setSimParameters(simPars);
    converged=TensorNetwork.solve(E0,dE);
    //returns 0 if dE is smaller than some threshold parameter
    if(converged==0){
      std::cout<<"SIMULATION CONVERGED\n";
      break;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::stageD(int const nStage){
  return (nStage+1)*DMax;
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::stageNSweeps(int const nStage){
  //Also, use only a few sweeps in the warmup phase
  if(nStage==(simPars.nStages-1)){
    return nSweepsMax;
  }
  return 2;
}

double Qsystem::stageTolInitial(int const nStage){
  if(nStage==(simPars.nStages-1)){
    return tolInitialMax;
  }
  double warmupTol=1e-2;
  return warmupTol*pow(tolInitialMax/warmupTol,static_cast<float>(nStage)/(simPars.nStages-1));
}

//---------------------------------------------------------------------------------------------------//

void Qsystem::measure(mpo<mpsEntryType> *const MPOperator, double &expectationValue, int iEigen){
  TensorNetwork.measure(MPOperator,expectationValue,iEigen);
}

//---------------------------------------------------------------------------------------------------//

void Qsystem::measureLocalOperators(localMpo<mpsEntryType> *const localMPOperator, std::vector<mpsEntryType > &result, int iEigen){
  TensorNetwork.measureLocalOperators(localMPOperator,result,iEigen);
}
