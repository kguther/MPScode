#include <time.h>
#include <iostream>
#include <math.h>
#include "Qsystem.h"
#include "network.h"
#include "parameters.h"

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
    }
    TensorNetwork.setSimParameters(simPars);
    converged=TensorNetwork.solve(E0,dE);
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

int Qsystem::measure(mpo<lapack_complex_double> *const MPOperator, double &expectationValue, int iEigen){
  return TensorNetwork.measure(MPOperator,expectationValue,iEigen);
}

//---------------------------------------------------------------------------------------------------//

int Qsystem::measureLocalOperators(localMpo<lapack_complex_double> *const localMPOperator, std::vector<std::complex<double> > &result, int iEigen){
  return TensorNetwork.measureLocalOperators(localMPOperator,result,iEigen);
}
