#ifndef MPS_CLASS_QSYSTEM
#define MPS_CLASS_QSYSTEM

#include <complex>
#include <vector>
#include "parameters.h"
#include "network.h"
#include "mpo.h"
#include "localMpo.h"

//---------------------------------------------------------------------------------------------------//
// This class is used to execute the simulation with its main task being the adaption of bond
// dimension D over the course of the run. 
// It is basically an extension of the network class
//---------------------------------------------------------------------------------------------------//

class Qsystem{
 public:
  Qsystem(problemParameters &pars, simulationParameters &simPars);
  network TensorNetwork;
  std::vector<double> E0, dE;
  void getGroundState();
  void measure(mpo<std::complex<double> > *const MPOperator, double &expectationValue, int iEigen=0);
  void measureLocalOperators(localMpo<std::complex<double> > *const localMPOperator, std::vector<std::complex<double> > &result, int iEigen=0);
 protected:
  problemParameters pars;
  simulationParameters simPars;
 private:
  int DMax, nSweepsMax;
  double tolInitialMax;
  int stageD(int const nStage);
  int stageNSweeps(int const nStage);
  double stageTolInitial(int const nStage);
};

#endif
