#ifndef MPS_CLASS_QSYSTEM
#define MPS_CLASS_QSYSTEM

#include <complex>
#include "parameters.h"
#include "network.h"
#include "mpo.h"

//---------------------------------------------------------------------------------------------------//
// This class is used to execute the simulation with its main task being the adaption of bond
// dimension D over the course of the run. Also, this is where measurements are going to be implemented.
//---------------------------------------------------------------------------------------------------//

class Qsystem{
 public:
  Qsystem(problemParameters pars, simulationParameters simPars);
  ~Qsystem();
  network TensorNetwork;
  double *E0;
  int getGroundState();
  int measure(mpo<std::complex<double> > &MPOperator, double &expectationValue);
 private:
  problemParameters pars;
  simulationParameters simPars;
  int DMax, nSweepsMax;
  double tolInitialMax;
  int stageD(int const nStage);
  int stageNSweeps(int const nStage);
  double stageTolInitial(int const nStage);
};

#endif
