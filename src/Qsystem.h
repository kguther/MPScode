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
// dimension D over the course of the run. Also, this is where measurements are going to be implemented.
//---------------------------------------------------------------------------------------------------//

class Qsystem{
 public:
  Qsystem(problemParameters &pars, simulationParameters &simPars);
  ~Qsystem();
  network TensorNetwork;
  std::vector<double> E0, dE;
  int getGroundState();
  int measure(mpo<std::complex<double> > *const MPOperator, double &expectationValue, mps *const MPState=0);
  int measureLocal(localMpo<std::complex<double> > *const localMPOperator, std::vector<std::complex<double> > &result, mps *const MPState=0);
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
