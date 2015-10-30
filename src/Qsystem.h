#ifndef MPS_CLASS_QSYSTEM
#define MPS_CLASS_QSYSTEM

#include <complex>
#include "parameters.h"
#include "network.h"

//---------------------------------------------------------------------------------------------------//
// This class is used to execute the simulation with its main task being the adaption of bond
// dimension D over the course of the run. Also, this is where measurements are going to be implemented.
//---------------------------------------------------------------------------------------------------//

class Qsystem{
 public:
  Qsystem(parameters pars);
  network TensorNetwork;
  double E0;
  int getGroundState();
 private:
  parameters pars;
  int stageD(int nStage);
  int stageNSweeps(int nStage);
};

#endif
