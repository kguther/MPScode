#ifndef MPS_CLASS_QSYSTEM
#define MPS_CLASS_QSYSTEM

#include <complex>
#include "parameters.h"
#include "network.h"

class Qsystem{
 public:
  Qsystem(parameters pars);
  network TensorNetwork;
  lapack_complex_double *****MPO;
  double E0;
  int getGroundState();
 private:
  parameters pars;
  int stageD(int nStage);
  int stageNSweeps(int nStage);
};

#endif
