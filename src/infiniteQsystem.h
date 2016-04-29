#ifndef IDMRG_ALGORITHM_TARGETING_TI
#define IDMRG_ALGORITHM_TARGETING_TI

#include "parameters.h"
#include "timps.h"
#include "mpo.h"
#include <complex>

class infiniteQsystem{
 public:
  infiniteQsystem(problemParameters parsIn, simulationParameters simParsIn, mpo<std::complex<double> > const &HIn);
  void solve();
 private:
  problemParameters pars;
  simulationParameters simPars;
  timps networkState;
  int unitCellSize;
  mpo<std::complex<double> > networkH;
};

#endif
