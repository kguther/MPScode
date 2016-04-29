#ifndef IDMRG_INITIAL_STATE_GROWER
#define IDMRG_INITIAL_STATE_GROWER

#include "imps.h"
#include "parameters.h"
#include "mpo.h"
#include <complex>

class initStateGrow{
 public:
  initStateGrow(problemParameters const &parsIn, simulationParameters const &simParsIn, mpo<std::complex<double> > const &H);
  void prepareInitialState(mps &target);
  void exportState(mps &target);
 private:
  int const initialLength;
  problemParameters pars;
  simulationParameters simPars;
  imps networkState;
  mpo<std::complex<double> > networkH;
  void growSystem();
};

#endif
