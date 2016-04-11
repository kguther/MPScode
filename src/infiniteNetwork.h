#ifndef NETWORK_CLASS_FOR_IDMRG
#define NETWORK_CLASS_FOR_IDMRG

#include "parameters.h"
#include "imps.h"
#include "quantumNumber.h"
#include "dimensionTable.h"
#include "uncachedMeasurement.h"
#include <arcomp.h>

class infiniteNetwork{
 public:
  infiniteNetwork(problemParameters parsIn, simulationParameters simParsIn);
  void growSystem(int L);
  void growTLSystem();
  void iDMRGStep();
  void addSite();
  void statePrediction(arcomplex<double> *target);
  int optimize(arcomplex<double> *target);
  void updateMPS(arcomplex<double> *source);
  void exportState(mps &target);
  mpo<lapack_complex_double> networkH;
 private:
  problemParameters pars;
  simulationParameters simPars;
  dimensionTable dimInfo;
  imps networkState;
  int i;
  std::vector<quantumNumber> conservedQNs;
  std::vector<double> diags, diagsm;
  //Beware that iDMRG builds up a regular system -> only three MPO matrices are referred - in particular is a MPO length of at least 3 required
  uncachedMeasurement pCtr;
  int explicitIndex(int iBlock, int j, int k);
};

#endif
