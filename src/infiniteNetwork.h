#ifndef NETWORK_CLASS_FOR_IDMRG
#define NETWORK_CLASS_FOR_IDMRG

#include "parameters.h"
#include "imps.h"
#include "quantumNumber.h"
#include "dimensionTable.h"
#include "uncachedMeasurement.h"
#include <arcomp.h>

struct sortData{
  double lambda;
  std::complex<int> QN;
  int index;
};

bool compareSortData(sortData const &a, sortData const &b);

class infiniteNetwork{
 public:
  infiniteNetwork(problemParameters const &parsIn, simulationParameters const &simParsIn, imps *MPState);
  void growTLSystem();
  void iDMRGStep();
  void addDiags();
  void statePrediction(arcomplex<double> *target);
  int optimize(arcomplex<double> *target);
  void updateMPS(arcomplex<double> *source);
  imps* getState();
  mpo<lapack_complex_double> networkH;
 private:
  int i;
  problemParameters pars;
  simulationParameters simPars;
  dimensionTable dimInfo;
  imps *networkState;
  std::vector<std::complex<int> > optLocalQNs;
  std::vector<double> diags, diagsm;
  //Beware that iDMRG builds up a regular system -> only three MPO matrices are referred - in particular is a MPO length of at least 3 required
  uncachedMeasurement pCtr;
  int explicitIndex(int iBlock, int j, int k);
  void addSite();
};

#endif
