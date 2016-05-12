#ifndef NETWORK_CLASS_FOR_IDMRG
#define NETWORK_CLASS_FOR_IDMRG

#include "parameters.h"
#include "impBase.h"
#include "dimensionTable.h"
#include "uncachedMeasurement.h"
#include "baseTensor.h"
#include "truncation.h"
#include <arcomp.h>
#include <vector>

void verifyCompression(arcomplex<double> *cVector, int dim);

class infiniteNetwork{
 public:
  infiniteNetwork(problemParameters const &parsIn, simulationParameters const &simParsIn, impBase *MPState);
  void growTLSystem();
  void iDMRGStep();
  void addDiags();
  int statePrediction(arcomplex<double> *target);
  void qnEnforcedPrediction(arcomplex<double> *target);
  int optimize(arcomplex<double> *target);
  void updateMPS(arcomplex<double> *source);
  void setPCtr(std::vector<arcomplex<double> > const &R, std::vector<arcomplex<double> > const &L);
  impBase* getState();
  mpo<arcomplex<double> > networkH;
 private:
  int i;
  int firstStep;
  problemParameters pars;
  simulationParameters simPars;
  dimensionTable dimInfo;
  impBase *networkState;
  std::vector<std::complex<int> > optLocalQNsL, optLocalQNsR;
  std::vector<double> diags, diagsm;
  //Beware that iDMRG builds up a regular system -> only three MPO matrices are referred - in particular is a MPO length of at least 3 required
  uncachedMeasurement pCtr;
  baseTensor<arcomplex<double> > aBuf, bBuf;
  int explicitIndex(int iBlock, int j, int k);
  void addSite();
};

#endif
