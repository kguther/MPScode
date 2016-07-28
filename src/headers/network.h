#ifndef NETWORK
#define NETWORK

#include <complex>
#include <arcomp.h>
#include <vector>
#include "parameters.h"
#include "mpo.h"
#include "localMpo.h"
#include "mps.h"
#include "iterativeMeasurement.h"
#include "projector.h"
#include "quantumNumber.h"
#include "dimensionTable.h"

//---------------------------------------------------------------------------------------------------//
// The network class contains all information required for a run of the simulation, that is, the whole
// MPS, the Hamiltonian in MPO representation and the partial contractions.
//---------------------------------------------------------------------------------------------------//

class network{
 public:
  network();
  //can throw a critical_error in case normalization fails completely
  network(problemParameters const &inputpars, simulationParameters const &inputsimPars);

  //These functions are the main interface for getting results

  //solves for the first nEigen eigenstates of networkH. On exit, lambda are the eigenenergies and deltaLambda the variances of energy.
  int solve(std::vector<double> &lambda, std::vector<double> &deltaLambda);
  //Gets the expectation value of some operator MPOperator (input as MPO) for the iEigen-th state if obtained. On exit, expValue contains the result. 
  void measure(mpo<arcomplex<double> > *const MPOperator, double &expValue, int iEigen=0);
  //Does the same thing as measure() but takes a local operator and moves it across the system, measuring at each site
  void measureLocalOperators(localMpo<arcomplex<double> > *const MPOperator, std::vector<arcomplex<double> > &expValue, int iEigen=0);
  //Besides observables, also the entanglement spectrum and entropy can be obtained from the MPS
  void getEntanglement(std::vector<double> &S, std::vector<std::vector<double> > &spectrum, int iEigen=0);
  //MPO needs to be initialized externally
  void setNetworkH(mpo<arcomplex<double> > const &newH){networkH=newH;}

//---------------------------------------------------------------------------------------------------//

//Functions below are more advanced and might not be required in all applications

  void loadNetworkState(mps const &source);
  void exportNetworkState(mps &target);
  void resetConvergence();
  void resetState();
  void quantumNumberVec(std::vector<quantumNumber> *target){target=&conservedQNs;}
  dimensionTable& dimTable() {return networkDimInfo;}
  int setSimParameters(simulationParameters const &newPars);
  mpo<arcomplex<double> > const& getNetworkH() const {return networkH;}
  int locd(int const i);
  //This is only for consistency checks
  mpo<arcomplex<double> > *check, *checkParity;
 private:
  //Order dependent, do not change
  problemParameters pars;
  simulationParameters simPars;
  int D,L,Dw;
  projector excitedStateP;
  dimensionTable networkDimInfo;
  int lDL, lDR, ld, lDwR, lDwL;
  mps networkState;
  mpo<arcomplex<double> > networkH;
  std::vector<int> nConverged;
  double shift, alpha;
  std::vector<quantumNumber> conservedQNs;
  iterativeMeasurement pCtr;
  arcomplex<double> expectationValue;
  //most of these methods are auxiliary functions
  int pctrIndex(int ai, int bi, int aip){return aip+bi*D+ai*D*Dw;}
  int stateIndex(int si, int ai, int aim){return aim+lDL*ai+lDL*lDR*si;}
  int optimize(int i, int maxIter, double tol, double &iolambda);
  int locDMax(int i);
  int gotoNextEigen();
  int setParameterD(int Dnew);
  double convergenceCheck();
  double getCurrentEnergy(int i);
  void getNewAlpha(int i, double &lambda, double prevLambda);
  void normalize(int i, int direction, int enrichment=0);
  void refineQNLabels(int i, std::vector<std::complex<int> > const &source);
  void sweep(double maxIter, double tol, double &lambda);
  void leftEnrichment(int i);
  void rightEnrichment(int i);
  void leftEnrichmentBlockwise(int i);
  void rightEnrichmentBlockwise(int i);
  void calcHSqrExpectationValue(double &ioHsqr);
  void getPExpressionLeft(int i, arcomplex<double> *pExpr);
  void getPExpressionRight(int i, arcomplex<double> *pExpr);
  void getLocalDimensions(int i);
  //For exception handling
  void resetSweep();
  int reSweep;
  //This one is only for consistency checks
  int checkQN();
  int checkEqualWeightState();
  void checkContractions(int i);
  arcomplex<double> *backupCtr;
};

#endif
