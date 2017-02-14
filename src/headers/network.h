#ifndef NETWORK
#define NETWORK

#include <vector>
#include "parameters.h"
#include "templates/mpo.h"
#include "templates/localMpo.h"
#include "mps.h"
#include "iterativeMeasurement.h"
#include "projector.h"
#include "quantumNumber.h"
#include "dimensionTable.h"
#include "mpstype.h"

//---------------------------------------------------------------------------------------------------//
// The network class contains all information required for a run of the simulation, that is, the whole
// MPS, the Hamiltonian in MPO representation and the partial contractions.
//---------------------------------------------------------------------------------------------------//

class network{
 public:
  network();
  network(problemParameters const &inputpars, simulationParameters const &inputsimPars);

  //These functions are the main interface for getting results

  //solves for the first nEigen eigenstates of networkH. On exit, lambda are the eigenenergies and deltaLambda the variances of energy.
  //can throw a critical_error in case normalization fails completely - this should never be the case
  int solve(std::vector<double> &lambda, std::vector<double> &deltaLambda);

  //Gets the expectation value of some operator MPOperator (input as MPO) for the iEigen-th state if obtained. On exit, expValue contains the result. 
  void measure(mpo<mpsEntryType > *const MPOperator, double &expValue, int iEigen=0);

  //Does the same thing as measure() but takes a local operator and moves it across the system, measuring at each site
  void measureLocalOperators(localMpo<mpsEntryType > *const MPOperator, std::vector<mpsEntryType > &expValue, int iEigen=0);

  //Besides observables, also the entanglement spectrum and entropy can be obtained from the MPS
  void getEntanglement(std::vector<double> &S, std::vector<std::vector<double> > &spectrum, int iEigen=0);

  //MPO needs to be initialized externally
  void setNetworkH(mpo<mpsEntryType > const &newH){networkH=newH;}

//---------------------------------------------------------------------------------------------------//

//Functions below are more advanced and might not be required in all applications

  //load and export does precisely what is expected
  void loadNetworkState(mps const &source);
  void exportNetworkState(mps &target);

  //There are internal flags for marking convergence - one should reset them if the same network shall do another run
  void resetConvergence();

  //Also, resetting the state can (!) be a good idea if subsequent runs are executed
  void resetState();
  //Similarly, resetting the quantum number might sometimes be necesary
  void setQuantumNumber(std::vector<std::complex<int> > const &targetQNs, std::vector<std::vector<std::complex<int> > > const &localQNs);

  //More obscure stuff for special applications
  void quantumNumberVec(std::vector<quantumNumber> *target){target=&conservedQNs;}
  dimensionTable& dimTable() {return networkDimInfo;}
  int setSimParameters(simulationParameters const &newPars);
  mpo<mpsEntryType > const& getNetworkH() const {return networkH;}
  int locd(int const i);

  //This is only for consistency checks
  mpo<mpsEntryType > *check, *checkParity;

//---------------------------------------------------------------------------------------------------//

 private:
  //Order dependent, do not change
  problemParameters pars;
  simulationParameters simPars;
  int D,L,Dw;
  projector excitedStateP;
  dimensionTable networkDimInfo;
  int lDL, lDR, ld, lDwR, lDwL;
  mps networkState;
  mpo<mpsEntryType > networkH;
  std::vector<int> nConverged;
  double shift, alpha;
  std::vector<quantumNumber> conservedQNs;
  iterativeMeasurement pCtr;
  mpsEntryType expectationValue;
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
  void getPExpressionLeft(int i, mpsEntryType *pExpr);
  void getPExpressionRight(int i, mpsEntryType *pExpr);
  void getLocalDimensions(int i);
  //For exception handling
  void resetSweep();
  int reSweep;
  //This one is only for consistency checks
  int checkQN();
  int checkEqualWeightState();
  void checkContractions(int i);
  mpsEntryType *backupCtr;
};

#endif
