#ifndef NETWORK
#define NETWORK

#include <complex>
#include <arcomp.h>
#include <vector>
#include <string>
#include "mkl_complex_defined.h"
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
  network(problemParameters const &inputpars, simulationParameters const &inputsimPars);
  int solve(std::vector<double> &lambda, std::vector<double> &deltaLambda);
  int measure(mpo<lapack_complex_double> *const MPOperator, double &expValue, int iEigen=0);
  int measureLocalOperators(localMpo<lapack_complex_double> *const MPOperator, std::vector<lapack_complex_double> &expValue, int iEigen=0);
  void getEntanglement(std::vector<double> &S, std::vector<std::vector<double> > &spectrum, int iEigen=0);
  void loadNetworkState(mps const &source);
  void resetConvergence();
  void getInitState();
  void quantumNumberVec(std::vector<quantumNumber> *target){target=&conservedQNs;}
  dimensionTable& dimTable() {return networkDimInfo;}
  int setSimParameters(simulationParameters const &newPars);
  //MPO needs to be initialized externally
  mpo<lapack_complex_double> networkH;
  int locd(int const i);
  //This is only for consistency checks
  void leftNormalizationMatrixFull();
  mpo<lapack_complex_double> *check, *checkParity;
 private:
  //Order dependent, do not change
  problemParameters pars;
  simulationParameters simPars;
  int D,L,Dw;
  projector excitedStateP;
  dimensionTable networkDimInfo;
  int lDL, lDR, ld, lDwR, lDwL;
  mps networkState;
  std::vector<int> nConverged;
  double shift, alpha;
  std::vector<quantumNumber> conservedQNs;
  iterativeMeasurement pCtr;
  lapack_complex_double expectationValue;
  //most of these methods are auxiliary functions
  int pctrIndex(int ai, int bi, int aip){return aip+bi*D+ai*D*Dw;}
  int stateIndex(int si, int ai, int aim){return aim+lDL*ai+lDL*lDR*si;}
  int optimize(int i, int maxIter, double tol, double &iolambda);
  int solveSiteEigenProb(int i, int maxIter, double tol, double &iolambda);
  int solveSiteEigenProbQNC(int i, int maxIter, double tol, double &iolambda);
  int locDMax(int i);
  int gotoNextEigen();
  int setParameterD(int Dnew);
  double convergenceCheck();
  double getCurrentEnergy(int i);
  void getNewAlpha(int i, double lambda, double prevLambda);
  void normalize(int i, int direction, int enrichment=0);
  void sweep(double maxIter, double tol, double &lambda);
  void leftEnrichment(int i);
  void rightEnrichment(int i);
  void leftEnrichmentBlockwise(int i);
  void rightEnrichmentBlockwise(int i);
  void calcHSqrExpectationValue(double &ioHsqr);
  void getPExpressionLeft(int i, lapack_complex_double *pExpr);
  void getPExpressionRight(int i, lapack_complex_double *pExpr);
  void getLocalDimensions(int i);
  //This one is only for consistency checks
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
  int checkQN();
  int checkEqualWeightState();
  void checkContractions(int i);
  lapack_complex_double *backupCtr;
};

#endif
