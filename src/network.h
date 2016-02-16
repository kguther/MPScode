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
  ~network();
  int solve(std::vector<double> &lambda, std::vector<double> &deltaLambda);
  int measure(mpo<lapack_complex_double> *const MPOperator, double &expValue, int iEigen=0);
  int measureLocalOperators(localMpo<lapack_complex_double> *const MPOperator, std::vector<lapack_complex_double> &expValue, int iEigen=0);
  void getEntanglement(std::vector<double> &S, std::vector<std::vector<double> > &spectrum, int iEigen=0);
  void initialize(problemParameters const &inputpars, simulationParameters const &inputSimPars);
  void loadNetworkState(mps const &source);
  void resetConvergence();
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
  network(network const &cpynet);//Copying networks is better avoided to save memory
  network& operator=(network const &cpynet);//Use the generate function instead, assignment is dangerous for networks with different parameters 
  mps networkState;
  projector excitedStateP;
  problemParameters pars;
  simulationParameters simPars;
  dimensionTable networkDimInfo;
  int D,L,Dw,icrit;
  int lDL, lDR, ld, lDwR, lDwL;
  int *nConverged;
  double shift, alpha;
  std::vector<quantumNumber> conservedQNs;
  iterativeMeasurement pCtr;
  lapack_complex_double expectationValue;
  //most of these methods are auxiliary functions
  int pctrIndex(int const ai, int const bi, int const aip){return aip+bi*D+ai*D*Dw;}
  int optimize(int const i, int const maxIter, double const tol, double &iolambda);
  int solveSiteEigenProb(int const i, int const maxIter, double const tol, double &iolambda);
  int solveSiteEigenProbQNC(int const i, int const maxIter, double const tol, double &iolambda);
  int locDMax(int const i);
  int gotoNextEigen();
  int setParameterD(int Dnew);
  double convergenceCheck();
  double getCurrentEnergy(int const i);
  void getNewAlpha(int const i, double const lambda, double const prevLambda);
  void normalize(int const i, int const direction, int const enrichment=0);
  void sweep(double const maxIter, double const tol, double &lambda);
  void leftEnrichment(int const i);
  void rightEnrichment(int const i);
  void leftEnrichmentBlockwise(int const i);
  void rightEnrichmentBlockwise(int const i);
  void calcHSqrExpectationValue(double &ioHsqr);
  void getPExpressionLeft(int const i, lapack_complex_double *pExpr);
  void getPExpressionRight(int const i, lapack_complex_double *pExpr);
  void getLocalDimensions(int const i);
  //This one is only for consistency checks
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
  int checkQN();
  int checkEqualWeightState();
  void checkContractions(int const i);
  lapack_complex_double *backupCtr;
};

#endif
