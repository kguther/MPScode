#ifndef NETWORK
#define NETWORK

#include <complex>
#include <arcomp.h>
#include "mkl_complex_defined.h"
#include "parameters.h"
#include "mpo.h"
#include "pContraction.h"
#include "mps.h"
#include "iterativeMeasurement.h"
#include "projector.h"
#include "overlap.h"

//---------------------------------------------------------------------------------------------------//
// The network class contains all information required for a run of the simulation, that is, the whole
// MPS, the Hamiltonian in MPO representation and the partial contractions.
//---------------------------------------------------------------------------------------------------//

class network{
 public:
  network();
  network(problemParameters inputpars, simulationParameters inputsimPars);
  ~network();
  int solve(double *lambda);
  int measure(mpo<lapack_complex_double> *MPOperator, double &expValue);
  void initialize(problemParameters inputpars, simulationParameters inputSimPars);
  void loadNetworkState(mps &source);
  int setSimParameters(simulationParameters newPars);
  //MPO needs to be initialized externally
  mpo<lapack_complex_double> networkH;
  //This one is only for consistency checks
  void leftNormalizationMatrixFull();
  mpo<lapack_complex_double> *check;
 private:
  network(network const &cpynet);//Copying and assigning networks is better avoided because it would work in a quite unintuitive way (the content of the array structures had to be copied, but the member itself must not be copied) and would be computationally quite expensive - but might be useful if one wanted to genuinely increase D during run (also: add a delete function for manual deletion)
  network& operator=(network const &cpynet);//Use the generate function instead, assignment is dangerous for networks with different parameters 
  mps networkState;
  projector excitedStateP;
  problemParameters pars;
  simulationParameters simPars;
  int d,D,L,Dw,icrit;
  int lDL, lDR, ld, lDwR, lDwL;
  int *nConverged;
  double shift;
  iterativeMeasurement pCtr;
  lapack_complex_double expectationValue;
  //most of these methods are auxiliary functions
  int pctrIndex(int const ai, int const bi, int const aip){return aip+bi*D+ai*D*Dw;}
  int optimize(int const i, int const maxIter, double const tol, double &iolambda);
  int locd(int const i);
  int locDMax(int const i);
  int gotoNextEigen();
  int setParameterD(int Dnew);
  double convergenceCheck();
  void sweep(double const maxIter, double const tol, double const alpha ,double &lambda);
  void leftEnrichment(double const alpha, int const i);
  void rightEnrichment(double const alpha, int const i);
  void calcHSqrExpectationValue(double &ioHsqr);
  void getPExpressionLeft(int const i, lapack_complex_double *pExpr);
  void getPExpressionRight(int const i, lapack_complex_double *pExpr);
  void getLocalDimensions(int const i);
  //This one is only for consistency checks
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
};

#endif
