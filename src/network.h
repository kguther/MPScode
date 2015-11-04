#ifndef NETWORK
#define NETWORK

#include <complex>
#include <vector>
#include <lapacke.h>
#include <arcomp.h>
#include "parameters.h"
#include "mpo.h"
#include "pContraction.h"
#include "mps.h"

//---------------------------------------------------------------------------------------------------//
// The network class contains all information required for a run of the simulation, that is, the whole
// MPS, the Hamiltonian in MPO representation and the partial contractions. It does not, however, 
// include measurements.
//---------------------------------------------------------------------------------------------------//

class network{
 public:
  network();
  network(parameters inputpars);
  double solve();
  void initialize(parameters inputpars);
  void setParameterNSweeps(int Nnew);
  int setParameterD(int Dnew);
  int calcCtrFull(int const direction);
  //MPO needs to be initialized externally
  mpo<lapack_complex_double> networkH;
  //This one is only for consistency checks
  void leftNormalizationMatrixFull();
 private:
  network(network const &cpynet);//Copying and assigning networks is better avoided because it would work in a quite unintuitive way (the content of the array structures had to be copied, but the member itself must not be copied) and would be computationally quite expensive - but might be useful if one wanted to genuinely increase D during run (also: add a delete function for manual deletion)
  network& operator=(network const &cpynet);//Use the generate function instead, assignment is dangerous for networks with different parameters 
  //most of these methods are auxiliary functions
  int pctrIndex(int const ai, int const bi, int const aip){return aip+bi*D+ai*D*Dw;}
  int optimize(int const i, double &iolambda);
  int locd(int const i);
  int locDMax(int const i);
  void leftEnrichment(int const i);
  void rightEnrichment(int const i);
  void calcCtrIterLeft(int const position); //iteratively builds up the partial contraction of the left side during a sweep
  void calcCtrIterRight(int const position); //and this does the same for the right side (implementation with two methods is way faster than with one)
  void getPExpressionLeft(int const i, lapack_complex_double *pExpr);
  void getPExpressionRight(int const i, lapack_complex_double *pExpr);
  void getLocalDimensions(int const i);
  mps networkState;
  parameters pars;
  int d,D,L,Dw,nSweeps,icrit;
  int lDL, lDR, ld, lDwR, lDwL;
  pContraction<lapack_complex_double> Lctr;
  pContraction<lapack_complex_double> Rctr;
  //This one is only for consistency checks
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
};

#endif
