#ifndef NETWORK
#define NETWORK

#include <complex>
#include <vector>
#include <lapacke.h>
#include <arcomp.h>
#include "parameters.h"
#include "mpo.h"
#include "pContraction.h"

//---------------------------------------------------------------------------------------------------//
// The network class contains all information required for a run of the simulation, that is, the whole
//---------------------------------------------------------------------------------------------------//

class network{
 public:
  network();
  network(parameters inputpars);
  ~network();
  double solve();
  void generate(parameters inputpars);
  void setParameterNSweeps(int Nnew);
  int setParameterD(int Dnew);
  void leftNormalizationMatrixFull();
  int calcCtrFull(const int direction);
  lapack_complex_double ****networkState;
  mpo<lapack_complex_double> networkH;
  pContraction<lapack_complex_double> Lctr;
  pContraction<lapack_complex_double> Rctr;
 private:
  network(network const &cpynet);//Copying and assigning networks is better avoided because it would work in a quite unintuitive way (the content of the array structures had to be copied, but the member itself must not be copied) and would be computationally quite expensive - but might be useful if one wanted to genuinely increase D during run (also: add a delete function for manual deletion)
  network& operator=(network const &cpynet);//Use the generate function instead, assignment is dangerous for networks with different parameters 
  void initialize(parameters pars);
  int optimize(int const i, double &iolambda);
  int locDimR(int const i);
  int locDimL(int const i);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  int pctrIndex(int ai, int bi, int aip){return aip+bi*D+ai*D*Dw;}
  void normalizeFinal(int const i);
  void calcCtrIterLeft(const int position); //iteratively builds up the partial contraction of the left side during a sweep
  void calcCtrIterRight(const int position); //and this does the same for the right side (implementation with two methods is way faster than with one)
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
  parameters pars;
  int d,D,L,Dw,nSweeps,icrit;
};

#endif
