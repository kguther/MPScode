#ifndef NETWORK
#define NETWORK

#include <complex>
#include <vector>
#include <lapacke.h>
#include <arcomp.h>
#include "parameters.h"

using namespace std;

class network{
 public:
  network(parameters inputpars);
  ~network();
  double solve();
  int locDimR(int const i);
  int locDimL(int const i);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
  int calcCtrFull(lapack_complex_double ****Pctr, const int direction);//computes partial contractions for preparation of a sweep
  void calcCtrIter(lapack_complex_double ****Pctr, const int direction, const int position); //iteratively builds up the partial contraction during a sweep (for the side not initialized) start at position==(L-1) for R expression and position==0 for L expression
  lapack_complex_double ****networkState;
  lapack_complex_double *****networkH;
 private:
  network();
  network(network const &cpynet);
  network& operator=(network const &cpynet);
  double optimize(int const i);
  parameters pars;
  int d,D,L,Dw,N,icrit;
  lapack_complex_double ****Lctr;
  lapack_complex_double ****Rctr;
};

#endif
