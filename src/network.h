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
  double optimize(int i);
  int locDimR(int i);
  int locDimL(int i);
  int leftNormalizeState(int i);
  int rightNormalizeState(int i);
  void normalizeFinal(int i);
  int calcCtrFull(lapack_complex_double ****Pctr, const int direction);//computes partial contractions for preparation of a sweep
  void calcCtrIter(lapack_complex_double ****Pctr, const int direction, const int position); //iteratively builds up the partial contraction during a sweep (for the side not initialized) start at position==(L-1) for R expression and position==0 for L expression
  lapack_complex_double ****networkState;
  lapack_complex_double *****networkH;
 private:
  network();
  network(network const &cpynet);
  network& operator=(network const &cpynet);
  void buildMatrix(int i, vector<int> *irow, vector<int> *pcol, int *nzel, vector<arcomplex<double> > *H);
  parameters pars;
  int d,D,L,Dw,N,icrit;
  lapack_complex_double ****Lctr;
  lapack_complex_double ****Rctr;
};

#endif
