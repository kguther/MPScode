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
  void leftNormalizationMatrixFull();
  int calcCtrFull(lapack_complex_double ****Pctr, const int direction);//computes partial contractions for preparation of a sweep
  lapack_complex_double ****networkState;
  lapack_complex_double *****networkH;
 private:
  network();                              //Creating a network without parameters does not make much sense
  network(network const &cpynet);         //Copying and assigning networks is better avoided because it would work in a quite unintuitive way (the content of the array structures had to be copied, but the member itself must not be copied) and would be computationally quite expensive
  network& operator=(network const &cpynet);
  int optimize(int const i, double *iolambda);
  void calcCtrIterLeft(lapack_complex_double ****Pctr, const int position); //iteratively builds up the partial contraction of the left side during a sweep
  void calcCtrIterRight(lapack_complex_double ****Pctr, const int position); //and this does the same for the right side (implementation with two methods is way faster than with one)
  void leftNormalizationMatrixIter(int i, lapack_complex_double *psi);
  parameters pars;
  int d,D,L,Dw,N,icrit;
  lapack_complex_double ****Lctr;
  lapack_complex_double ****Rctr;
};

#endif
