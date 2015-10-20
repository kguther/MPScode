#ifndef SINGLE_SITE_OPTIMIZER_EXTERNAL_ARPACK
#define SINGLE_SITE_OPTIMIZER_EXTERNAL_ARPACK

#include <arcomp.h>

class siteoptimizer{
 public:
  siteoptimizer(int dimensionin, int nzelin, arcomplex<double> *Hin, int *irowin, int *pcolin);
  void solveEigen(arcomplex<double> *plambda, arcomplex<double> *currentM);
 private:
  int *irow;
  int *pcol;
  int nzel, dimension;
  arcomplex<double> *H;
};

#endif
