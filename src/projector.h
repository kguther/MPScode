#ifndef NETWORK_PROJECTOR
#define NETWORK_PROJECTOR

#include <lapacke.h>
#include "mps.h"
#include "overlap.h"
#include "siteArray.h"

class projector{
 public:
  projector();
  ~projector();
  void initialize(int const nEigsin);
  void setParameterD(int const Dnew);
  void loadScalarProducts(mps *variationalState, int const iEigen);
  void updateScalarProducts(int const i, int const direction);
  void getProjector(int const i);
  void project(lapack_complex_double *vec, int const i);
  mps *orthoStates;
  overlap *scalarProducts;
  int nCurrentEigen;
 private:
  siteArray<lapack_complex_double> auxiliaryMatrix;
  void getGramMatrix(lapack_complex_double *gram, int const i);
  void getLocalDimensions(int const i);
  int ld, lDL, lDR;
  int nEigs, nRelevantEigens;
  lapack_complex_double *projectionMatrix;
  int vecIndex(int si, int ai, int aim){return aim+ai*lDL+si*lDL*lDR;}
};

#endif
