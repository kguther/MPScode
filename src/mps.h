#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include <lapacke.h>
#include <complex>

class mps{
 public:
  mps();
  mps(int d, int D, int L);
  ~mps();
  int setParameterD(int Dnew);
  lapack_complex_double& global_access(int i, int si, int ai, int aim){return state_array_access_structure[i][si][ai][aim];}
  void subMatrixStart(lapack_complex_double *&pStart, int i, int si=0){pStart=state_array_access_structure[i][si][0];}
  void generate(int din, int Din, int Lin);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
  int locDimR(int const i);
  int locDimL(int const i);
  int locd(int const i);
  int locDMax(int const i);
 private:
  int d,D,L,icrit;
  lapack_complex_double ****state_array_access_structure;
  void initialize(int din, int Din, int Lin);
  void getIcrit();
};

#endif
