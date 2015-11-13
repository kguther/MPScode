#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include <lapacke.h>
#include <complex>

class mps{
 public:
  mps();
  mps(int d, int D, int L);
  ~mps();
  void mpsCpy(mps &source);
  int setParameterD(int Dnew);
  lapack_complex_double& global_access(int const i, int const si, int const ai, int const aim){return state_array_access_structure[i][si][ai][aim];}
  void subMatrixStart(lapack_complex_double *&pStart, int const i, int const si=0){pStart=state_array_access_structure[i][si][0];}
  void generate(int const din, int const Din, int const Lin);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
  int locDimR(int const i);
  int locDimL(int const i);
  int locd(int const i);
  int locDMax(int const i);
  int maxDim() const {return D;}
  int siteDim() const {return d;}
  int length() const {return L;}
 private:
  mps(mps const &cpymps);
  mps& operator=(mps const &cpymps);
  int d,D,L,icrit;
  lapack_complex_double ****state_array_access_structure;
  void initialize(int const din, int const Din, int const Lin);
  void createInitialState();
  void getIcrit();
};

#endif
