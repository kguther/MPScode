#ifndef MATRIX_PRODUCT_STATE_ARRAY
#define MATRIX_PRODUCT_STATE_ARRAY

#include "mkl_complex_defined.h"
#include "arraycreation.h"
#include "dimensionTable.h"

class stateArray{
 public:
  stateArray();
  stateArray(int const din, int const Din, int const Lin);
  ~stateArray();
  void mpsCpy(stateArray &source);
  int setParameterD(int Dnew);
  lapack_complex_double& global_access(int const i, int const si, int const ai, int const aim){return state_array_access_structure[i][si][ai][aim];}
  void subMatrixStart(lapack_complex_double *&pStart, int const i, int const si=0){pStart=state_array_access_structure[i][si][0];}
  virtual void generate(int const din, int const Din, int const Lin);
  int locDimR(int const i);
  int locDimL(int const i);
  int locd(int const i);
  int locDMax(int const i);
  int maxDim() const {return D;}
  int siteDim() const {return d;}
  int length() const {return L;}
  dimensionTable dimInfo;
 private:
  stateArray(stateArray const &cpyState);
  stateArray& operator=(stateArray const &cpyState);
 protected:
  int d,D,L,icrit;
  lapack_complex_double ****state_array_access_structure;
  virtual void initialize(int const din, int const Din, int const Lin);
  void getIcrit();
  void createStateArray(int const din, int const Din, int const Lin, lapack_complex_double *****array);
};

#endif
