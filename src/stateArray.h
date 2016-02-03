#ifndef MATRIX_PRODUCT_STATE_ARRAY
#define MATRIX_PRODUCT_STATE_ARRAY

#include "mkl_complex_defined.h"
#include "arraycreation.h"
#include "dimensionTable.h"
#include "localHSpaces.h"

class stateArray{
 public:
  stateArray();
  stateArray(dimensionTable &dimInfoIn);
  ~stateArray();
  virtual void mpsCpy(stateArray &source);
  int setParameterD(int Dnew);
  lapack_complex_double& global_access(int const i, int const si, int const ai, int const aim){return state_array_access_structure[i][si][ai][aim];}
  void subMatrixStart(lapack_complex_double *&pStart, int const i, int const si=0){pStart=state_array_access_structure[i][si][0];}
  virtual void generate(dimensionTable &dimInfoIn);
  int locDimR(int const i);
  int locDimL(int const i);
  int locd(int const i);
  int locDMax(int const i);
  int maxDim() const {return dimInfo.D();}
  int siteDim() const {return dimInfo.d();}
  int length() const {return dimInfo.L();}
  dimensionTable dimInfo;
 private:
  stateArray(stateArray const &cpyState);
  stateArray& operator=(stateArray const &cpyState);
 protected:
  int D,L;
  int icrit;
  lapack_complex_double ****state_array_access_structure;
  virtual void initialize(dimensionTable &dimInfoIn);
  void getIcrit();
  void createStateArray(int const Din, int const Lin, lapack_complex_double *****array);
};

#endif
