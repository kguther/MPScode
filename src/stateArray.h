#ifndef MATRIX_PRODUCT_STATE_ARRAY
#define MATRIX_PRODUCT_STATE_ARRAY

#include "mkl_complex_defined.h"
#include "arraycreation.h"
#include "dimensionTable.h"
#include "localHSpaces.h"

class stateArray{
 public:
  stateArray();
  stateArray(dimensionTable const &dimInfoIn);
  virtual ~stateArray();
  virtual void mpsCpy(stateArray const &source);
  int setParameterD(int Dnew);
  lapack_complex_double& global_access(int i, int si, int ai, int aim){return state_array_access_structure[i][si][ai][aim];}
  const lapack_complex_double& global_access(int i, int si, int ai, int aim) const {return state_array_access_structure[i][si][ai][aim];}
  void subMatrixStart(lapack_complex_double *&pStart, int i, int si=0){pStart=state_array_access_structure[i][si][0];}
  virtual void generate(dimensionTable const &dimInfoIn);
  int locDimR(int i) const;
  int locDimL(int i) const;
  int locd(int i) const;
  int locDMax(int i) const;
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
  virtual void initialize(dimensionTable const &dimInfoIn);
  void getIcrit();
  void createStateArray(int const Din, int const Lin, lapack_complex_double *****array);
};

#endif
