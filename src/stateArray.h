#ifndef MATRIX_PRODUCT_STATE_ARRAY
#define MATRIX_PRODUCT_STATE_ARRAY

#include <vector>
#include "mkl_complex_defined.h"
#include "arraycreation.h"
#include "dimensionTable.h"
#include "localHSpaces.h"
#include "baseTensor.h"

class stateArray{
 public:
  stateArray();
  stateArray(dimensionTable const &dimInfoIn);
  virtual void mpsCpy(stateArray const &source);
  virtual int setParameterD(int Dnew);
  virtual int setParameterL(int Lnew);
  const lapack_complex_double& operator() (int i, int si, int ai, int aim) const{std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  lapack_complex_double& operator()(int i, int si, int ai, int aim){std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  lapack_complex_double* operator()(int i, int si=0){lapack_complex_double *target; stateArrayAccessStructure[i].getPtr(target,si); return target;}
  lapack_complex_double& global_access(int i, int si, int ai, int aim){std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  const lapack_complex_double& global_access(int i, int si, int ai, int aim) const {std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  void subMatrixStart(lapack_complex_double *&pStart, int i, int si=0){stateArrayAccessStructure[i].getPtr(pStart,si);}
  void subMatrixStart(lapack_complex_double const*&pStart, int i, int si=0)const {stateArrayAccessStructure[i].getPtr(pStart,si);}
  int locDimR(int i) const;
  int locDimL(int i) const;
  int locd(int i) const;
  int locDMax(int i) const;
  int maxDim() const {return dimInfo.D();}
  int siteDim() const {return dimInfo.d();}
  int length() const {return dimInfo.L();}
  virtual void initialize(dimensionTable const &dimInfoIn);
  dimensionTable dimInfo;
 private:
  stateArray(stateArray const &cpyState);
  stateArray& operator=(stateArray const &cpyState);
 protected:
  int D,L;
  std::vector<baseTensor<lapack_complex_double> > stateArrayAccessStructure;
  void createStateArray(int const Lin);
  std::vector<int> getIndexVec(int si, int ai, int aim) const;
};

#endif
