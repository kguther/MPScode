#ifndef MATRIX_PRODUCT_STATE_ARRAY
#define MATRIX_PRODUCT_STATE_ARRAY

#include <complex>
#include <vector>
#include "dimensionTable.h"
#include "localHSpaces.h"
#include "templates/baseTensor.h"

//---------------------------------------------------------------------------------------------------//
// Base class for MPS-type objects. The state array only contains an array of tensors with global 
// access functions, submatrix access functions and a dimension table.
//---------------------------------------------------------------------------------------------------//

class stateArray{
 public:
  stateArray();
  stateArray(dimensionTable const &dimInfoIn);
  virtual int setParameterD(int Dnew);
  virtual int setParameterL(int Lnew);
  //The access operators grant direct access to matrix entries, but are very slow, they should only be used in uncritical functions.
  const std::complex<double>& operator() (int i, int si, int ai, int aim) const{std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  std::complex<double>& operator()(int i, int si, int ai, int aim){std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  std::complex<double>* operator()(int i, int si=0){std::complex<double> *target; stateArrayAccessStructure[i].getPtr(target,si); return target;}
  std::complex<double>& global_access(int i, int si, int ai, int aim){std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  const std::complex<double>& global_access(int i, int si, int ai, int aim) const {std::vector<int> indices=getIndexVec(si,ai,aim); return stateArrayAccessStructure[i](indices);}
  //The subMatrixStart functions set pointers to the internal subarrays of the tensor array. Use for fast access to matrix elements of the stateArray (much faster than global_access)
  virtual void subMatrixStart(std::complex<double> *&pStart, int i, int si=0){stateArrayAccessStructure[i].getPtr(pStart,si);}
  virtual void subMatrixStart(std::complex<double> const*&pStart, int i, int si=0)const {stateArrayAccessStructure[i].getPtr(pStart,si);}
  void setStateArray(stateArray const &source){stateArrayAccessStructure=source.stateArrayAccessStructure;}
  int locDimR(int i) const;
  int locDimL(int i) const;
  int locd(int i) const;
  int locDMax(int i) const;
  int maxDim() const {return dimInfo.D();}
  int siteDim() const {return dimInfo.d();}
  int length() const {return dimInfo.L();}
  virtual void initialize(dimensionTable const &dimInfoIn);
  dimensionTable const& getDimInfo()const {return dimInfo;}
 private:
  std::vector<baseTensor<std::complex<double> > > stateArrayAccessStructure;
  void createStateArray(int const Lin);
  std::vector<int> getIndexVec(int si, int ai, int aim) const;
 protected:
  int D,L;
  dimensionTable dimInfo;
  baseTensor<std::complex<double> > const &getStateArrayEntry(int i){return stateArrayAccessStructure[i];}
  void setStateArrayEntry(int i, baseTensor<std::complex<double> > const &source){stateArrayAccessStructure[i]=source;}
};

#endif
