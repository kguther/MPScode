#ifndef MATRIX_PRODUCT_STATE_VARIABLE_SIZE
#define MATRIX_PRODUCT_STATE_VARIABLE_SIZE

#include "mps.h"
#include "twositeQNOrderMatrix.h"
#include "impBase.h"

//---------------------------------------------------------------------------------------------------//
// The iMPS is a MPS with the additional functionality of growing, making it suited for the iDMRG
// algorithm. Also, it supports a two-site optimization.
//---------------------------------------------------------------------------------------------------//


class imps: public impBase, public mps{
 public:
  imps();
  imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin);
  virtual void addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN);
  void exportState(mps &target);
  void importState(mps const &source);
  virtual void subMatrixStart(lapack_complex_double *&pStart, int i, int si=0){mps::subMatrixStart(pStart,i,si);}
  virtual int refineQN(int i, std::vector<std::complex<int> > const &source);
  virtual int currentSite()const {return dimInfo.L()/2-1;}
  virtual twositeQNOrderMatrix const& centralIndexTable(){return centralIndexTableVar;}
 private:
  twositeQNOrderMatrix centralIndexTableVar;
};

#endif
