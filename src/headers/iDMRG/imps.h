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
  virtual ~imps(){}
  imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin);
  imps(mps const &source);
  virtual int addSite(int Lnew, int i);
  void exportState(mps &target);
  virtual void subMatrixStart(arcomplex<double> *&pStart, int i, int si=0){mps::subMatrixStart(pStart,i,si);}
  virtual int refineQN(int i, std::vector<std::complex<int> > const &leftSideLabels, std::vector<std::complex<int> > const &rightSideLabels, std::vector<std::complex<int> > const &targetQN);
  virtual int currentSite()const {return dimInfo.L()/2-1;}
  virtual int internalSite()const {return currentSite();}
  virtual twositeQNOrderMatrix const& centralIndexTable() const{return centralIndexTableVar;}
  basisQNOrderMatrix const& indexTable() const{return indexTableVar;}
  virtual siteQNOrderMatrix const& leftIndexTable()const {return indexTableVar.getLocalIndexTable(currentSite());}
virtual siteQNOrderMatrix const& rightIndexTable()const {return indexTableVar.getLocalIndexTable(1+currentSite());}
  virtual pseudoQuantumNumber* getConservedQNs(int iQN){return &(conservedQNs[iQN]);}
  virtual baseTensor<arcomplex<double> > const& getSiteTensor(int i){return getStateArrayEntry(i);}
 private:
  twositeQNOrderMatrix centralIndexTableVar;
  std::vector<std::complex<int> > internalTargetQNBuffer, internalSourceBuffer;
};

#endif
