#ifndef TRANSLATION_INVARIANT_MPS
#define TRANSLATION_INVARIANT_MPS

//---------------------------------------------------------------------------------------------------//
// MPS class for memory-efficient iDMRG targeting translationally invariant states.
//---------------------------------------------------------------------------------------------------//

#include "baseTensor.h"
#include "impBase.h"
#include "dimensionTable.h"
#include "manualQuantumNumber.h"
#include <vector>
#include <complex>
#include "twositeQNOrderMatrix.h"
#include "basisQNOrderMatrix.h"

class timps: public impBase{
 public:
 timps():unitCellSize(0){}
  timps(int unitCellSizeIn, dimensionTable const &dimInfoBaseIn, std::vector<quantumNumber> const &conservedQNsIn);
  virtual void subMatrixStart(std::complex<double> *&pStart, int i, int si=0);
  virtual int addSite(int Lnew, int i);
  virtual int refineQN(int i, std::vector<std::complex<int> > const &leftSideLabels, std::vector<std::complex<int> > const &rightSideLabels, std::vector<std::complex<int> > const &targetQN);
  virtual int currentSite() const{return (dimInfoBase.L()-1)/2;}
  virtual int internalSite() const{return position;}
  virtual twositeQNOrderMatrix const& centralIndexTable() const{return centralIndexTableVar;}
  virtual siteQNOrderMatrix const& leftIndexTable() const{return leftTable;}
  virtual siteQNOrderMatrix const& rightIndexTable() const{return rightTable;}
  virtual pseudoQuantumNumber* getConservedQNs(int iQN){return &(conservedQNs[iQN]);}
  virtual baseTensor<std::complex<double> > const& getSiteTensor(int i){return unitCell[convertPosition(i)];}
 private:
  int unitCellSize;
  twositeQNOrderMatrix centralIndexTableVar;
  siteQNOrderMatrix leftTable, rightTable;
  std::vector<manualQuantumNumber> conservedQNs;
  std::vector<manualQuantumNumber> reducedIndexLabelQNs;
  std::vector<baseTensor<std::complex<double> > > unitCell;
  int position;
  int convertPosition(int i){return position-currentSite()+i;}
  void setUpSingleSiteTables();
  void setUpTwositeTable();
  dimensionTable internalDimInfo;
};

#endif
