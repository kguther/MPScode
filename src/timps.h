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
  virtual int addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN, std::vector<std::complex<int> > const &source);
  virtual int refineQN(int i, std::vector<std::complex<int> > const& newQN);
  virtual int currentSite() const{return (dimInfoBase.L()-1)/2;}
  virtual int internalSite() const{return position;}
  virtual twositeQNOrderMatrix const& centralIndexTable() const{return centralIndexTableVar;}
  virtual basisQNOrderMatrix const& indexTable() const{return indexTableVar;}
  virtual pseudoQuantumNumber* getConservedQNs(int iQN){return &(conservedQNs[iQN]);}
  virtual baseTensor<std::complex<double> > const& getSiteTensor(int i){return unitCell[convertPosition(i)];}
 private:
  int unitCellSize;
  twositeQNOrderMatrix centralIndexTableVar;
  basisQNOrderMatrix indexTableVar;
  std::vector<manualQuantumNumber> conservedQNs;
  std::vector<baseTensor<std::complex<double> > > unitCell;
  int position;
  int convertPosition(int i){return position-currentSite()+i;}
  void setUpTables();
  dimensionTable internalDimInfo;
};

#endif
