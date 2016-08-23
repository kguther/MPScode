#ifndef HEFF_MATRIX_CLASS_QN_EXPLOITING
#define HEFF_MATRIX_CLASS_QN_EXPLOITING

#include <vector>
#include "optHMatrix.h"
#include "siteQNOrderMatrix.h"
#include "templates/tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// The blockHMatrix class is an expansion of the optHMatrix class which is designed to exploit good
// quantum numbers. Therefore, it accesses an index table (i.e. a basisQNOrderMatrix object) containing
// the block structure of the MPS matrices. Usually, this is the index table of the corresponding MPS
//---------------------------------------------------------------------------------------------------//

//In case of boredom: remove inheritance from blockHMatrix class - this only makes things complicated

class blockHMatrix: private optHMatrix{
 public: 
  blockHMatrix(std::complex<double> *R, std::complex<double> *L, mpo<std::complex<double> > const *Hin, dimensionTable &dimInfo, int Dwin, int iIn, siteQNOrderMatrix const *indexTable, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin, int const cached=1);
  //this is the essential function of the blockMatrix
  void MultMvBlocked(std::complex<double> *v, std::complex<double> *w);
  //LP means left pairing, because the sigma-index is grouped with the left index of the matrix
  void MultMvBlockedLP(std::complex<double> *v, std::complex<double> *w);
  void storageCompress(std::complex<double> *v, std::complex<double> *vCompressed);
  void storageExpand(std::complex<double> *v, std::complex<double> *vExpanded);
  void prepareInput(std::complex<double> *startingVector);
  void readOutput(std::complex<double> *outputVector);
  //the sparse version ONLY works for nearest-neighbour hamiltonians as it is optimized for their MPO representation
  void buildSparseHBlocked();
  virtual int dim()const {return dimension;}
  std::complex<double>* getCompressedVector() {return &(compressedVector[0]);}
 private:
  int explicitMv;
  std::vector<int> blockOffset;
  std::vector<std::complex<double> > sparseMatrix;
  std::vector<std::complex<double> > compressedVector;
  std::vector<int> rowPtr, colIndices;
  tmpContainer<std::complex<double> > innerContainer;
  tmpContainer<std::complex<double> > outerContainer;
  siteQNOrderMatrix const *indexTable;
  int vecBlockIndexLP(int const iBlock, int const j, int const k){return k+j*indexTable->lBlockSizeLP(iBlock)+blockOffset[iBlock];}
  void excitedStateProject(std::complex<double> *v);
  std::complex<double> HEffEntry(int const si, int const aim, int const ai, int const sip, int const aimp, int const aip);
};

#endif
