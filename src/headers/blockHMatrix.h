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
  blockHMatrix(mpsEntryType *R, mpsEntryType *L, mpo<mpsEntryType > const *Hin, dimensionTable &dimInfo, int Dwin, int iIn, siteQNOrderMatrix const *indexTable, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin, int const cached=1);
  //this is the essential function of the blockMatrix
  void MultMvBlocked(mpsEntryType *v, mpsEntryType *w);
  //LP means left pairing, because the sigma-index is grouped with the left index of the matrix
  void MultMvBlockedLP(mpsEntryType *v, mpsEntryType *w);
  void storageCompress(mpsEntryType *v, mpsEntryType *vCompressed);
  void storageExpand(mpsEntryType *v, mpsEntryType *vExpanded);
  void prepareInput(mpsEntryType *startingVector);
  void readOutput(mpsEntryType *outputVector);
  //the sparse version ONLY works for nearest-neighbour hamiltonians as it is optimized for their MPO representation
  void buildSparseHBlocked();
  virtual int dim()const {return dimension;}
  mpsEntryType* getCompressedVector() {return &(compressedVector[0]);}
 private:
  int explicitMv;
  std::vector<int> blockOffset;
  std::vector<mpsEntryType > sparseMatrix;
  std::vector<mpsEntryType > compressedVector;
  std::vector<int> rowPtr, colIndices;
  tmpContainer<mpsEntryType > innerContainer;
  tmpContainer<mpsEntryType > outerContainer;
  siteQNOrderMatrix const *indexTable;
  int vecBlockIndexLP(int const iBlock, int const j, int const k){return k+j*indexTable->lBlockSizeLP(iBlock)+blockOffset[iBlock];}
  void excitedStateProject(mpsEntryType *v);
  mpsEntryType HEffEntry(int const si, int const aim, int const ai, int const sip, int const aimp, int const aip);
};

#endif
