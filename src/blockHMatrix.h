#ifndef HEFF_MATRIX_CLASS_QN_EXPLOITING
#define HEFF_MATRIX_CLASS_QN_EXPLOITING

#include <vector>
#include "optHMatrix.h"
#include "basisQNOrderMatrix.h"

//---------------------------------------------------------------------------------------------------//
// The blockHMatrix class is an expansion of the optHMatrix class which is designed to exploit good
// quantum numbers. Therefore, it accesses an index table (i.e. a basisQNOrderMatrix object) containing
// the block structure of the MPS matrices. Usually, this is the index table of the corresponding MPS
//---------------------------------------------------------------------------------------------------//

class blockHMatrix: private optHMatrix{
 public: 
  blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, mpo<arcomplex<double> > *Hin, dimensionTable &dimInfo, int Dwin, int iIn, basisQNOrderMatrix const *indexTable, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin, int const cached=1);
  void MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w);
  void MultMvBlockedLP(arcomplex<double> *v, arcomplex<double> *w);
  void storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed);
  void storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded);
  void prepareInput(arcomplex<double> *startingVector);
  void readOutput(arcomplex<double> *outputVector);
  void buildSparseHBlocked();
  virtual int dim()const {return dimension;}
  arcomplex<double>* getCompressedVector() {return &(compressedVector[0]);}
 private:
  int explicitMv;
  std::vector<int> blockOffset;
  std::vector<arcomplex<double> > sparseMatrix;
  std::vector<arcomplex<double> > compressedVector;
  std::vector<int> rowPtr, colIndices;
  mpo<arcomplex<double> > *HMPO;
  basisQNOrderMatrix const *indexTable;
  int vecBlockIndexLP(int const iBlock, int const j, int const k){return k+j*indexTable->lBlockSizeLP(i,iBlock)+blockOffset[iBlock];}
  void excitedStateProject(arcomplex<double> *v);
  arcomplex<double> HEffEntry(int const si, int const aim, int const ai, int const sip, int const aimp, int const aip);
};

#endif
