#ifndef HEFF_MATRIX_CLASS_QN_EXPLOITING
#define HEFF_MATRIX_CLASS_QN_EXPLOITING

#include <vector>
#include "optHMatrix.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"

//---------------------------------------------------------------------------------------------------//
// The blockHMatrix class is an expansion of the optHMatrix class which is designed to exploit good
// quantum numbers. Therefore, it accesses an index table (i.e. a basisQNOrderMatrix object) containing
// the block structure of the MPS matrices. Usually, this is the index table of the corresponding MPS
//---------------------------------------------------------------------------------------------------//

class blockHMatrix: public optHMatrix{
 public: 
  blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, int sweepDirectionIn, basisQNOrderMatrix *indexTable, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin);
  ~blockHMatrix();
  void MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w);
  void storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed);
  void storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded);
  void prepareInput(arcomplex<double> *startingVector);
  void readOutput(arcomplex<double> *outputVector);
  arcomplex<double> *compressedVector;
 private:
  std::vector<quantumNumber> *conservedQNsB;
  std::vector<int> blockOffset;
  int sweepDirection;
  basisQNOrderMatrix *indexTable;
  int vecBlockIndex(int const iBlock, int const j, int const k){return k+j*indexTable->lBlockSizeLP(i,iBlock)+blockOffset[iBlock];}
  void excitedStateProject(arcomplex<double> *v, int const i);
};

#endif
