#ifndef HEFF_MATRIX_CLASS_QN_EXPLOITING
#define HEFF_MATRIX_CLASS_QN_EXPLOITING

#include <vector>
#include "optHMatrix.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"

class blockHMatrix: public optHMatrix{
 public: 
  blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, basisQNOrderMatrix *indexTable, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin);
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  void storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed);
 private:
  std::vector<quantumNumber> *conservedQNsB;
  std::vector<int> blockOffset;
  basisQNOrderMatrix *indexTable;
  int vecBlockIndex(int const iBlock, int const j, int const k){return k+j*indexTable->lBlockSizeLP(i,iBlock)+blockOffset[iBlock];}
  void excitedStateProject(arcomplex<double> *v, int const i);
  void storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded);
};

#endif
