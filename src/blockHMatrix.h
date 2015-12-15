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
 private:
  std::vector<quantumNumber> *conservedQNs;
  basisQNOrderMatrix *indexTable;

};

#endif
