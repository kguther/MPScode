#include "blockHMatrix.h"

blockHMatrix::  blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, basisQNOrderMatrix *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein),
  conservedQNs(conservedQNsin)
{}


