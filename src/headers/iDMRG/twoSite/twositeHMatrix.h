#ifndef HEFF_MATRIX_CLASS_TWOSITE_QN_EXPLOITING
#define HEFF_MATRIX_CLASS_TWOSITE_QN_EXPLOITING

#include <arcomp.h>
#include <vector>
#include "mpo.h"
#include "mps.h"
#include "twositeQNOrderMatrix.h"
#include "dimensionTable.h"
#include "projector.h"

//---------------------------------------------------------------------------------------------------//
// STORAGE SCHEME FOR TWOSITE MATRICES: (sip,si,air,aim)
//---------------------------------------------------------------------------------------------------//

class twositeHMatrix{
 public:
  twositeHMatrix(arcomplex<double> *R, arcomplex<double> *L, mpo<arcomplex<double> > *Hin, int HPosL, int HPosR, dimensionTable const &dimInfo, twositeQNOrderMatrix const* indexTable, double shift=0.0, projector *excitedStateP=0);
  void MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w);
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  void readOutput(arcomplex<double> *outputVector);
  void prepareInput(arcomplex<double> *inputVector);
  void storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed);
  void storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded);
  int dim()const {return dimension;}
  int dimFull()const {return ld*ld*lDL*lDRR;}
  arcomplex<double>* getCompressedVector() {return &(compressedVector[0]);}
 private:
  int dimension;
  int i, ld, ldp, lDL, lDRR;
  int lDwL, lDwRR;
  int D, Dw, HMatrixPosL, HMatrixPosR;
  double shift;
  std::vector<int> blockOffset;
  dimensionTable dimInfo;
  mpo<arcomplex<double> > *HMPO;
  arcomplex<double> *Lctr, *Rctr;
  std::vector<arcomplex<double> > compressedVector, W;
  twositeQNOrderMatrix const* indexTable;
  int vecBlockIndex(int iBlock, int j, int k){return k+j*indexTable->lBlockSize(iBlock)+blockOffset[iBlock];}
  int vecIndex(int sip, int si, int air, int aim){return aim+si*lDL*lDRR+air*lDL+sip*ld*lDL*lDRR;}
  int hIndex(int si, int sip, int sit, int sitp, int bir, int bim){return bim+bir*lDwL+sitp*lDwL*lDwRR+sit*ldp*lDwL*lDwRR+sip*ldp*ldp*lDwL*lDwRR+si*lDwL*lDwRR*ld*ldp*ldp;}
  int ctrIndex(int ai, int bi, int aip){return aip+bi*D+ai*D*Dw;}
};

#endif
