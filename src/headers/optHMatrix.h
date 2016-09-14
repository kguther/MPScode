#ifndef HEFF_MATRIX_CLASS
#define HEFF_MATRIX_CLASS

#include <vector>
#include "projector.h"
#include "dimensionTable.h"
#include "quantumNumber.h"
#include "templates/mpo.h"
#include "mpstype.h"

//---------------------------------------------------------------------------------------------------//
// This class is used as an interface to provide the matrix-vector product for ARPACK++
// It also contains a projector for excited state search. The projection does nothing in ground state
// search.
//---------------------------------------------------------------------------------------------------//

class optHMatrix{
 public:
  optHMatrix(mpsEntryType *R, mpsEntryType *L, mpo<mpsEntryType > const *Hin, dimensionTable const &dimInfo, int Dwin, int iIn, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin);
  virtual ~optHMatrix(){}
  virtual void MultMv(mpsEntryType *v, mpsEntryType *w);
  void MultMvQNConserving(mpsEntryType *v, mpsEntryType *w);
  virtual int dim() const {return dimension;}
 protected:
  mpsEntryType *Rctr;
  mpsEntryType *Lctr;
  mpo<mpsEntryType > const *HMPO;
  double shift;
  projector *P;
  std::vector<quantumNumber> *conservedQNs;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,dimension,i;
  void projectQN(mpsEntryType *v);
  int ctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*Dw*D;}//Partial contractions are of uniform size, so no usage of local dimension here
  int hIndex(int const si, int const sip, int const bi, int const bim) {return bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si;}//same is true for MPO
  int vecIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
};

#endif

