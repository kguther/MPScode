#ifndef HEFF_MATRIX_CLASS
#define HEFF_MATRIX_CLASS

#include <complex>
#include <vector>
#include "projector.h"
#include "dimensionTable.h"
#include "quantumNumber.h"
#include "templates/mpo.h"

//---------------------------------------------------------------------------------------------------//
// This class is used as an interface to provide the matrix-vector product for ARPACK++
// It also contains a projector for excited state search. The projection does nothing in ground state
// search.
//---------------------------------------------------------------------------------------------------//

class optHMatrix{
 public:
  optHMatrix(std::complex<double> *R, std::complex<double> *L, mpo<std::complex<double> > const *Hin, dimensionTable const &dimInfo, int Dwin, int iIn, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin);
  virtual ~optHMatrix(){}
  virtual void MultMv(std::complex<double> *v, std::complex<double> *w);
  void MultMvQNConserving(std::complex<double> *v, std::complex<double> *w);
  virtual int dim() const {return dimension;}
 protected:
  std::complex<double> *Rctr;
  std::complex<double> *Lctr;
  mpo<std::complex<double> > const *HMPO;
  double shift;
  projector *P;
  std::vector<quantumNumber> *conservedQNs;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,dimension,i;
  void projectQN(std::complex<double> *v);
  int ctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*Dw*D;}//Partial contractions are of uniform size, so no usage of local dimension here
  int hIndex(int const si, int const sip, int const bi, int const bim) {return bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si;}//same is true for MPO
  int vecIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
};

#endif

