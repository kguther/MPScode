#ifndef HEFF_MATRIX_CLASS
#define HEFF_MATRIX_CLASS

#include <arcomp.h>
#include "projector.h"
#include "parameters.h"
#include "quantumNumber.h"

//---------------------------------------------------------------------------------------------------//
// This class is used as an interface to provide the matrix-vector product for ARPACK++
// It also contains a projector for excited state search. The projection does nothing in ground state
// search.
//---------------------------------------------------------------------------------------------------//

class optHMatrix{
 public:
  optHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, problemParameters pars, int Din, int iIn, projector *excitedStateP, double shift, int const nQNsin=0, quantumNumber *conservedQNsin=0);
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  void MultMvQNConserving(arcomplex<double> *v, arcomplex<double> *w);
  int dim() const {return dimension;}
 private:
  arcomplex<double> *Rctr;
  arcomplex<double> *Lctr;
  arcomplex<double> *H;
  double shift;
  projector *P;
  quantumNumber *conservedQNs;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,icrit,dimension,nQNs,i;
  void projectQN(arcomplex<double> *v);
  int ctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*Dw*D;}//Partial contractions are of uniform size, so no usage of local dimension here
  int hIndex(int const si, int const sip, int const bi, int const bim) {return bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si;}//same is true for MPO
  int vecIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
};

#endif

