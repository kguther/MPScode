#include <arcomp.h>
#include "projector.h"
#include "parameters.h"

//---------------------------------------------------------------------------------------------------//
// This class is used as an interface to provide the matrix-vector product for ARPACK++
// It also contains a projector for excited state search. The projection does nothing in ground state
// search.
//---------------------------------------------------------------------------------------------------//

class optHMatrix{
 public:
  optHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, problemParameters pars, int Din, int i, projector *excitedStateP, double shift);
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  int dim() const {return dimension;}
 private:
  arcomplex<double> *Rctr;
  arcomplex<double> *Lctr;
  arcomplex<double> *H;
  projector *P;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,icrit,dimension,currentSite;
  double shift;
  int ctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*Dw*D;}//Partial contractions are of uniform size, so no usage of local dimension here
  int hIndex(int const si, int const sip, int const bi, int const bim) {return bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si;}//same is true for MPO
  int vecIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
};

