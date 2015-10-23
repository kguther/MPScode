#include <arcomp.h>
#include "parameters.h"

class optHMatrix{
 public:
  optHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, parameters pars, int i);
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  int dim() const {return dimension;}
 private:
  arcomplex<double> *Rctr;
  arcomplex<double> *Lctr;
  arcomplex<double> *H;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,icrit,dimension;
  int ctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*Dw*D;}//Partial contractions are of uniform size, so no usage of local dimension here
  int hIndex(int const si, int const sip, int const bi, int const bim) {return bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si;}//same is true for MPO
  int containerIndexInner(int const si, int const aim, int const ai, int const bi) {return bi+lDwR*ai+lDwR*lDR*aim+lDwR*lDR*lDL*si;}
  int containerIndexOuter(int const si, int const bim, int const ai, int const aim) {return aim+ai*lDL+bim*lDL*lDR+si*lDwL*lDL*lDR;}
  int vecIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
};

