#include <arcomp.h>
#include "parameters.h"

class optHMatrix{
 public:
  optHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, parameters pars, int i);  //Use R[i][0][0] and local dimensions as input instead of site index
  void MultMv(arcomplex<double> *v, arcomplex<double> *w);
  int dim() const {return dimension;}
 private:
  arcomplex<double> *Rctr;
  arcomplex<double> *Lctr;
  arcomplex<double> *H;
  int d,D,L,Dw,lDR,lDL,lDwR,lDwL,icrit,dimension;
  int ctrIndex(int const ai, int const bi, int const aip);
  int hIndex(int const si, int const sip, int const bi, int const bim);
  int containerIndexInner(int const si, int const aim, int const ai, int const bi);
  int containerIndexOuter(int const si, int const bim, int const ai, int const aim);
  int vecIndex(int const si, int const ai, int const aim);
};
