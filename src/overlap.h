#ifndef MPS_OVERLAP
#define MPS_OVERLAP

#include <lapacke.h>
#include "mps.h"

class overlap{
 public:
  overlap();
  ~overlap();
  void loadMPS(mps *psi, mps *phi);
  void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  lapack_complex_double& Lctr_access(int const i, int const aim, int const aimp){return Lctr[aimp+aim*D+i*D*D];}
  lapack_complex_double& Rctr_access(int const i, int const ai, int const aip){return Rctr[aip+ai*D+i*D*D];}
 private:
  int L, int D, int d;
  mps *psi;
  mps *phi;
  lapack_complex_double *Lctr;
  lapack_complex_double *Rctr;
  void subContractionStartLeft(lapack_complex_double *&pStart, int i);
  void subContractionStartRight(lapack_complex_double *&pStart, int i);
  int pCtrLocalIndex(int const ai, int const aip){return ai+D*aip;} 
};

#endif
