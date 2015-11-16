#ifndef MPS_OVERLAP
#define MPS_OVERLAP

#include <lapacke.h>
#include "mps.h"

//---------------------------------------------------------------------------------------------------//
// The overlap class is the framework for computation of the scalar product of two mps. Note that it
// makes a difference whether one computes the scalar product of an mps with the result of the application
// of some mpo onto that mps or the expectation value of some mpo using the measurement class since
// the latter has more efficient buffering.
//---------------------------------------------------------------------------------------------------//

//FOR NOW, ONLY MPS WITH THE SAME BOND DIMENSION CAN BE MULTIPLIED. EXTENSION MIGHT BE COMING, BUT IS NOT NECESSARY NOW.

class overlap{
 public:
  overlap();
  ~overlap();
  void loadMPS(mps *psi, mps *phi);
 private:
  int L, D, d;
  mps *psi;
  mps *phi;
  lapack_complex_double *Lctr;
  lapack_complex_double *Rctr;
  lapack_complex_double& Lctr_access(int const i, int const aim, int const aimp){return Lctr[aimp+aim*D+i*D*D];}
  lapack_complex_double& Rctr_access(int const i, int const ai, int const aip){return Rctr[aip+ai*D+i*D*D];}
  void subContractionStartLeft(lapack_complex_double *&pStart, int i);
  void subContractionStartRight(lapack_complex_double *&pStart, int i);
  void calcCtrIterLeft(int const i);
  void calcCtrIterRight(int const i);
  int pCtrLocalIndex(int const ai, int const aip){return ai+D*aip;} 
};

#endif
