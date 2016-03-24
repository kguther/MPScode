#ifndef MPS_OVERLAP
#define MPS_OVERLAP

#include "mkl_complex_defined.h"
#include "mps.h"

//---------------------------------------------------------------------------------------------------//
// The overlap class is the framework for computation of the scalar product of two mps. Note that it
// makes a difference whether one computes the scalar product of an mps with the result of the application
// of some mpo onto that mps or the expectation value of some mpo using the measurement class since
// the latter has more efficient buffering.
//---------------------------------------------------------------------------------------------------//

//WITH THIS CLASS, ONLY MPS WITH THE SAME BOND DIMENSION CAN BE MULTIPLIED. CHILD CLASS FOR CONVERTING MPS TO FIT IN BOND DIMENSION MIGHT BE COMING, BUT IS NOT NECESSARY NOW.

class overlap{
 public:
  overlap();
  overlap(overlap const &source);
  ~overlap();
  overlap& operator=(overlap const &source);
  stateArray F;
  //Updates are done with respect to phi, i.e. psi is expected to be constant during updating
  void loadMPS(mps const*const psi, mps const*const phi);
  lapack_complex_double fullOverlap();
  lapack_complex_double getFullOverlap();
  void stepLeft(int const i);
  void stepRight(int const i);
  lapack_complex_double applyF(lapack_complex_double *vec, int i);
 private:
  int L, D, d;
  mps const *psi;
  mps const *phi;
  lapack_complex_double *Lctr;
  lapack_complex_double *Rctr;
  lapack_complex_double& Lctr_access(int i, int aim, int aimp){return Lctr[aimp+aim*D+i*D*D];}
  lapack_complex_double& Rctr_access(int i, int ai, int aip){return Rctr[aip+ai*D+i*D*D];}
  void subContractionStartLeft(lapack_complex_double *&pStart, int i);
  void subContractionStartRight(lapack_complex_double *&pStart, int i);
  void calcCtrIterLeft(int i);
  void calcCtrIterRight(int i);
  void calcCtrIterLeftQNOpt(int i, lapack_complex_double const*const source, lapack_complex_double *const target);
  void calcCtrIterRightQNOpt(int i, lapack_complex_double const*const source, lapack_complex_double *const target);
  void getF();
  //Important: During sweeping, only the F matrix of the last updated site can be used, and only the mps site matrices of the last updated site should be manipulated. This ensures that the overlap is always up to date. Of course, use the corresponding step for updating (i.e. update the correct direction)
  void updateF(int const i);
  void ovCpy(overlap const &source);
  int pCtrLocalIndex(int aip, int ai){return aip+D*ai;} 
};

#endif
