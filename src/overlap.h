#ifndef MPS_OVERLAP
#define MPS_OVERLAP

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
  overlap(mps const*const psi, mps const*const phi);
  ~overlap();
  overlap& operator=(overlap const &source);
  stateArray const& getF() const{return F;}
  //Updates are done with respect to phi, i.e. psi is expected to be constant during updating
  void loadMPS(mps const*const psi, mps const*const phi);
  arcomplex<double> fullOverlap();
  arcomplex<double> getFullOverlap();
  void stepLeft(int const i);
  void stepRight(int const i);
  arcomplex<double> applyF(arcomplex<double> *vec, int i);
 private:
  int L, D, d;
  mps const *psi;
  mps const *phi;
  arcomplex<double> *Lctr;
  arcomplex<double> *Rctr;
  arcomplex<double>& Lctr_access(int i, int aim, int aimp){return Lctr[aimp+aim*D+i*D*D];}
  arcomplex<double>& Rctr_access(int i, int ai, int aip){return Rctr[aip+ai*D+i*D*D];}
  stateArray F;
  void subContractionStartLeft(arcomplex<double> *&pStart, int i);
  void subContractionStartRight(arcomplex<double> *&pStart, int i);
  void calcCtrIterLeft(int i);
  void calcCtrIterRight(int i);
  void calcCtrIterLeftQNOpt(int i, arcomplex<double> const*const source, arcomplex<double> *const target);
  void calcCtrIterRightQNOpt(int i, arcomplex<double> const*const source, arcomplex<double> *const target);
  void calcF();
  //Important: During sweeping, only the F matrix of the last updated site can be used, and only the mps site matrices of the last updated site should be manipulated. This ensures that the overlap is always up to date. Of course, use the corresponding step for updating (i.e. update the correct direction)
  void updateF(int const i);
  void ovCpy(overlap const &source);
  int pCtrLocalIndex(int aip, int ai){return aip+D*ai;} 
};

#endif
