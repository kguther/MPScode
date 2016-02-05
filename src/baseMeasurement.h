#ifndef BASE_MEASUREMENT
#define BASE_MEASUREMENT

#include "mkl_complex_defined.h"
#include "mpo.h"
#include "mps.h"
#include "pContraction.h"
#include "tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// The baseMeasurement class is the parent class for the measurements used in the network class. 
// It contains a partial right-side contraction and an iterate function which adds one site to the 
// partial contraction from the left. It can not be directly accessed, only via the child classes.
//---------------------------------------------------------------------------------------------------//

class baseMeasurement{
 public: 
  //Result of iteration is stored in target
  void calcCtrIterRightBase(int const i, lapack_complex_double *target);
  void calcOuterContainerRight(int const i, tmpContainer<lapack_complex_double> &outerContainer);
  pContraction<lapack_complex_double> Rctr;
 protected:
  baseMeasurement();
  baseMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  mpo<lapack_complex_double> *MPOperator;
  mps *MPState;
  void getLocalDimensions(int const i);
  void initializeBase(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  void setupMeasurement(mpo<lapack_complex_double> *MPOperator, mps *MPState);
  int lDwL, lDwR, lDL, lDR, ld, D, Dw;
  int pctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*D*Dw;}
  int stateIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
  int operatorIndex(int const si, int const sip, int const bi, int const bim) {return bim+bi*Dw+sip*Dw*Dw+si*ld*Dw*Dw;}
  void calcCtrIterRightBaseQNOpt(int const i, lapack_complex_double *target);
  void calcOuterContainerRightQNOpt(int const i, tmpContainer<lapack_complex_double> &outerContainer);
};

#endif
