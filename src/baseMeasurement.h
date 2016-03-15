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
 protected:
  //Result of iteration is stored in target
  void calcOuterContainerLeft(int const i, lapack_complex_double *const source, tmpContainer<lapack_complex_double> &outerContainer);
  void calcCtrIterLeftBase(int const i, lapack_complex_double *const source, lapack_complex_double *const targetPctr);
  baseMeasurement();
  baseMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  mpo<lapack_complex_double> *MPOperator;
  mps *MPState;
  void getLocalDimensions(int const i);
  void initializeBase();
  void setupMeasurement(mpo<lapack_complex_double> *const MPOperator, mps *const MPState);
  int lDwL, lDwR, lDL, lDR, ld, D, Dw;
  int pctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*D*Dw;}
  int stateIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
  int operatorIndex(int const si, int const sip, int const bi, int const bim) {return bim+bi*Dw+sip*Dw*Dw+si*ld*Dw*Dw;}
  void calcCtrIterLeftBaseQNOpt(int const i, lapack_complex_double *const source, lapack_complex_double *const targetPctr);
  void calcOuterContainerLeftQNOpt(int const i, lapack_complex_double *const source, tmpContainer<lapack_complex_double> &outerContainer);
};

#endif
