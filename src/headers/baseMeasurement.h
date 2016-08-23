#ifndef BASE_MEASUREMENT
#define BASE_MEASUREMENT

#include <complex>
#include "templates/mpo.h"
#include "mps.h"
#include "templates/contractor.h"

//---------------------------------------------------------------------------------------------------//
// The baseMeasurement class is the parent class for the measurements used in the network class. 
// It contains a partial right-side contraction and an iterate function which adds one site to the 
// partial contraction from the left. It can not be directly accessed, only via the child classes.
//---------------------------------------------------------------------------------------------------//

class baseMeasurement{
 protected:
  //Result of iteration is stored in target
  virtual void calcOuterContainerLeft(int const i, std::complex<double> *const source, tmpContainer<std::complex<double> > &outerContainer);
  virtual void calcOuterContainerRight(int i, std::complex<double> *const source, tmpContainer<std::complex<double> > &outerContainer);
  virtual void calcCtrIterLeftBase(int const i, std::complex<double> *const source, std::complex<double> *const targetPctr);
  virtual void calcCtrIterRightBase(int const i, std::complex<double> *const source, std::complex<double> *const target);
  baseMeasurement();
  baseMeasurement(mpo<std::complex<double> > *const MPOperator, mps *const MPState);
  mpo<std::complex<double> > *MPOperator;
  mps *MPState;
  void setupMeasurement(mpo<std::complex<double> > *const MPOperator, mps *const MPState);
  void initializeBase();
  void getLocalDimensions(int i);
  int lDwL, lDwR, lDL, lDR, ld, D, Dw;
  int pctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*D*Dw;}
  int stateIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
  int operatorIndex(int const si, int const sip, int const bi, int const bim) {return bim+bi*Dw+sip*Dw*Dw+si*ld*Dw*Dw;}
 private:
  contractor calcer;
};

#endif
