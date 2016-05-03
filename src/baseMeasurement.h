#ifndef BASE_MEASUREMENT
#define BASE_MEASUREMENT

#include <arcomp.h>
#include "mpo.h"
#include "mps.h"
#include "contractor.h"

//---------------------------------------------------------------------------------------------------//
// The baseMeasurement class is the parent class for the measurements used in the network class. 
// It contains a partial right-side contraction and an iterate function which adds one site to the 
// partial contraction from the left. It can not be directly accessed, only via the child classes.
//---------------------------------------------------------------------------------------------------//

class baseMeasurement{
 protected:
  //Result of iteration is stored in target
  virtual void calcOuterContainerLeft(int const i, arcomplex<double> *const source, tmpContainer<arcomplex<double> > &outerContainer);
  virtual void calcCtrIterLeftBase(int const i, arcomplex<double> *const source, arcomplex<double> *const targetPctr);
  virtual void calcCtrIterRightBase(int const i, arcomplex<double> *const source, arcomplex<double> *const target);
  baseMeasurement();
  baseMeasurement(mpo<arcomplex<double> > *const MPOperator, mps *const MPState);
  mpo<arcomplex<double> > *MPOperator;
  mps *MPState;
  void setupMeasurement(mpo<arcomplex<double> > *const MPOperator, mps *const MPState);
  void initializeBase();
  void calcOuterContainerRightQNOpt(int const i, arcomplex<double> *const source,  tmpContainer<arcomplex<double> > &outerContainer);
  void getLocalDimensions(int i);
  int lDwL, lDwR, lDL, lDR, ld, D, Dw;
  int pctrIndex(int const ai, int const bi, int const aip) {return aip+bi*D+ai*D*Dw;}
  int stateIndex(int const si, int const ai, int const aim) {return aim+ai*lDL+si*lDL*lDR;}
  int operatorIndex(int const si, int const sip, int const bi, int const bim) {return bim+bi*Dw+sip*Dw*Dw+si*ld*Dw*Dw;}
 private:
  void calcCtrIterLeftBaseQNOpt(int const i, arcomplex<double> *const source, arcomplex<double> *const targetPctr);
  void calcOuterContainerLeftQNOpt(int const i, arcomplex<double> *const source, tmpContainer<arcomplex<double> > &outerContainer);
  void calcCtrIterRightBaseQNOpt(int const i, arcomplex<double> *const source, arcomplex<double> *const target);
  contractor calcer;
};

#endif
