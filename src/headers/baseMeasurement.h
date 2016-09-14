#ifndef BASE_MEASUREMENT
#define BASE_MEASUREMENT

#include "templates/mpo.h"
#include "mps.h"
#include "templates/contractor.h"
#include "mpstype.h"

//---------------------------------------------------------------------------------------------------//
// The baseMeasurement class is the parent class for the measurements used in the network class. 
// It contains a partial right-side contraction and an iterate function which adds one site to the 
// partial contraction from the left. It can not be directly accessed, only via the child classes.
//---------------------------------------------------------------------------------------------------//

class baseMeasurement{
 protected:
  //Result of iteration is stored in target
  virtual void calcOuterContainerLeft(int const i, mpsEntryType *const source, tmpContainer<mpsEntryType > &outerContainer);
  virtual void calcOuterContainerRight(int i, mpsEntryType *const source, tmpContainer<mpsEntryType > &outerContainer);
  virtual void calcCtrIterLeftBase(int const i, mpsEntryType *const source, mpsEntryType *const targetPctr);
  virtual void calcCtrIterRightBase(int const i, mpsEntryType *const source, mpsEntryType *const target);
  baseMeasurement();
  baseMeasurement(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
  mpo<mpsEntryType > *MPOperator;
  mps *MPState;
  void setupMeasurement(mpo<mpsEntryType > *const MPOperator, mps *const MPState);
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
