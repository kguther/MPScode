#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

#include "dimensionTable.h"

//-------------------------------------------------------------------------------------------//
// The quantumNumber class contains the quantum number  labels for the MPS matrix indices and
// can check whether an entry has to be zero due to QN constraints.
//-------------------------------------------------------------------------------------------//

class quantumNumber{
 public:
  quantumNumber();
  ~quantumNumber();
  void initialize(dimensionTable &dimInfoin, int const Nin, int *QNlocin, int mult=0);
  int qnConstraint(int const i, int const si, int const ai, int const aim);
  int qnCriterium(int const i, int const si, int const ai, int const aim);
  int QNLabel(int const i, int const ai);
  int QNLabel(int const si);
  int QNValue() const {return N;}
  void setParameterD(int const Dnew);
  int parityType() const {return parityNumber;}
 private:
  dimensionTable dimInfo;
  int *leftLabel;
  int N;
  int iLRSwap;
  int QNlocMax;
  int QNlocMin;
  int parityNumber;
  int *QNloc;
  void initializeLabelList();
  int QNLowerCheck(int i, int ai);
  int QNUpperCheck(int i, int ai);
  int exactLabel(int const i, int const ai);
  int truncLabel(int const i, int const ai);
  int groupOperation(int const label, int const labelp);
  int groupInverse(int const label);
  int integerParity(int const n);
};

#endif
