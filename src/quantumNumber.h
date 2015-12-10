#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

#include "dimensionTable.h"

class quantumNumber{
 public:
  quantumNumber();
  ~quantumNumber();
  void initialize(dimensionTable &dimInfoin, int const Nin, int *QNlocin);
  int qnCriterium(int const i, int const si, int const ai, int const aim);
  int QNLabel(int const i, int const ai);
  int QNLabel(int const si);
  int QNLowerCheck(int i, int ai);
  int QNUpperCheck(int i, int ai);
  void setParameterD(int const Dnew);
  int *indexLabel;
 private:
  dimensionTable dimInfo;
  int *leftLabel, *rightLabel;
  int N;
  int QNlocMax;
  int QNlocMin;
  int *QNloc;
  void initializeLabelList();
  int exactLabel(int const i, int const ai);
};

#endif
