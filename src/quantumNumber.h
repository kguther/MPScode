#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

class quantumNumber{
 public:
  quantumNumber();
  ~quantumNumber();
  void initialize(int const din, int const Lin, int const Nin, int *QNlocin);
  int qnCriterium(int const i, int const si, int const ai, int const aim);
  int QNLabel(int const i, int const ai);
  int QNLabel(int const si);
  int QNLowerCheck(int i, int ai);
  int QNUpperCheck(int i, int ai);
  int *leftLabel;
 private:
  int L,N,d;
  int QNlocMax;
  int QNlocMin;
  int *QNloc;
  int *rightLabel;
  void initializeLabelList();
};

#endif
