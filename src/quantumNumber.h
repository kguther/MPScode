#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

class quantumNumber{
 public:
  quantumNumber();
  void initialize(int const din, int const Lin, int const Nin, int *QNlocin);
  int QNLabel(int const i, int const ai, int const lDR);
  int QNLabel(int const si);
 private:
  int L,N,d;
  int *QNloc;
  int QNlocMax;
};

#endif
