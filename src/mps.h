#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include "stateArray.h"

class mps: public stateArray{
 public:
  mps();
  mps(int const d, int const D, int const L);
  void generate(int const din, int const Din, int const Lin);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
 private:
  int totalQN;
  int conservedQN;
  void createInitialState();
  int QN(int const ai);
};

#endif
