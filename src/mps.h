#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include "stateArray.h"
#include "quantumNumber.h"

class mps: public stateArray{
 public:
  mps();
  mps(int const d, int const D, int const L, quantumNumber *conservedQNsin);
  void generate(int const din, int const Din, int const Lin, quantumNumber *conservedQNsin);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
 private:
  quantumNumber *conservedQNs;
  void createInitialState();
};

#endif
