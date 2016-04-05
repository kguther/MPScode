#ifndef MATRIX_PRODUCT_STATE_VARIABLE_SIZE
#define MATRIX_PRODUCT_STATE_VARIABLE_SIZE

#include "mps.h"

class imps: public mps{
 public:
  imps();
  imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin);
  void addSite(std::complex<int> *targetQN);
  void exportState(mps &target);
};

#endif
