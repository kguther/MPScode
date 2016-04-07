#ifndef MATRIX_PRODUCT_STATE_VARIABLE_SIZE
#define MATRIX_PRODUCT_STATE_VARIABLE_SIZE

#include "mps.h"
#include "twositeQNOrderMatrix.h"

class imps: public mps{
 public:
  imps(dimensionTable const &dimInfo, std::vector<quantumNumber> const &conservedQNsin);
  void addSite(std::vector<std::complex<int> > const &targetQN);
  void exportState(mps &target);
  twositeQNOrderMatrix centralIndexTable;
};

#endif
