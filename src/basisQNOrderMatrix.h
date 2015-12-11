#ifndef MATRIX_TO_CREATE_QN_BLOCK_ORDERING
#define MATRIX_TO_CREATE_QN_BLOCK_ORDERING

#include "mkl_complex_defined.h"
#include "quantumNumber.h"
#include <vector>

struct multInt{
  int aim;
  int si;
};

class basisQNOrderMatrix{
 public:
  basisQNOrderMatrix(dimensionTable &dimin, std::vector<quantumNumber> *conservedQNsin);
  int blockStructure(int const i, std::vector<std::vector<int> > &aiIndices, std::vector<std::vector<multInt> > &siaimIndices);
 private:
  std::vector<quantumNumber> *conservedQNs;
  dimensionTable dimInfo;
};

#endif
