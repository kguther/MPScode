#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include <vector>
#include "stateArray.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"

class mps: public stateArray{
 public:
  mps();
  mps(dimensionTable &dimInfoIn, std::vector<quantumNumber> *conservedQNsin);
  ~mps();
  void generate(dimensionTable &dimInfoIn, std::vector<quantumNumber> *conservedQNsin);
  void setToExactGroundState();
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
  void restoreQN(int const i);
  basisQNOrderMatrix indexTable;
 private:
  int nQNs;
  std::vector<quantumNumber> *conservedQNs;
  void createInitialState();
  int leftNormalizeStateBlockwise(int const i);
  int rightNormalizeStateBlockwise(int const i);
  void convertIndicesLP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim);
  void convertIndicesRP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim);
  lapack_complex_double exactGroundStateEntry(int const i, int const si, int const ai, int const aim);
};

#endif
