#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include <vector>
#include "stateArray.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"

class mps: public stateArray{
 public:
  mps();
  mps(int const d, int const D, int const L, std::vector<quantumNumber> *conservedQNsin);
  ~mps();
  void generate(int const din, int const Din, int const Lin, std::vector<quantumNumber> *conservedQNsin);
  int leftNormalizeState(int const i);
  int rightNormalizeState(int const i);
  void normalizeFinal(int const i);
  void restoreQN(int const i);
 private:
  int nQNs;
  std::vector<quantumNumber> *conservedQNs;
  std::vector<std::vector<int> > *aiBlockIndicesLP;
  std::vector<std::vector<multInt> > *siaimBlockIndicesLP;
  std::vector<std::vector<int> > *aimBlockIndicesRP;
  std::vector<std::vector<multInt> > *siaiBlockIndicesRP;
  void createInitialState();
  int leftNormalizeStateBlockwise(int const i);
  int rightNormalizeStateBlockwise(int const i);
  void convertIndicesLP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim);
  void convertIndicesRP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim);
};

#endif
