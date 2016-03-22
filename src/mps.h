#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include <vector>
#include "stateArray.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"

class mps: public stateArray{
 public:
  mps();
  mps(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin);
  mps(mps const &source);
  mps& operator=(mps const &source);
  void generate(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin);
  void mpsCpy(mps const &source);
  int leftNormalizeState(int i);
  int rightNormalizeState(int i);
  int setParameterD(int Dnew);
  void normalizeFinal(int i);
  void restoreQN(int i);
  void getEntanglementSpectrum(int i, double &S, std::vector<double> &spectrum);
  void getEntanglementEntropy(std::vector<double> &S, std::vector<std::vector<double> > &spectra);
  basisQNOrderMatrix indexTable;
 private:
  int nQNs;
  std::vector<quantumNumber> conservedQNs;
  void createInitialState();
  void setUpQNs(std::vector<quantumNumber> const &conservedQNs);
  int leftNormalizeStateBlockwise(int  i);
  int rightNormalizeStateBlockwise(int i);
  void convertIndicesLP(int i, int j, int k, int iBlock, int &si, int &ai, int &aim);
  void convertIndicesRP(int i, int j, int k, int iBlock, int &si, int &ai, int &aim);
  lapack_complex_double exactGroundStateEntry(int i, int si, int ai, int aim);
  void getEntanglementSpectrumOC(int i, double &S, std::vector<double> &spectrum);
};

#endif
