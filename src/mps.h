#ifndef MATRIX_PRODUCT_STATE
#define MATRIX_PRODUCT_STATE

#include "stateArray.h"
#include "quantumNumber.h"
#include "basisQNOrderMatrix.h"
#include "siteQNOrderMatrix.h"
#include <complex>

class mps: public stateArray{
 public:
  mps();
  mps(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin);
  void generate(dimensionTable const &dimInfoIn, std::vector<quantumNumber> const &conservedQNsin);
  void leftNormalizeState(int i);
  void rightNormalizeState(int i);
  int setParameterD(int Dnew);
  int setParameterL(int Lnew);
  void normalizeFinal(int i);
  void restoreQN(int i);
  void getEntanglementSpectrum(int i, double &S, std::vector<double> &spectrum);
  void getEntanglementEntropy(std::vector<double> &S, std::vector<std::vector<double> > &spectra);
  void adaptLabels(int i, int direction);
  void refineQNLabels(int i, int iQN, std::vector<std::complex<int> > const &source);
  int setUpQNs(std::vector<quantumNumber> const &conservedQNs);
  std::vector<quantumNumber>& getConservedQNs(){return conservedQNs;}
  std::vector<quantumNumber> const& getConservedQNs() const{return conservedQNs;}
  basisQNOrderMatrix const& indexTable()const {return indexTableVar;}
 protected:
  int nQNs;
  std::vector<quantumNumber> conservedQNs;
  void createInitialState();
  void leftNormalizeStateBlockwise(int  i);
  void rightNormalizeStateBlockwise(int i);
  void convertIndicesLP(siteQNOrderMatrix const& localIndexTable, int j, int k, int iBlock, int &si, int &ai, int &aim);
  void convertIndicesRP(siteQNOrderMatrix const& localIndexTable, int j, int k, int iBlock, int &si, int &ai, int &aim);
  int loadIndexTables();
  void loadIndexTablesNoexcept();
  basisQNOrderMatrix indexTableVar;
 private:
  void getEntanglementSpectrumOC(int i, double &S, std::vector<double> &spectrum);
  void leftNormalizePrimitive(int i);
  void rightNormalizePrimitive(int i);
};

void printQNLabels(mps const &test);

#endif
