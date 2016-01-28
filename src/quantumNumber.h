#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

#include <vector>
#include <complex>
#include "dimensionTable.h"

//-------------------------------------------------------------------------------------------//
// The quantumNumber class contains the quantum number  labels for the MPS matrix indices and
// can check whether an entry has to be zero due to QN constraints.
//-------------------------------------------------------------------------------------------//

class quantumNumber{
 public:
  quantumNumber();
  ~quantumNumber();
  void initialize(dimensionTable &dimInfoin, std::complex<int> const Nin, std::complex<int> *QNlocin);
  int qnConstraint(int const i, int const si, int const ai, int const aim);
  std::complex<int> QNLabel(int const i, int const ai);
  std::complex<int> QNLabel(int const si);
  std::complex<int> QNValue() const {return N;}
  int primaryIndex(int const i, int const ai);
  void setParameterD(int const Dnew);
 private:
  dimensionTable dimInfo;
  std::complex<int> *leftLabel;
  std::complex<int> *rightLabel;
  std::complex<int> *indexLabel;
  std::vector<std::vector<int> > primaryIndices;
  std::complex<int> N;
  std::complex<int> *QNloc;
  void initializeLabelList();
  void initializeLabelListLP();
  void initializeLabelListRP();
  int initializeLabelList(int const i, int const direction=-1);
  void gatherBlocks(int const i, std::vector<std::vector<int> > &ai, std::vector<std::complex<int> > &qnLabels, int const direction);
  std::complex<int> exactLabel(int const i, int const ai);
  std::complex<int> truncLabel(int const i, int const ai);
  std::complex<int> QNLabelLP(int const i, int const ai);
  std::complex<int> QNLabelRP(int const i, int const ai);
  std::complex<int> pre(std::complex<int> a, int const direction);
  int integerParity(int const n);
  std::complex<int> groupOperation(std::complex<int> a, std::complex<int> b, int const pre=1);
};

#endif
