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
  std::vector<std::vector<int> > primaryIndices;
  std::complex<int> N;
  std::complex<int> *QNloc;
  void initializeLabelList();
  void initializeLabelList(int const i);
  void gatherBlocks(int const i, std::vector<std::vector<int> > &ai, std::vector<std::vector<int> > &si);
  std::complex<int> exactLabel(int const i, int const ai);
  std::complex<int> truncLabel(int const i, int const ai);
  int integerParity(int const n);
  std::complex<int> groupOperation(std::complex<int> a, std::complex<int> b);
};

#endif
