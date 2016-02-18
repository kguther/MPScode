#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

#include <vector>
#include <complex>
#include "dimensionTable.h"

//-------------------------------------------------------------------------------------------//
// The quantumNumber class contains the quantum number  labels for the MPS matrix indices and
// can check whether an entry has to be zero due to QN constraints.
// We use std::complex<int> as quantum numbers, where the real part can be used as a U(1)-QN
// and the imaginary part as a Z_2-QN (if either is not desired, assign the physical indices
// the neutral element). This makes it easier to take interdependencies of number of particles
// and fermionic parity into account.
// IMPORTANT: THE QN CONSTRAINT IS SUFFICIENT, BUT NOT NECESSARY TO HAVE A GOOD QN. THIS IS
// WHY A QN CONSERVING ALGORITHM IS A BIT TRICKY.
//-------------------------------------------------------------------------------------------//

class quantumNumber{
 public:
  quantumNumber();
  void initialize(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::complex<int> const *const QNlocin);
  int qnConstraint(int i, int si, int ai, int aim);
  std::complex<int> QNLabel(int i, int ai);
  std::complex<int> QNLabel(int si);
  std::complex<int> QNValue() const {return N;}
  int primaryIndex(int i, int ai);
  void setParameterD(int Dnew);
 private:
  dimensionTable dimInfo;
  std::vector<std::complex<int> > leftLabel;
  std::vector<std::complex<int> > rightLabel;
  std::vector<std::complex<int> > indexLabel;
  std::vector<std::vector<int> > primaryIndices;
  std::complex<int> N;
  std::vector<std::complex<int> > QNloc;
  void initializeLabelList();
  void initializeLabelListLP();
  void initializeLabelListRP();
  int initializeLabelList(int i, int direction=-1);
  void gatherBlocks(int i, std::vector<int> &ai, std::vector<std::complex<int> > &qnLabels, int direction);
  std::complex<int> exactLabel(int i, int ai);
  std::complex<int> truncLabel(int i, int ai);
  std::complex<int> QNLabelLP(int i, int ai);
  std::complex<int> QNLabelRP(int i, int ai);
  std::complex<int> pre(std::complex<int> a, int direction);
  int integerParity(int n);
  std::complex<int> groupOperation(std::complex<int> const &a, std::complex<int> const &b, int pre=1);
};

#endif
