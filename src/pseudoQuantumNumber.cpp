#include "pseudoQuantumNumber.h"

pseudoQuantumNumber::pseudoQuantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin):
  dimInfo(dimInfoin),
  N(Nin),
  QNloc(QNlocin)
{}

//---------------------------------------------------------------------------------------------------//

std::complex<int> pseudoQuantumNumber::QNLabel(int si)const {
  return QNloc[si];
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> pseudoQuantumNumber::QNLabel(int i, int ai)const {
  return indexLabel[ai+(i+1)*dimInfo.D()];
}

//---------------------------------------------------------------------------------------------------//

std::complex<int> pseudoQuantumNumber::groupOperation(std::complex<int> const &a, std::complex<int> const &b, int const pre)const {
  //Defines the real part as the U(1) part and the imaginary as the Z_2 part of a quantum number
  std::complex<int> result;
  result.real(real(a)+pre*real(b));
  result.imag(imag(a)*imag(b));
  return result;
}
