#ifndef BASE_CLASS_FOR_QUANTUM_NUMBERS
#define BASE_CLASS_FOR_QUANTUM_NUMBERS

#include "dimensionTable.h"
#include <vector>
#include <complex>

//---------------------------------------------------------------------------------------------------//
// Base class for anything quantumNumber-like. Contains the core functionality of quantumNumbers: 
// There are label functions and a group operation. How the labeling is constructed is an issue of
// the child classes
//---------------------------------------------------------------------------------------------------//


class pseudoQuantumNumber{
 public:
  pseudoQuantumNumber(){}
  pseudoQuantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin);
  std::complex<int> QNLabel(int i, int ai)const;
  std::complex<int> QNLabel(int si)const;
  std::complex<int> QNValue() const{return N;}
  std::complex<int> groupOperation(std::complex<int> const &a, std::complex<int> const &b, int pre=1)const;
  std::vector<std::complex<int> > localQNValue() const {return QNloc;}
  std::vector<std::complex<int> > const& indexLabelAccess() const{return indexLabel;}
 protected:
  dimensionTable dimInfo;
  std::vector<std::complex<int> > indexLabel;
  std::complex<int> N;
  std::vector<std::complex<int> > QNloc;
};

#endif
