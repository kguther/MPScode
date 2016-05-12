#ifndef QUANTUM_NUMBER_INFO
#define QUANTUM_NUMBER_INFO

#include <vector>
#include <complex>
#include "dimensionTable.h"
#include "pseudoQuantumNumber.h"

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

class quantumNumber: public pseudoQuantumNumber{
 public:
  quantumNumber();
  quantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin);
  int qnConstraint(int i, int si, int ai, int aim);
  int primaryIndex(int i, int ai);
  int setParameterD(int Dnew);
  int setParameterL(int Lnew);
  int grow(int L, int i, std::complex<int> const &targetQN, std::vector<std::complex<int> > const &source);
  //Refine sets the left-side index labels on bond i to source
  int refine(int i, std::vector<std::complex<int> > const &source);
  virtual int validQN(int i, std::complex<int> const &label)const;
  //The failed flag is set when the given quantum number cannot be reached in the given system or the initialization failed for some other reason. It is set to zero else.
  int failed;
 private:
  std::vector<std::complex<int> > leftLabel;
  std::vector<std::complex<int> > rightLabel;
  std::vector<std::vector<int> > primaryIndices;
  int initializeLabelList();
  int initializeLabelListLP();
  int initializeLabelListRP();
  int initializeLabelList(int i, int direction=-1);
  void gatherBlocks(int i, std::vector<int> &ai, std::vector<std::complex<int> > &qnLabels, int direction);
  std::complex<int> exactLabel(int i, int ai);
  std::complex<int> truncLabel(int i, int ai);
  std::complex<int> QNLabelLP(int i, int ai);
  std::complex<int> QNLabelRP(int i, int ai);
  std::complex<int> pre(std::complex<int> a, int direction);
  int integerParity(int n) const;
};

template<typename T>
void reduceMaximum(std::vector<T> &vec, T const &red){
  if(!vec.empty()){
    int pmax=0;
    for(int m=0;m<vec.size();++m){
      if(vec[m]>vec[pmax]){
	pmax=m;
      }
    }
    if(vec[pmax]>red)
      vec[pmax]-=red;
  }
}

template<typename T>
void enforceSum(std::vector<T> &vec, T const &vsum){
  if(!vec.empty()){
    T check=vec[0];
    for(int m=1;m<vec.size();++m){
      check+=vec[m];
    }
    if(check!=vsum){
      reduceMaximum(vec,check-vsum);
    }
  }
}

#endif
