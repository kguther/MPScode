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
  quantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin);
  int qnConstraint(int i, int si, int ai, int aim);
  std::complex<int> QNLabel(int i, int ai)const;
  std::complex<int> QNLabel(int si)const;
  std::complex<int> QNValue() const {return N;}
  int primaryIndex(int i, int ai);
  int setParameterD(int Dnew);
  int setParameterL(int Lnew);
  std::vector<std::complex<int> > localQNValue() const {return QNloc;}
  //The failed flag is set when the given quantum number cannot be reached in the given system or the initialization failed for some other reason. It is set to zero else.
  int failed;
  std::complex<int> groupOperation(std::complex<int> const &a, std::complex<int> const &b, int pre=1)const;
 private:
  dimensionTable dimInfo;
  std::vector<std::complex<int> > leftLabel;
  std::vector<std::complex<int> > rightLabel;
  std::vector<std::complex<int> > indexLabel;
  std::vector<std::vector<int> > primaryIndices;
  std::complex<int> N;
  std::vector<std::complex<int> > QNloc;
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
  int integerParity(int n);
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
