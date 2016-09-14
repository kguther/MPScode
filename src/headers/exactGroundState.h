#ifndef EXACT_POINT_GROUND_STATE
#define EXACT_POINT_GROUND_STATE

#include "mps.h"
#include "quantumNumber.h"

//---------------------------------------------------------------------------------------------------//
// This class can write the ground state at the exact point into any target MPS given, using the
// target's dimension table.
//---------------------------------------------------------------------------------------------------//


class exactGroundState{
 public:
  exactGroundState(std::complex<int> N);
  void writeExactGroundState(mps &target);
 private:
  void generateExactState(mps &target);
  mpsEntryType exactGroundStateEntry(int i, int si, int ai, int aim);
  std::complex<int> QNValue;
  std::vector<quantumNumber> QNsVec;
};

#endif
