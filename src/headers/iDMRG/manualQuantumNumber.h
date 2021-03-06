#ifndef MANUAL_QUANTUM_NUMBER_CLASS
#define MANUAL_QUANTUM_NUMBER_CLASS

#include "pseudoQuantumNumber.h"

//---------------------------------------------------------------------------------------------------//
// The explicit version of setting the quantumNumber labeling scheme: This class allows for manual
// setting of the indexLabels using external functions.
//---------------------------------------------------------------------------------------------------//


class manualQuantumNumber: public pseudoQuantumNumber{
 public:
 manualQuantumNumber():pseudoQuantumNumber(){}
 manualQuantumNumber(dimensionTable const &dimInfoin, std::complex<int> const &Nin, std::vector<std::complex<int> > const &QNlocin):
  pseudoQuantumNumber(dimInfoin,Nin,QNlocin)
    {
      indexLabel.resize((1+dimInfo.L())*dimInfo.D());
    }
  void setTargetQN(std::complex<int> const &targetQN){N=targetQN;}
  std::vector<std::complex<int> >& indexLabelAccess(){return indexLabel;}
};

#endif
