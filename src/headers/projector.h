#ifndef NETWORK_PROJECTOR
#define NETWORK_PROJECTOR

#include "mkl_complex_defined.h"
#include "mps.h"
#include "overlap.h"

//---------------------------------------------------------------------------------------------------//
// The projector class contains much more than just the projector onto some excited state space.
// It is used for computation of excited states, therefore, it contains a mps array which can be 
// used to store the computed states. 
//
// Also, it contains an array of overlaps to compute the scalar products of some external state with
// the previously stored ones. From this overlap, the projector onto the space of interest is computed.
//
// The project(..) method can project some on-site matrices onto the space orthogonal to all previously
// stored states. It is automatically used in the matrix class used for energy optimization.
//---------------------------------------------------------------------------------------------------//

class projector{
 public:
  projector(){}
  projector(int nEigsin);
  projector(projector const &source);
  ~projector();
  projector& operator=(projector const &source);
  void setParameterD(int Dnew);
  void setParameterL(int Lnew);
  void loadScalarProducts(mps const &variationalState, int iEigen);
  void updateScalarProducts(int i, int direction);
  int getProjector(int i);
  void project(mpsEntryType *vec, int i);
  void storeCurrentState(mps const &source);
  void getStoredState(mps *&target, int iEigen);
  int loadNextState(mps &target, int iEigen);
  int loadNextState(mps &target);
  void storeOrthoState(mps const &source, int iEigen);
  mpsEntryType fullOverlap(int k);
  int nEigen()const {return nCurrentEigen;}
 private:
  mps *orthoStates;
  overlap *scalarProducts;
  int nCurrentEigen;
  baseTensor<mpsEntryType> auxiliaryMatrix;
  void getGramMatrix(mpsEntryType *gram, int i);
  void getLocalDimensions(int i);
  void pCpy(projector const &source);
  int ld, lDL, lDR;
  int nEigs, nRelevantEigens;
  //mpsEntryType *projectionMatrix;
  int vecIndex(int si, int ai, int aim){return aim+ai*lDL+si*lDL*lDR;}
};

#endif
