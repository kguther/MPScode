#ifndef NETWORK_PROJECTOR
#define NETWORK_PROJECTOR

#include "mkl_complex_defined.h"
#include "mps.h"
#include "overlap.h"
#include "siteArray.h"

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
  projector();
  ~projector();
  void initialize(int const nEigsin);
  void setParameterD(int const Dnew);
  void loadScalarProducts(mps *const variationalState, int const iEigen);
  void updateScalarProducts(int const i, int const direction);
  int getProjector(int const i);
  void project(lapack_complex_double *vec, int const i);
  void storeCurrentState(mps &source);
  void getStoredState(mps *target, int const iEigen);
  int loadNextState(mps &target, int const iEigen);
  int loadNextState(mps &target);
  void storeOrthoState(mps const &source, int const iEigen);
  lapack_complex_double fullOverlap(int const k);
  int nEigen()const {return nCurrentEigen;}
 private:
  mps *orthoStates;
  overlap *scalarProducts;
  int nCurrentEigen;
  siteArray<lapack_complex_double> auxiliaryMatrix;
  void getGramMatrix(lapack_complex_double *gram, int const i);
  void getLocalDimensions(int const i);
  int ld, lDL, lDR;
  int nEigs, nRelevantEigens;
  //lapack_complex_double *projectionMatrix;
  int vecIndex(int si, int ai, int aim){return aim+ai*lDL+si*lDL*lDR;}
};

#endif
