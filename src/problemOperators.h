#ifndef COLLECTION_OF_RELEVANT_OPERATORS
#define COLLECTION_OF_RELEVANT_OPERATORS

#include "localMpo.h"
#include "network.h"
#include "delta.h"

//---------------------------------------------------------------------------------------------------//
// Here, operators relevant for our system are stored. The following functions generate 
// important operators like the Hamiltonian or correlations functions in MPO representation.
//---------------------------------------------------------------------------------------------------//

double aMatrix(int const a, int const b);
double bMatrix(int const a, int const b);
int writePhasedSecondOrder(localMpo<std::complex<double> > &gamma, double theta);
int writeHamiltonian(network &sys, double J, double g, double W, double deltaP=0);
int writeHamiltonianSingleParticleHopping(network &sys, double J, double g, double W, std::complex<double> t, double deltaP=0);
double disorder(double deltaP);

#endif
