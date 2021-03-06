#ifndef COLLECTION_OF_RELEVANT_OPERATORS
#define COLLECTION_OF_RELEVANT_OPERATORS

#include "templates/localMpo.h"
#include "network.h"
#include "delta.h"
#include "mpstype.h"
#include <complex>

//---------------------------------------------------------------------------------------------------//
// Here, operators relevant for our system are stored. The following functions generate 
// important operators like the Hamiltonian or correlations functions in MPO representation.
//---------------------------------------------------------------------------------------------------//

double aMatrix(int const a, int const b);
double bMatrix(int const a, int const b);
int writePhasedSecondOrder(localMpo<std::complex<double> > &gamma, double theta);
int writeHamiltonian(network &sys, double J, double g, double W, std::complex<double> t ,double deltaP=0, int tSite=-1);
double disorder(double deltaP);
mpsEntryType tLocalScale(int i);
mpsEntryType tSingleSite(int i, int targetSite);

#endif
