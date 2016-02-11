#ifndef COLLECTION_OF_RELEVANT_OPERATORS
#define COLLECTION_OF_RELEVANT_OPERATORS

#include "network.h"

//---------------------------------------------------------------------------------------------------//
// Here, operators relevant for our system are stored. The following functions generate 
// important operators like the Hamiltonian or correlations functions in MPO representation.
//---------------------------------------------------------------------------------------------------//

double aMatrix(int const a, int const b);
double bMatrix(int const a, int const b);
int writeHamiltonian(network &sys, double const J, double const g);
int delta(int const a, int const b);


#endif
