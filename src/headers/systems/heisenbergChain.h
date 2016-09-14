#ifndef HEISENBERG_CHAIN
#define HEISENBERG_CHAIN

#include "mpstype.h"
#include "templates/mpo.h"

void generateBoseHubbard(int d, double V, double U, mpo<mpsEntryType > &H);
void generateHubbardHamiltonian(double t, double U, mpo<mpsEntryType > &H);
void generateHeisenbergHamiltonian(mpo<mpsEntryType > &H);
void generateFFHamiltonian(mpo<mpsEntryType > &H);
#endif
