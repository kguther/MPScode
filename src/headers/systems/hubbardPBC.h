#ifndef HUBBARD_PBC_MPO
#define HUBBARD_PBC_MPO

#include "templates/mpo.h"
#include "mpstype.h"

void generateHubbardHamiltonianPBC(double t, double U, mpo<mpsEntryType > &H);

#endif
