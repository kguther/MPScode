#ifndef HEISENBERG_CHAIN
#define HEISENBERG_CHAIN

#include <complex>

void generateBoseHubbard(int d, double V, double U, mpo<std::complex<double> > &H);
void generateHubbardHamiltonian(double t, double U, mpo<std::complex<double> > &H);
void generateHeisenbergHamiltonian(mpo<std::complex<double> > &H);
#endif