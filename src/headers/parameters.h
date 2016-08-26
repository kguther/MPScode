#ifndef PARAMETER_CLASSES
#define PARAMETER_CLASSES

#include <vector>
#include <complex>
#include "localHSpaces.h"

//parameters to be supplied: 
//din: local hilbert space dimension (implicit type conversion from int can be used -> an int can be supplied)
//Lin: system size (number of sites)
//DWin: bond dimension of the hamiltonian
//nEigsin: number of eigenvalues to be computed
//NumberQNs: number of quantum numbers to be fixed (actually, number of pars U(1) + Z_2 QNs)
//IN CASE NO Z_2 SYMMETRY IS USED: JUST TAKE INTs AS QUANTUM NUMBER VALUES (IMPLICIT TYPE CONVERSION)
//QNconservedin: values of the quantum numbers - real part is U(1) part and imaginary part is Z_2 part
//QNListin: quantum numbers of the phyiscal indices
//tReal/tImag: deprecated


class problemParameters{
 public:
  problemParameters(){}
  problemParameters(localHSpaces const &din, int Lin, int Dwin, int nEigsin=1, int NumberQNs=0, std::complex<int> *QNconservedin=0, std::complex<int> * QNListin=0, double tReal=0.0, double tImag=0.0);
  localHSpaces d;
  int L, Dw, nEigs,nQNs;
  std::vector<double> filling;
  std::complex<double> t;
  std::vector<std::complex<int> > QNconserved;
  std::vector<std::vector<std::complex<int> > > QNLocalList;
};

//---------------------------------------------------------------------------------------------------//

//THE CONVERGENCE CHECK MAY OR MAY NOT FAIL FOR EXCITED STATES, SINCE THEY CAN CONTAIN NOTABLE CONTRIBUTIONS FROM THE GROUND STATE AT SOME POINT IN THE ALGORITHM, LEADING TO A FALSE POSITIVE.
//THEREFORE, SET accin=0 (OR SOME VERY LOW VALUE) FOR EXCITED STATE SEARCH
//POSSIBLE FIX: THRESHOLD VALUE FOR OVERLAP WITH GROUNDSTATE -> EXCLUDES STATES TOO SIMILAR TO THE GROUND STATE

//Din: virtual bond dimension of the MPS
//Nin: number of sweeps to be executed
//nStagesin: deprecated, but can be used with the QSystem class
//admixture: initial coefficient for subspace expansion
//accin: threshold value for the energy variance (below this, the run is considered converged and stops)
//tolMinin: minimal tolerance for the eigensolver
//tolInitialin: initial tolerance for the eigensolver

class simulationParameters{
 public:
 simulationParameters(int Din=100, int Nin=4, int nStagesin=2, double admixture=1e-2 ,double accin=1e-8, double tolMinin=1e-8, double tolInitialin=1e-4): D(Din), nSweeps(Nin), nStages(nStagesin),devAccuracy(accin),tolMin(tolMinin),tolInitial(tolInitialin), alpha(admixture){}
  int D, nSweeps, nStages;
  double devAccuracy, tolMin, tolInitial, alpha;
};

#endif
