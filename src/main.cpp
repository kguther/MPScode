#include <iostream>
#include <time.h>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include "mkl_complex_defined.h"
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "Qsystem.h"
#include "localHSpaces.h"
#include "problemOperators.h"
#include "simulation.h"

using namespace std;

void testNormalization();
void sysSolve(int const L, int const N, int const alpha, int const nEigens=1, int const D=1);
void sysSolve(double const J, double const g);
void getScaling(double const J, double const g);

int main(int argc, char *argv[]){
  double J,g;
  if(argc!=3){
    J=1;
    g=0;
  }
  else{
    J=atof(argv[1]);
    g=atof(argv[2]);
  }
  sysSolve(J,g);
  return 0;
}


//-------------------------------------------------------------------------------------------//
// EVALUATION FUNCTIONS
//-------------------------------------------------------------------------------------------//

void sysSolve(double const J, double const g){
  std::string fileName="first/run_1";
  int const nEigens=1;
  int const L=20;
  int const N=20;
  int const D=1;
  int const numPts=1;
  int const nQuantumNumbers=1;
  int const minimalD=(2*N>4)?2*N:4;
  int const usedD=(D>minimalD)?D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(N,-1)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,12,1,1e-3,1e-4,1e-7,1e-4);

  simulation sim(pars,simPars,J,g,numPts,fileName);
  int parityQNs[4]={1,-1,-1,1};

  //The required bond dimension for the perturbed system seems to be greater than that of the unperturbed system
  localMpo<lapack_complex_double> greensFunction(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> densityCorrelation(pars.d.maxd(),1,L,1,0);
  localMpo<lapack_complex_double> localDensity(pars.d.maxd(),1,L,1,0);
  localMpo<lapack_complex_double> interChainCorrelation(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> superconductingOrder(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> interChainDensityCorrelation(pars.d.maxd(),1,L,1,0);
  std::string gFName="Intrachain correlation";
  std::string dCName="Intrachain density correlation";
  std::string lDName="Local density";
  std::string iCDCName="Interchain density correlation";
  std::string iCCName="Interchain correlation";
  std::string scName="Superconducting order parameter";
  //Define some interesting operators in MPO representation. These are mostly correlation functions which are product operators and therefore have Dw=1
  for(int i=0;i<L;++i){
    for(int si=0;si<pars.d.maxd();++si){
      for(int sip=0;sip<pars.d.maxd();++sip){
	greensFunction.global_access(i,si,sip,0,0)=delta(si,sip);
        densityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	localDensity.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	superconductingOrder.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
      }
    }
  }
  for(int si=0;si<pars.d.maxd();++si){
    for(int sip=0;sip<pars.d.maxd();++sip){
      greensFunction.global_access(1,si,sip,0,0)=bMatrix(sip,si);
      greensFunction.global_access(0,si,sip,0,0)=bMatrix(si,sip)*(delta(sip,2)-delta(sip,3));
      densityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      densityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,2)+delta(si,3));
      localDensity.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainCorrelation.global_access(1,si,sip,0,0)=bMatrix(sip,si);
      interChainCorrelation.global_access(0,si,sip,0,0)=aMatrix(si,sip)*(delta(sip,1)-delta(sip,3));
      superconductingOrder.global_access(1,si,sip,0,0)=aMatrix(si,sip);
      superconductingOrder.global_access(0,si,sip,0,0)=aMatrix(si,sip)*(delta(sip,1)-delta(sip,3)); 
    }
  }
  sim.setLocalMeasurement(greensFunction,gFName);
  //The hamiltonian is subchain parity conserving, thus, the expectation vaue of interChainCorrelation is zero
  //sim.setLocalMeasurement(interChainCorrelation,iCCName);
  sim.setLocalMeasurement(densityCorrelation,dCName);
  sim.setLocalMeasurement(interChainDensityCorrelation,iCDCName);
  sim.setLocalMeasurement(localDensity,lDName);
  //The hamiltonian is particle number conserving, thus, the expectation value of superconductingOrder is zero
  //sim.setLocalMeasurement(superconductingOrder,scName);
  clock_t curtime;
  curtime=clock();
  sim.run();
  curtime=clock()-curtime;
  std::cout<<"\nTotal simulation took "<<(float)curtime/CLOCKS_PER_SEC<<" seconds\n\n";
  cout<<setprecision(21);
  for(int mi=0;mi<nEigens;++mi){
    cout<<"Obtained energy of state "<<mi<<" as: "<<sim.E0[mi]<<endl;
  }
}

//-------------------------------------------------------------------------------------------//

void getScaling(int const J, int const g, double const rho, int const odd, int const par){
  std::string fileName="scaling/run_1";
  std::vector<double> energies, accs;
  int const nEigens=2;
  int const L0=30;
  int const LMax=120;
  int N=0;
  int const D=1;
  int const numPts=1;
  int const nQuantumNumbers=1;
  int minimalD=(2*N>4)?2*N:4;
  int usedD=(D>minimalD)?D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(N,-1)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,1,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,12,1,1e-3,1e-4,1e-7,1e-4);
  simulation sim;
  for(int L=L0;L<=LMax;L+=5){
    pars.L=L;
    N=L*rho+odd;
    minimalD=(2*N>4)?2*N:4;
    usedD=(D>minimalD)?D:minimalD;
    simPars.D=usedD;
    QNValue[1]=std::complex<int>(N,par);
    pars.QNconserved=QNValue;
    sim.generate(pars,simPars,J,g,numPts,fileName);
  }
}
