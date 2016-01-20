#include <iostream>
#include <time.h>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <math.h>
#include "mkl_complex_defined.h"
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "Qsystem.h"
#include "localHSpaces.h"
#include "problemOperators.h"

using namespace std;

void testNormalization();
void testSolve();
void testMatrix();
int delta(int i, int j);

//-----------------------------------------------------------------//
//HUGE TESTING REALM (CURRENTLY: QUANTUM NUMBERS)
//-----------------------------------------------------------------//

int main(int argc, char *argv[]){
  testSolve();
  return 0;
}


//-------------------------------------------------------------------------------------------//
// DEBUG FUNCTIONS
//-------------------------------------------------------------------------------------------//

void testSolve(){
  double eigVal;
  double const mEl=1;
  int const nEigens=1;
  int const L=25;
  int const N=1;
  int const nQuantumNumbers=2;
  int const minimalD=(2*N>4)?2*N:4;
  int hInfo;
  int QNValue[2]={N,N};
  int QNList[8]={0,1,1,2,1,1,-1,-1};
  //Due to poor planning, subchain parity QNs work a bit odd. The QNValue has to be the total particle number, while the parityNumber gives the subchain parity, i.e. it is to be set +-1
  int parityNumber[2]={0,1};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,12,nEigens,nQuantumNumbers,QNValue,QNList,parityNumber);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(minimalD,1,1,0,1e-4,1e-8,1e-4);
  Qsystem sys(pars,simPars);
  hInfo=writeHamiltonian(&sys,1,1);
  if(hInfo){
    std::cout<<"Invalid bond dimension for the construction of H. Terminating process.\n";
    exit(1);
  }
  double matEls;
  mpo<lapack_complex_double> particleNumber(pars.d.maxd(),2,L);
  localMpo<lapack_complex_double> subChainParity(pars.d.maxd(),1,L);
  for(int i=0;i<L;++i){
    for(int bi=0;bi<1;++bi){
      for(int bim=0;bim<1;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip)*(delta(si,1)+delta(si,0)-delta(si,2)-delta(si,3));
	    subChainParity.global_access(i,si,sip,bi,bim)=matEls;
	  }
	}
      }
    }
  }
  for(int i=0;i<L;++i){
    for(int bi=0;bi<2;++bi){
      for(int bim=0;bim<2;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip);
	    if(i!=0 && i!=L-1 && bi==1 && bim==0){
	      matEls=0;
	    }
	    if(bi==0 && bim==particleNumber.locDimL(i)-1){
	      matEls*=(delta(si,1)+delta(si,2)+2*delta(si,3));
	    }
	    particleNumber.global_access(i,si,sip,bi,bim)=matEls;
	  }
	}
      }
    }
  }
  double spinQN;
  sys.TensorNetwork.check=&particleNumber;
  sys.TensorNetwork.checkParity=&subChainParity;
  sys.getGroundState();
  cout<<setprecision(21);
  for(int mi=0;mi<nEigens;++mi){
    cout<<"Obtained energy of state "<<mi<<" as: "<<sys.E0[mi]<<endl;
  }
  sys.measure(particleNumber,spinQN);
  cout<<"Final total particle number: "<<spinQN<<endl;
}
