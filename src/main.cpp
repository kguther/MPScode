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
  int const L=15;
  int const nQuantumNumbers=1;
  int hInfo;
  int QNValue[2]={L,1};
  int QNList[8]={0,1,1,2,1,1,-1,1};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(100,4,3,0,1e-4,1e-8,1e-4);
  Qsystem sys(pars,simPars);
  hInfo=writeHamiltonian(&sys,1,1);
  if(hInfo){
    std::cout<<"Invalid bond dimension for the construction of H. Terminating process.\n";
    exit(1);
  }
  /*int lDwR, lDwL, Dw;
  Dw=pars.Dw;
  for(int i=0;i<pars.L;i++){
    if(i==0){
      lDwL=1;
    }
    else{
      lDwL=Dw;
    }
    if(i==(pars.L-1)){
      lDwR=1;
    }
    else{
      lDwR=Dw;
    }
    for(int s=0;s<pars.d.maxd();s++){
      for(int sp=0;sp<pars.d.maxd();sp++){
	for(int bi=0;bi<lDwR;bi++){
	  for(int bim=0;bim<lDwL;bim++){
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 4:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0;
		  break;
		case 0:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=delta(s,sp);
		  break;
		case 1:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=mEl*delta(s,sp+1);
		  break;
		case 2:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=mEl*delta(s,sp-1);
		  break;
		case 3:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=(s-0.5)*delta(s,sp);
		  break;
		default:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0;
		}
	      }
	      else{
	        sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0;
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 1:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0.5*mEl*delta(s,sp-1);
		  break;
		case 2:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0.5*mEl*delta(s,sp+1);
		  break;
		case 3:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=(s-0.5)*delta(s,sp);
		  break;
		case 4:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=delta(s,sp);
		  break;
		default:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0;
		}
	      }
	      else{
	        sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0;
	      }
	    }
	  }
	}
      }
    }
  }
  */
  double matEls;
  mpo<lapack_complex_double> particleNumber(4,2,L);
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
  //Note that the inital state is not normalized, the result of this measurement does not make sense therefore
  sys.measure(particleNumber,spinQN);
  cout<<"Initial total spin: "<<spinQN<<endl;
  sys.TensorNetwork.check=&particleNumber;
  sys.getGroundState();
  cout<<setprecision(21);
  for(int mi=0;mi<nEigens;++mi){
    cout<<"Obtained energy of state "<<mi<<" as: "<<sys.E0[mi]<<endl;
  }
  sys.measure(particleNumber,spinQN);
  cout<<"Final total spin: "<<spinQN<<endl;
}
