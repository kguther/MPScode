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
  int const L=12;
  int const nQuantumNumbers=1;
  int QNValue[1]={L};
  int QNList[2]={1,-1};
  problemParameters pars(2,L,5,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(100,1,1,1e-3,1e-4,1e-8,1e-4);
  Qsystem sys(pars,simPars);
  int lDwR, lDwL, Dw;
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
    for(int s=0;s<pars.d;s++){
      for(int sp=0;sp<pars.d;sp++){
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
  double matEls;
  mpo<lapack_complex_double> spin(2,2,L);
  for(int i=0;i<L;++i){
    for(int bi=0;bi<2;++bi){
      for(int bim=0;bim<2;++bim){
	for(int si=0;si<pars.d;++si){
	  for(int sip=0;sip<pars.d;++sip){
	    matEls=delta(si,sip);
	    if(i!=0 && i!=L-1 && bi==1 && bim==0){
	      matEls=0;
	    }
	    if(bi==0 && bim==spin.locDimL(i)-1 && si==1){
	      matEls*=-1;
	    }
	    spin.global_access(i,si,sip,bi,bim)=matEls;
	  }
	}
      }
    }
  }
  double spinQN;
  sys.measure(spin,spinQN);
  cout<<"Initial total spin: "<<spinQN<<endl;
  sys.TensorNetwork.check=&spin;
  sys.getGroundState();
  cout<<setprecision(21);
  for(int mi=0;mi<nEigens;++mi){
    cout<<"Obtained energy of state "<<mi<<" as: "<<sys.E0[mi]<<endl;
  }
  sys.measure(spin,spinQN);
  cout<<"Final total spin: "<<spinQN<<endl;
}

int delta(int i, int j){
  if(i==j){
    return 1;
  }
  return 0;
}
