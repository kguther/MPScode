#include "delta.h"
#include "heisenbergChain.h"
#include <complex>
#include <cmath>

double upMatrix(int a, int b);
double downMatrix(int a, int b);
double bosonMatrix(int a, int b);

//bose-Hubbard with NN-Interaction

void generateBoseHubbard(int d, double V, double U, mpo<std::complex<double> > &H){
  int lDwL, lDwR;
  int const L=H.length();
  for(int i=0;i<L;++i){
    lDwL=H.locDimL(i);
    lDwR=H.locDimR(i);
    for(int si=0;si<d;++si){
      for(int sip=0;sip<d;++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bim=0;bim<lDwL;++bim){
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 0:
		  H(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		case 1:
		  H(i,si,sip,bi,bim)=bosonMatrix(si,sip);
		  break;
		case 2:
		  H(i,si,sip,bi,bim)=bosonMatrix(sip,si);
		  break;
		case 3:
		  H(i,si,sip,bi,bim)=delta(si,sip)*si;
		  break;
		case 4:
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*si*(si-1);
		  break; 
		default:
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*si*(si-1);
		}
		else{
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 0:
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*si*(si-1);
		  break;
		case 1:
		  H(i,si,sip,bi,bim)=-bosonMatrix(sip,si);
		  break;
		case 2:
		  H(i,si,sip,bi,bim)=-bosonMatrix(si,sip);
		  break;
		case 3:
		  H(i,si,sip,bi,bim)=V*delta(si,sip)*si;
		  break;
		case 4:
		  H(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		default:
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	      else{
		H(i,si,sip,bi,bim)=0.0;
	      }
	    }
	  }
	}
      }
    }
  }
}

void generateHubbardHamiltonian(double t, double U, mpo<std::complex<double> > &H){
  int lDwL, lDwR;
  int const d=4;
  int const L=H.length();
  for(int i=0;i<L;++i){
    lDwL=H.locDimL(i);
    lDwR=H.locDimR(i);
    for(int si=0;si<d;++si){
      for(int sip=0;sip<d;++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bim=0;bim<lDwL;++bim){
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 0:
		  H(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		case 1:
		  H(i,si,sip,bi,bim)=upMatrix(si,sip);
		  break;
		case 2:
		  H(i,si,sip,bi,bim)=downMatrix(si,sip);
		  break;
		case 3:
		  H(i,si,sip,bi,bim)=upMatrix(sip,si);
		  break;
		case 4:
		  H(i,si,sip,bi,bim)=downMatrix(sip,si);
		  break;
		case 5:
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*delta(si,3);
		  break; 
		default:
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*delta(si,3);
		}
		else{
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 0:
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*delta(si,3);
		  break;
		case 1:
		  H(i,si,sip,bi,bim)=t*upMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
		  break;
		case 2:
		  H(i,si,sip,bi,bim)=t*downMatrix(sip,si)*(delta(sip,0)-delta(sip,1));
		  break;
		case 3:
		  H(i,si,sip,bi,bim)=-t*upMatrix(si,sip)*(delta(sip,3)-delta(sip,1));
		  break;
		case 4:
		  H(i,si,sip,bi,bim)=-t*downMatrix(si,sip)*(delta(sip,3)-delta(sip,2));
		  break;
		case 5:
		  H(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		default:
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	      else{
		H(i,si,sip,bi,bim)=0.0;
	      }
	    }
	  }
	}
      }
    }
  }
}


void generateHeisenbergHamiltonian(mpo<std::complex<double> > &H){
  int lDwR, lDwL;
  int const Dw=5;
  int const d=2;
  int const L=H.length();
  for(int i=0;i<L;i++){
    if(i==0){
      lDwL=1;
    }
    else{
      lDwL=Dw;
    }
    if(i==(L-1)){
      lDwR=1;
    }
    else{
      lDwR=Dw;
    }
    double mEl=1.0;
    for(int s=0;s<d;s++){
      for(int sp=0;sp<d;sp++){
	for(int bi=0;bi<lDwR;bi++){
	  for(int bim=0;bim<lDwL;bim++){
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 4:
		  H.global_access(i,s,sp,bi,bim)=0;
		  break;
		case 0:
		  H.global_access(i,s,sp,bi,bim)=delta(s,sp);
		  break;
		case 1:
		  H.global_access(i,s,sp,bi,bim)=mEl*delta(s,sp+1);
		  break;
		case 2:
		  H.global_access(i,s,sp,bi,bim)=mEl*delta(s,sp-1);
		  break;
		case 3:
		  H.global_access(i,s,sp,bi,bim)=(s-0.5)*delta(s,sp);
		  break;
		default:
		  H.global_access(i,s,sp,bi,bim)=0;
		}
	      }
	      else{
	        H.global_access(i,s,sp,bi,bim)=0;
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 1:
		  H.global_access(i,s,sp,bi,bim)=0.5*mEl*delta(s,sp-1);
		  break;
		case 2:
		  H.global_access(i,s,sp,bi,bim)=0.5*mEl*delta(s,sp+1);
		  break;
		case 3:
		  H.global_access(i,s,sp,bi,bim)=(s-0.5)*delta(s,sp);
		  break;
		case 4:
		  H.global_access(i,s,sp,bi,bim)=delta(s,sp);
		  break;
		default:
		  H.global_access(i,s,sp,bi,bim)=0;
		}
	      }
	      else{
		H.global_access(i,s,sp,bi,bim)=0;
	      }
	    }
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------------------------------------//
// For simplicity, these functions return the matrix elements of identity and the fermionic
// site operators a and b, repsectively (a for up, b for down)
//-------------------------------------------------------------------------------------------//


double upMatrix(int a, int b){
  if(a==0 && b==1){
    return 1.0;
  }
  if(a==2 && b==3){
    return 1.0;
  }
  return 0.0;
}

//-------------------------------------------------------------------------------------------//

double downMatrix(int a, int b){
  if(a==0 && b==2){
    return 1.0;
  }
  if(a==1 && b==3){
    return -1.0;
  }
  return 0.0;
}

//-------------------------------------------------------------------------------------------//

double bosonMatrix(int a, int b){
  if(a==(b-1)){
    return sqrt(b);
  }
  return 0.0;
}
/*
int main(int argc, char *argv[]){
  int const nQuantumNumbers=1;
  int D=40;
  int const L=10;
  int const N=L;
  int const up=5;
  int const nSweeps=10;
  int const d=4;
  int const Dw=5;
  double const U=1;
  double const t=1;
  std::complex<int> QNValue[2]={std::complex<int>(N,0),std::complex<int>(up,0)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(2,1),std::complex<int>(3,1),std::complex<int>(0,1),std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,1)};
  localHSpaces localHilbertSpaceDims(d);
  mpo<arcomplex<double> > Hubbard(d,Dw,L);
  generateBoseHubbard(d,t,U,Hubbard);
  problemParameters pars(localHilbertSpaceDims,L,Dw,1,nQuantumNumbers,QNValue,QNList);
  simulationParameters simPars(D,nSweeps,1,1e-3,1e-7,1e-8,1e-3);
  network sys(pars,simPars);
  sys.setNetworkH(Hubbard);
  mpo<std::complex<double> > particleNumber(d,2,L);
  mpo<std::complex<double> > spin(d,2,L);
  double matEls, spinEls;
  for(int i=0;i<pars.L;++i){
    for(int bi=0;bi<2;++bi){
      for(int bim=0;bim<2;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip);
	    spinEls=delta(si,sip);
	    if(i!=0 && i!=L-1 && bi==1 && bim==0){
	      matEls=0.0;
	      spinEls=0.0;
	    }
	    if(bi==0 && bim==particleNumber.locDimL(i)-1){
	      matEls*=(delta(si,1)+delta(si,2)+2*delta(si,3));
	      spinEls*=(delta(si,2)+delta(si,3));
	    }
	    particleNumber.global_access(i,si,sip,bi,bim)=matEls;
	    spin.global_access(i,si,sip,bi,bim)=spinEls;
	  }
	}
      }
    }
  }
  sys.check=&particleNumber;
  sys.checkParity=&spin;
  std::vector<double> E0,dE;
  sys.solve(E0,dE);
  return 0;
}
*/
/*
int main(int argc, char *argv[]){
  int const nQuantumNumbers=1;
  int D=200;
  int const L=100;
  int const nUp=50;
  int const nSweeps=10;
  int const d=2;
  int const Dw=5;
  double const U=1;
  double const t=1;
  std::complex<int> QNValue[2]={std::complex<int>(nUp,0)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1)};
  localHSpaces localHilbertSpaceDims(d);
  mpo<arcomplex<double> > Heisenberg(d,Dw,L);
  generateHeisenbergHamiltonian(Heisenberg);
  problemParameters pars(localHilbertSpaceDims,L,Dw,1,nQuantumNumbers,QNValue,QNList);
  simulationParameters simPars(D,nSweeps,1,1e-3,1e-7,1e-8,1e-3);
  network sys(pars,simPars);
  sys.setNetworkH(Heisenberg);
  mpo<std::complex<double> > particleNumber(d,2,L);
  mpo<std::complex<double> > spin(d,2,L);
  double matEls, spinEls;
  for(int i=0;i<pars.L;++i){
    for(int bi=0;bi<2;++bi){
      for(int bim=0;bim<2;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip);
	    spinEls=delta(si,sip);
	    if(i!=0 && i!=L-1 && bi==1 && bim==0){
	      matEls=0.0;
	      spinEls=0.0;
	    }
	    if(bi==0 && bim==particleNumber.locDimL(i)-1){
	      matEls*=(delta(si,1)-delta(si,0));
	    }
	    particleNumber.global_access(i,si,sip,bi,bim)=matEls;
	    spin.global_access(i,si,sip,bi,bim)=spinEls;
	  }
	}
      }
    }
  }
  sys.check=&particleNumber;
  sys.checkParity=&spin;
  std::vector<double> E0,dE;
  sys.solve(E0,dE);
  return 0;
}
*/
