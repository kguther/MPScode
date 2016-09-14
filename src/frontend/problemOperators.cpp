#include "problemOperators.h"
#include "templates/mpo.h"
#include "delta.h"
#include "math.h"
#include <iostream> 
#include <random>
#include <cstdlib>
#include <ctime>

int writeHamiltonian(network &sys, double J, double g, double W, std::complex<double> t, double deltaP, int tSite){
  std::srand(std::time(0));
  mpo<mpsEntryType > bufH=sys.getNetworkH();
  int const Dw=bufH.maxDim();
  if(Dw!=12){
    //The minimal bond dimension of our Hamiltonian is 12, so the system better has a Hamiltonian of that bond dimension.
    return 2;
  }
  int const L=bufH.length();
  int lDwL, lDwR;
  double prefactor, JD, gD, preD;
#ifdef REAL_MPS_ENTRIES
  mpsEntryType tD=real(t);
  if(std::abs(imag(t))>1e-12){
    std::cout<<"Error: Imaginary hamiltonian matrix elements in run using real mps\n";
    exit(1);
  }
#else
  mpsEntryType tD=t;
#endif
  for(int i=0;i<L;++i){
    lDwL=bufH.locDimL(i);
    lDwR=bufH.locDimR(i);
    if(i==0 || (i==L-1)){
      //The chemical potential term counts each side twice, except for the boundaries, which are counted once.
      prefactor=1;
    }
    else{
      prefactor=2;
    }
    //disorder() is a stochastic function, therefore, it has to be called each time anew
    gD=g*(1+disorder(deltaP));
    if(tSite<0){
      tD=tD*tLocalScale(i);
    }
    else{
      tD=tD*tSingleSite(i,tSite);
    }
    JD=J*(1+disorder(deltaP));
    preD=W*(1+disorder(deltaP));
    for(int si=0;si<sys.locd(i);++si){
      for(int sip=0;sip<sys.locd(i);++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bim=0;bim<lDwL;++bim){
	    //Canoncial form of a local Hamiltonian in MPO representation
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 0:
		  bufH(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		case 1:
		  bufH(i,si,sip,bi,bim)=aMatrix(si,sip);
		  break;
		case 2:
		  bufH(i,si,sip,bi,bim)=bMatrix(si,sip);
		  break;
		case 3:
		  bufH(i,si,sip,bi,bim)=aMatrix(sip,si);
		  break;
		case 4:
		  bufH(i,si,sip,bi,bim)=bMatrix(sip,si);
		  break;
		case 5:
		  bufH(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  bufH(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  bufH(i,si,sip,bi,bim)=delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  bufH(i,si,sip,bi,bim)=delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  bufH(i,si,sip,bi,bim)=delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  bufH(i,si,sip,bi,bim)=delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  bufH(i,si,sip,bi,bim)=prefactor*JD*delta(si,sip)*(1-delta(si,0)+delta(si,3))+tD*delta(sip,2)*delta(si,1)+conj(tD)*delta(sip,1)*delta(si,2);
		  break;
		default:
		  bufH(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  bufH(i,si,sip,bi,bim)=prefactor*JD*delta(si,sip)*(1-delta(si,0)+delta(si,3))+tD*delta(sip,2)*delta(si,1)+conj(tD)*delta(sip,1)*delta(si,2);
		}
		else{
		  bufH(i,si,sip,bi,bim)=0;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 0:
		  bufH(i,si,sip,bi,bim)=prefactor*JD*delta(si,sip)*(1-delta(si,0)+delta(si,3))+tD*delta(sip,2)*delta(si,1)+conj(tD)*delta(sip,1)*delta(si,2);
		  break;		  
		case 1:
		  //THIS IS TRICKY: By construction, a and b anticommute, but only for the same site. For the nearest-neigbhour hopping, one has to take into account extra signs from anticommutation of operators on adjacent sites. All other terms are at least quadratic in the local fermionic operators, so this problem only occurs here. 
		  bufH(i,si,sip,bi,bim)=-aMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
		  break;
		case 2:
		  bufH(i,si,sip,bi,bim)=-bMatrix(sip,si)*(delta(sip,0)-delta(sip,1));
		  break;
		case 3:
		  bufH(i,si,sip,bi,bim)=aMatrix(si,sip)*(delta(sip,3)-delta(sip,1));
		  break;
		case 4:
		  bufH(i,si,sip,bi,bim)=bMatrix(si,sip)*(delta(sip,3)-delta(sip,2));
		  break;
		case 5:
		  bufH(i,si,sip,bi,bim)=-2*JD*delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  bufH(i,si,sip,bi,bim)=-2*JD*delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  bufH(i,si,sip,bi,bim)=gD*delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  bufH(i,si,sip,bi,bim)=gD*delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  bufH(i,si,sip,bi,bim)=-preD*delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  bufH(i,si,sip,bi,bim)=-preD*delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  bufH(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		default:
		  bufH(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		bufH(i,si,sip,bi,bim)=0;
	      }
	    }
	  }
	}
      }
    }
  }
  sys.setNetworkH(bufH);
  return 0;
}

//-------------------------------------------------------------------------------------------//

int writePhasedSecondOrder(localMpo<std::complex<double> > &gamma, double theta){
  std::complex<double> const iUnit=std::complex<double>(0.0,1.0);
  int const Dw=gamma.maxDim();
  if(Dw!=1){
    return 2;
  }
  int const L=gamma.length();
  int lDwL, lDwR;
  for(int i=0;i<L;++i){
    for(int si=0;si<gamma.maxlocd();++si){
      for(int sip=0;sip<gamma.maxlocd();++sip){
	gamma(i,si,sip,0,0)=delta(si,sip);
      }
    }
  }
  for(int si=0;si<gamma.maxlocd();++si){
    for(int sip=0;sip<gamma.maxlocd();++sip){
      gamma(gamma.currentSite(),si,sip,0,0)=(std::exp(iUnit*theta)*delta(si,1)*delta(sip,2)+std::exp(-iUnit*theta)*delta(sip,1)*delta(si,2));
      gamma(gamma.currentSite()-1,si,sip,0,0)=(std::exp(iUnit*theta)*delta(si,1)*delta(sip,2)+std::exp(-iUnit*theta)*delta(sip,1)*delta(si,2));
    }
  }
  return 0;
}

//-------------------------------------------------------------------------------------------//
// For simplicity, these functions return the matrix elements of identity and the fermionic
// site operators a and b, repsectively
//-------------------------------------------------------------------------------------------//


double aMatrix(int const a, int const b){
  if(a==0 && b==1){
    return 1;
  }
  if(a==2 && b==3){
    return 1;
  }
  return 0;
}

//-------------------------------------------------------------------------------------------//

double bMatrix(int const a, int const b){
  if(a==0 && b==2){
    return 1;
  }
  if(a==1 && b==3){
    return -1;
  }
  return 0;
}

//-------------------------------------------------------------------------------------------//

int delta(int const a, int const b){
  if(a==b){
    return 1;
  }
  return 0;
}

//-------------------------------------------------------------------------------------------//

double disorder(double deltaP){
  std::random_device rng;
  double const normalizer=static_cast<double>(rng.max());
  return deltaP*(rng()/normalizer-0.5)*2;
}

//-------------------------------------------------------------------------------------------//

mpsEntryType tLocalScale(int i){
  return 1.0+disorder(0.02);
}

mpsEntryType tSingleSite(int i, int targetSite){
  if(i==targetSite){
    return 1.0;
  }
  return 0.0;
}
