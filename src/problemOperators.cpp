#include "problemOperators.h"
#include "mpo.h"
#include "localMpo.h"
#include "delta.h"
#include "math.h"
#include <iostream> 
#include <cstdlib>
#include <ctime>

int writeHamiltonian(network &sys, double J, double g, double W, double deltaP){
  std::srand(std::time(0));
  int const Dw=sys.networkH.maxDim();
  if(Dw!=12 && Dw!=14){
    //The minimal bond dimension of our Hamiltonian is 12, so the system better has a Hamiltonian of that bond dimension.
    return 2;
  }
  std::cout<<deltaP<<std::endl;
  int const L=sys.networkH.length();
  int lDwL, lDwR;
  double prefactor;
  //The prefactor pre exists only for consistency checks and is the relative weight of the inter- and intrachain term
  double pre=W;
  for(int i=0;i<L;++i){
    lDwL=sys.networkH.locDimL(i);
    lDwR=sys.networkH.locDimR(i);
    if(i==0 || (i==L-1)){
      //The chemical potential term counts each side twice, except for the boundaries, which are counted once.
      prefactor=1;
    }
    else{
      prefactor=2;
    }
    for(int si=0;si<sys.locd(i);++si){
      for(int sip=0;sip<sys.locd(i);++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bim=0;bim<lDwL;++bim){
	    //Canoncial form of a local Hamiltonian in MPO representation
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 0:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		case 1:
		  sys.networkH(i,si,sip,bi,bim)=aMatrix(si,sip);
		  break;
		case 2:
		  sys.networkH(i,si,sip,bi,bim)=bMatrix(si,sip);
		  break;
		case 3:
		  sys.networkH(i,si,sip,bi,bim)=aMatrix(sip,si);
		  break;
		case 4:
		  sys.networkH(i,si,sip,bi,bim)=bMatrix(sip,si);
		  break;
		case 5:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  sys.networkH(i,si,sip,bi,bim)=prefactor*J*(1+disorder(deltaP))*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		  break;
		default:
		  sys.networkH(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  sys.networkH(i,si,sip,bi,bim)=prefactor*J*(1+disorder(deltaP))*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		}
		else{
		  sys.networkH(i,si,sip,bi,bim)=0;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 0:
		  sys.networkH(i,si,sip,bi,bim)=prefactor*J*(1+disorder(deltaP))*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		  break;		  
		case 1:
		  //THIS IS TRICKY: By construction, a and b anticommute, but only for the same site. For the nearest-neigbhour hopping, one has to take into account extra signs from anticommutation of operators on adjacent sites. All other terms are at least quadratic in the local fermionic operators, so this problem only occurs here. 
		  sys.networkH(i,si,sip,bi,bim)=-aMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
		  break;
		case 2:
		  sys.networkH(i,si,sip,bi,bim)=-bMatrix(sip,si)*(delta(sip,0)-delta(sip,1));
		  break;
		case 3:
		  sys.networkH(i,si,sip,bi,bim)=aMatrix(si,sip)*(delta(sip,3)-delta(sip,1));
		  break;
		case 4:
		  sys.networkH(i,si,sip,bi,bim)=bMatrix(si,sip)*(delta(sip,3)-delta(sip,2));
		  break;
		case 5:
		  sys.networkH(i,si,sip,bi,bim)=-2*J*(1+disorder(deltaP))*delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  sys.networkH(i,si,sip,bi,bim)=-2*J*(1+disorder(deltaP))*delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  sys.networkH(i,si,sip,bi,bim)=g*(1+disorder(deltaP))*pre*delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  sys.networkH(i,si,sip,bi,bim)=g*(1+disorder(deltaP))*pre*delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  sys.networkH(i,si,sip,bi,bim)=-pre*(1+disorder(deltaP))*delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  sys.networkH(i,si,sip,bi,bim)=-pre*(1+disorder(deltaP))*delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  sys.networkH(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		default:
		  sys.networkH(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		sys.networkH(i,si,sip,bi,bim)=0;
	      }
	    }
	  }
	}
      }
    }
  }
  return 0;
}

//-------------------------------------------------------------------------------------------//

int writeHamiltonianSingleParticleHopping(network &sys, double J, double g, double W, std::complex<double> t, double deltaP){
  int info;
  int const Dw=sys.networkH.maxDim();
  if(Dw!=14){
    return 2;
  }
  info=writeHamiltonian(sys,J,g,W,deltaP);
  int const L=sys.networkH.length();
  int lDwL, lDwR;
  for(int i=0;i<L;++i){
    lDwL=sys.networkH.locDimL(i);
    lDwR=sys.networkH.locDimR(i);
    for(int si=0;si<sys.locd(i);++si){
      for(int sip=0;sip<sys.locd(i);++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bim=0;bim<lDwL;++bim){
	    if(bi==0){
	      if(i!=0){
		switch(bim){
		case 12:
		  sys.networkH(i,si,sip,bi,bim)=aMatrix(sip,si)*(delta(si,0)-delta(si,2));
		  break;
		case 13:
		  sys.networkH(i,si,sip,bi,bim)=bMatrix(sip,si)*(delta(si,0)-delta(si,1));
		  break;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 12:
		  sys.networkH(i,si,sip,bi,bim)=t*bMatrix(si,sip);
		  break;
		case 13:
		  sys.networkH(i,si,sip,bi,bim)=t*aMatrix(si,sip);
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return info;
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
  double dval=static_cast<double>(rand())/RAND_MAX;
  return deltaP*dval;
}
