#include "problemOperators.h"
#include "mpo.h"
#include "localMpo.h"
#include "delta.h"
#include "math.h"
#include <iostream> 

int writeHamiltonian(network &sys, double J, double g, double W){
  int const Dw=sys.networkH.maxDim();
  if(Dw!=12){
    //The minimal bond dimension of our Hamiltonian is 12, so the system better has a Hamiltonian of that bond dimension.
    return 2;
  }
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
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		case 1:
		  sys.networkH.global_access(i,si,sip,bi,bim)=aMatrix(si,sip);
		  break;
		case 2:
		  sys.networkH.global_access(i,si,sip,bi,bim)=bMatrix(si,sip);
		  break;
		case 3:
		  sys.networkH.global_access(i,si,sip,bi,bim)=aMatrix(sip,si);
		  break;
		case 4:
		  sys.networkH.global_access(i,si,sip,bi,bim)=bMatrix(sip,si);
		  break;
		case 5:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  sys.networkH.global_access(i,si,sip,bi,bim)=prefactor*J*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		  break;
		default:
		  sys.networkH.global_access(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  sys.networkH.global_access(i,si,sip,bi,bim)=prefactor*J*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		}
		else{
		  sys.networkH.global_access(i,si,sip,bi,bim)=0;
		}
	      }
	    }
	    else{
	      if(bim==lDwL-1){
		switch(bi){
		case 0:
		  sys.networkH.global_access(i,si,sip,bi,bim)=prefactor*J*delta(si,sip)*(1-delta(si,0)+delta(si,3));
		  break;		  
		case 1:
		  //THIS IS TRICKY: By construction, a and b anticommute, but only for the same site. For the nearest-neigbhour hopping, one has to take into account extra signs from anticommutation of operators on adjacent sites. All other terms are at least quadratic in the local fermionic operators, so this problem only occurs here. 
		  sys.networkH.global_access(i,si,sip,bi,bim)=-aMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
		  break;
		case 2:
		  sys.networkH.global_access(i,si,sip,bi,bim)=-bMatrix(sip,si)*(delta(sip,0)-delta(sip,1));
		  break;
		case 3:
		  sys.networkH.global_access(i,si,sip,bi,bim)=aMatrix(si,sip)*(delta(sip,3)-delta(sip,1));
		  break;
		case 4:
		  sys.networkH.global_access(i,si,sip,bi,bim)=bMatrix(si,sip)*(delta(sip,3)-delta(sip,2));
		  break;
		case 5:
		  sys.networkH.global_access(i,si,sip,bi,bim)=-2*J*delta(si,sip)*(delta(si,1)+delta(si,3));
		  break;
		case 6:
		  sys.networkH.global_access(i,si,sip,bi,bim)=-2*J*delta(si,sip)*(delta(si,2)+delta(si,3));
		  break;
		case 7:
		  sys.networkH.global_access(i,si,sip,bi,bim)=g*pre*delta(si,sip)*delta(si,1);
		  break;
		case 8:
		  sys.networkH.global_access(i,si,sip,bi,bim)=g*pre*delta(si,sip)*delta(si,2);
		  break;
		case 9:
		  sys.networkH.global_access(i,si,sip,bi,bim)=-pre*delta(si,1)*delta(sip,2);
		  break;
		case 10:
		  sys.networkH.global_access(i,si,sip,bi,bim)=-pre*delta(si,2)*delta(sip,1);
		  break;
		case 11:
		  sys.networkH.global_access(i,si,sip,bi,bim)=delta(si,sip);
		  break;
		default:
		  sys.networkH.global_access(i,si,sip,bi,bim)=0;
		}
	      }
	      else{
		sys.networkH.global_access(i,si,sip,bi,bim)=0;
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
	gamma.global_access(i,si,sip,0,0)=delta(si,sip);
      }
    }
  }
  for(int si=0;si<gamma.maxlocd();++si){
    for(int sip=0;sip<gamma.maxlocd();++sip){
      gamma.global_access(gamma.currentSite(),si,sip,0,0)=(std::exp(iUnit*theta)*delta(si,1)*delta(sip,2)+std::exp(-iUnit*theta)*delta(sip,1)*delta(si,2));
      gamma.global_access(gamma.currentSite()-1,si,sip,0,0)=(std::exp(iUnit*theta)*delta(si,1)*delta(sip,2)+std::exp(-iUnit*theta)*delta(sip,1)*delta(si,2));
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
