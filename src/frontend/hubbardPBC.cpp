#include "delta.h"
#include "hubbardPBC.h"
#include "mpstype.h"

double upMatrixA(int a, int b);
double downMatrixA(int a, int b);
double upMatrixB(int a, int b);
double downMatrixB(int a, int b);

void generateHubbardHamiltonianPBC(double t, double U, mpo<mpsEntryType > &H){
  int lDwL, lDwR;
  int const d=8;
  // recall that a PBC system is effectively 2L long
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
		  H(i,si,sip,bi,bim)=upMatrixA(si,sip);
		  break;
		case 2:
		  H(i,si,sip,bi,bim)=downMatrixA(si,sip);
		  break;
		case 3:
		  H(i,si,sip,bi,bim)=upMatrixA(sip,si);
		  break;
		case 4:
		  H(i,si,sip,bi,bim)=downMatrixA(sip,si);
		  break;
		case 5:
		  H(i,si,sip,bi,bim)=upMatrixB(si,sip);
		  break;
		case 6:
		  H(i,si,sip,bi,bim)=downMatrixB(si,sip);
		  break;
		case 7:
		  H(i,si,sip,bi,bim)=upMatrixB(sip,si);
		  break;
		case 8:
		  H(i,si,sip,bi,bim)=downMatrixB(sip,si);
		  break;
		case 9:
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*(delta(si,3)+delta(si,6)+delta(si,9)+delta(si,12)+delta(si,13)+delta(si,14)+2*delta(si,15));
		  break; 
		default:
		  H(i,si,sip,bi,bim)=0.0;
		}
	      }
	      else{
		if(bim==lDwL-1){
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*(delta(si,3)+delta(si,6)+delta(si,9)+delta(si,12)+delta(si,13)+delta(si,14)+2*delta(si,15));
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
		  H(i,si,sip,bi,bim)=U*delta(si,sip)*(delta(si,3)+delta(si,6)+delta(si,9)+delta(si,12)+delta(si,13)+delta(si,14)+2*delta(si,15));
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

double upMatrixA(int a, int b){
}

double downMatrixA(int a, int b){
}

double upMatrixB(int a, int b){
}

double downMatrixB(int a, int b){
}
