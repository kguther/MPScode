#include <iostream>
#include <complex>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"

using namespace std;

void testNormalization();
void testLR();

//-----------------------------------------------------------------//
//HUGE TESTING REALM (CURRENTLY: NORMALIZATION PROCEDURE)
//-----------------------------------------------------------------//

int main(int argc, char *argv[]){
  testLR();
  return 0;
}


//-------------------------------------------------------------------------------------------//
// DEBUG FUNCTIONS
//-------------------------------------------------------------------------------------------//

void testNormalization(){ //CORRECT RESULT: (-1/sqrt(2)) and sqrt(2) as matrix entries in (0) and (1) matrices for leftnormalize and unity for rightnormalize
  lapack_complex_double ****state;
  parameters pars(2,100,2,5,10);
  int lD,rD;
  network system(pars);
  state=system.networkState;
  for(int i=0;i<pars.L;i++){
    lD=system.locDimL(i);
    rD=system.locDimR(i);
    for(int s=0;s<pars.d;s++){
      for(int ai=0;ai<rD;ai++){
	for(int aim=0;aim<lD;aim++){
	  if(ai==aim){
	    state[i][s][ai][aim]=1;
	  }
	  else{
	    state[i][s][ai][aim]=0;
	  }
	}
      }
      matrixprint(rD,lD,state[i][s][0]);
    }
  }
  cout<<"NORMALIZING\n";
  int d=pars.d;
  int D=pars.D;
  int L=pars.L;
  for(int i=L-1;i>0;i--){
    system.rightNormalizeState(i);
  }
  lapack_complex_double zone(1,0);
  lapack_complex_double *un=new lapack_complex_double[4];
  un[0]=-1;
  un[1]=0;
  un[2]=0;
  un[3]=-1;
  for(int i=0;i<L;i++){
    for(int si=0;si<d;si++){
      matrixprint(system.locDimL(i),system.locDimR(i),state[i][si][0]);
    }
  }
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,2,2,2,&zone,state[1][0][0],2,state[1][0][0],2,&zone,un,2);
  matrixprint(2,2,un);
}

void testLR(){ //CORRECT RESULT: 4
  lapack_complex_double *****array, ****state, ****contr;
  parameters pars(2,100,20,5,10);
  network system(pars);
  array=system.networkH;
  state=system.networkState;
  int lD, rD;
  int icrit=pars.L/2;
  for(int i=0;i<pars.L;i++){
    if(pow(pars.d,i+1)>pars.D){
      icrit=i;
      break;
    }
  }
  create4D(pars.L,pars.D,pars.Dw,pars.D,&contr);
  for(int i=0;i<pars.L;i++){
    for(int s=0;s<pars.d;s++){
      for(int sp=0;sp<pars.d;sp++){
	for(int bi=0;bi<pars.Dw;bi++){
	  for(int bip=0;bip<pars.Dw;bip++){
	    if(bi==bip){
	      array[i][s][sp][bi][bip]=1;
	    }
	    else{
	      array[i][s][sp][bi][bip]=0;
	    }
	  }
	}
      }
    }
  }
  for(int i=0;i<pars.L;i++){
    lD=locDimL(pars.d,pars.D,pars.L,i,icrit);
    rD=locDimR(pars.d,pars.D,pars.L,i,icrit);
    for(int s=0;s<pars.d;s++){
      for(int ai=0;ai<rD;ai++){
	for(int aim=0;aim<lD;aim++){
	  if(ai==aim){
	    state[i][s][ai][aim]=1;
	  }
	  else{
	    state[i][s][ai][aim]=0;
	  }
	}
      }
    }
  }
  system.calcCtrFull(contr,1);
  cout<<contr[pars.L-2][0][0][0]<<endl;
  delete4D(&contr);
  }
