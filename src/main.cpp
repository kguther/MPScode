#include <iostream>
#include <complex>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"

using namespace std;

void testNormalization();
void testLR();
void nancheck(network *n, int i);
void testSolve();
void testMatrix();

//-----------------------------------------------------------------//
//HUGE TESTING REALM (CURRENTLY: SINGLE SITE OPTIMIZATION)
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
    }
  }
  int d=pars.d;
  int D=pars.D;
  int L=pars.L;
  /*for(int i=L-1;i>0;i--){
    system.rightNormalizeState(i);
    }*/
  system.normalizeFinal(1);
  for(int i=0;i<(L-1);i++){
    system.leftNormalizeState(i);
    //nancheck(&system,i);
   }
  system.normalizeFinal(0);
  for(int i=L-1;i>0;i--){
    system.rightNormalizeState(i);
    //nancheck(&system,i);
  }
  system.normalizeFinal(1);
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
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,2,2,2,&zone,state[L-1][0][0],2,state[L-1][0][0],2,&zone,un,2);
  matrixprint(2,2,un);
  cout<<cblas_dznrm2(d*system.locDimR(0),system.networkState[0][0][0],1)<<endl;
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
  system.calcCtrFull(contr,-1);
  cout<<contr[1][0][0][0]<<endl;
  delete4D(&contr);
}

void testSolve(){
  double eigVal;
  lapack_complex_double *****array, ****state, ****contr;
  parameters pars(2,100,100,5,10);
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
  eigVal=system.solve();
  cout<<eigVal<<endl;
}

void testMatrix(){
  lapack_complex_double *****array, ****state, ****Rcontr, ****Lcontr;
  lapack_complex_double Ldummy=lapack_make_complex_double(1.0,0.0);
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
  create4D(pars.L,pars.D,pars.Dw,pars.D,&Lcontr);
  create4D(pars.L,pars.D,pars.Dw,pars.D,&Rcontr);
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
  Lcontr[0][0][0][0]=1;
  system.calcCtrFull(Rcontr,1);
  optHMatrix prb(Rcontr[0][0][0],&Ldummy,system.networkH[0][0][0][0],pars,0);
  prb.MultMv(state[0][0][0],state[0][0][0]);
}
  
