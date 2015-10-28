#include <iostream>
#include <time.h>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "Qsystem.h"

using namespace std;

void testNormalization();
void testLR();
void nancheck(network *n, int i);
void testSolve();
void testMatrix();
int delta(int i, int j);

//-----------------------------------------------------------------//
//HUGE TESTING REALM (CURRENTLY: SINGLE SITE OPTIMIZATION)
//-----------------------------------------------------------------//

int main(int argc, char *argv[]){
  testMatrix();
  return 0;
}


//-------------------------------------------------------------------------------------------//
// DEBUG FUNCTIONS
//-------------------------------------------------------------------------------------------//

void testSolve(){
  double eigVal;
  lapack_complex_double *****array, ****state;
  parameters pars(2,100,10,5,4,1,2);
  Qsystem sys(pars);
  int lD, rD, lDwR, lDwL, Dw;
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
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=delta(s,sp+1);
		  break;
		case 2:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=delta(s,sp-1);
		  break;
		case 3:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=2*(s-0.5)*delta(s,sp);
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
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0.5*delta(s,sp-1);
		  break;
		case 2:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=0.5*delta(s,sp+1);
		  break;
		case 3:
		  sys.TensorNetwork.networkH.global_access(i,s,sp,bi,bim)=2*(s-0.5)*delta(s,sp);
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
  //eigVal=sys.solve();
  sys.getGroundState();
  eigVal=sys.E0;
  cout<<"Obtained grounds state energy as: "<<setprecision(21)<<eigVal<<endl;
}

void testMatrix(){
  lapack_complex_double *****array, ****state, *RTerm, *LTerm, *HTerm;
  lapack_complex_double Ldummy=lapack_make_complex_double(1.0,0.0);
  parameters pars(2,200,20,5,10);
  network system(pars);
  state=system.networkState;
  int lD, rD;
  int icrit=pars.L/2;
  for(int i=0;i<pars.L/2;i++){
    if(pow(pars.d,i+1)>pars.D){
      icrit=i;
      break;
    }
  }
  for(int i=0;i<pars.L;i++){
    for(int s=0;s<pars.d;s++){
      for(int sp=0;sp<pars.d;sp++){
	for(int bi=0;bi<pars.Dw;bi++){
	  for(int bip=0;bip<pars.Dw;bip++){
	    if(bi==bip){
	      system.networkH.global_access(i,s,sp,bi,bip)=1;
	    }
	    else{
	      system.networkH.global_access(i,s,sp,bi,bip)=0;
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
  clock_t curtime;
  curtime=clock();
  system.calcCtrFull(-1);
  system.calcCtrFull(1);
  curtime=clock()-curtime;
  cout<<"Partial contraction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  matrixprint(1,pars.d*pars.d,state[0][0][0]);
  system.Lctr.subContractionStart(10,&LTerm);
  system.Rctr.subContractionStart(10,&RTerm);
  system.networkH.subMatrixStart(10,&HTerm);
  optHMatrix prb(RTerm,LTerm,HTerm,pars,10);
  curtime=clock();
  prb.MultMv(state[10][0][0],state[10][0][0]);
  curtime=clock()-curtime;
  cout<<"Matrix mutliplication took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  matrixprint(1,pars.d*pars.d,state[0][0][0]);
}

int delta(int i, int j){
  if(i==j){
    return 1;
  }
  return 0;
}

/*void testNormalization(){ //CORRECT RESULT: (-1/sqrt(2)) and sqrt(2) as matrix entries in (0) and (1) matrices for leftnormalize and unity for rightnormalize
  lapack_complex_double ****state;
  int dval=3;
  parameters pars(dval,5,8,5,10);
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
  int lDL=system.locDimL(2);
  int lDR=system.locDimR(2);
  for(int i=L-1;i>0;i--){
    system.rightNormalizeState(i);
  }
  system.normalizeFinal(1);
  for(int i=0;i<(L-1);i++){
    system.leftNormalizeState(i);
  }
  system.normalizeFinal(0);
  for(int i=L-1;i>0;i--){
    system.rightNormalizeState(i);
  }
  system.normalizeFinal(1);
  for(int i=0;i<(L-1);i++){
    system.leftNormalizeState(i);
  }
  system.normalizeFinal(0);
  lapack_complex_double *un=new lapack_complex_double[D*D];
  system.leftNormalizationMatrixFull();
  lapack_complex_double zone(1,0);
  for(int i=0;i<D;i++){
    for(int j=0;j<D;j++){
      un[j+i*D]=(i==j)?-1:0;
    }
  }
  for(int i=0;i<L;i++){
    for(int si=0;si<d;si++){
      //matrixprint(system.locDimL(i),system.locDimR(i),state[i][si][0]);
    }
  }
  for(int si=0;si<d;si++){
    transp(lDR,lDL,state[2][si][0]);
  }
  cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,lDR,lDR,d*lDL,&zone,state[2][0][0],lDR,state[2][0][0],lDR,&zone,un,lDR);
  matrixprint(D,D,un);
  cout<<cblas_dznrm2(d*system.locDimL(L-1),system.networkState[L-1][0][0],1)<<endl;
}

void testLR(){ //CORRECT RESULT: 4
  lapack_complex_double *****array, ****state, ****contr;
  parameters pars(2,100,2,5,10);
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
  for(int ai=0;ai<2;ai++){
    for(int bi=0;bi<5;bi++){
      for(int aip=0;aip<2;aip++){
	cout<<contr[pars.L-2][ai][bi][aip]<<endl;
      }
    }
  }
  delete4D(&contr);
}
*/
  
