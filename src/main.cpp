#include <iostream>
#include <complex>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "parameters.h"
#include "arrayprocessing.h"
#include "arraycreation.h"

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX, i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//----------------------------------------------------------------------------------------------//
// Main file containing the MPS ansatz solver itself
//----------------------------------------------------------------------------------------------//


#ifndef complex_type_plh
#define complex_type_plh lapack_complex_double
#endif

using namespace std;

int solve(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, const parameters pars, double *eigenvalue); //the actual numerics

int leftNormalizeState(int d, int D1, int D2, int D3, int i ,complex_type_plh ****cstate);
int rightNormalizeState(int d, int D1, int D2, int D3, int i, complex_type_plh ****cstate);

int calcCtrFull(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, complex_type_plh ****Pctr, const int direction, const parameters pars);                                                          //computes partial contractions for preparation of a sweep
void calcCtrIter(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, complex_type_plh ****Pctr, const int direction, const int position, const parameters pars); //iteratively builds up the partial contraction during a sweep (for the side not initialized) start at position==(L-1) for R expression and position==0 for L expression

void testNormalization();
void testLR();

//-----------------------------------------------------------------//
//HUGE TESTING REALM (CURRENTLY: GENERATION OF PARTIAL CONTRACTIONS) -- WORKS, BUT SLOW!! are R/L-expressions sparse?
//-----------------------------------------------------------------//

int main(int *argc, char *argv[]){
  testLR();
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// This is the main solver. Input are initial state, hamiltonian in MPO representation and parameter 
// set. Output are eigenstate (initial state is overwritten) and corresponding eigenvalue.
//---------------------------------------------------------------------------------------------------//

int solve(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, const parameters pars, double *eigenvalue){
  complex_type_plh ****Lctr; 
  complex_type_plh ****Rctr;
  complex_type_plh *Qcontainer, *Rcontainer;
  complex_type_plh normalization;
  int D=pars.D;
  int L=pars.L;
  int N=pars.N;
  int Dw=pars.Dw;
  int Dl=D;
  int d=pars.d;
  int icrit=L/2;
  for(int i=0;i<L;i++){                                                            //Undefined if chain is too short
    if(pow(d,i)>D){
      icrit=i-1;
      break;
    }
  }
  create4D(L,D,Dw,D,&Lctr);                                                        //Allocates a 4D array of dimension LxDxDwxD which is contigous in the last index
  create4D(L,D,Dw,D,&Rctr);                                                        //for the partial contractions Lctr and Rctr
  Lctr[0][0][0][0]=1;                                                              //initialization neccessary for iterative building of Lctr
  calcCtrFull(cstate,hamiltonian,Rctr,1,pars);                                     //Preparing for the first sweep 
  for(int nsweep=0;nsweep<N;nsweep++){
    for(int i=0;i<(L-1);i++){
      //Eigenvalue solver to actualize cstate
      leftNormalizeState(d,locDimL(d,D,L,i,icrit),locDimR(d,D,L,i,icrit),locDimR(d,D,L,i+1,icrit),i,cstate);
      calcCtrIter(cstate,hamiltonian,Lctr,-1,i+1,pars);
    }
    for(int i=L;i>0;i--){
      //Eigenvalue solver to actualize cstate
      rightNormalizeState(d,locDimR(d,D,L,i,icrit),locDimL(d,D,L,i,icrit),locDimL(d,D,L,i-1,icrit),i,cstate);
      calcCtrIter(cstate,hamiltonian,Rctr,1,i-1,pars);
    }
  }
  delete4D(&Lctr);                                                               //Frees the memory previously allocated for partial contractions
  delete4D(&Rctr);
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state. This is lots of numerical effort and probably the best 
// candidate for optimization
//---------------------------------------------------------------------------------------------------//

void calcCtrIter(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, complex_type_plh ****Pctr, const int direction, const int position, const parameters pars){ //direction==1 builds R expression, direction==-1 builds L expression
  int D, sitedim, Dw, Dm, Dwm, icrit;
  D=pars.D;
  sitedim=pars.d;
  Dw=pars.Dw;   
  Dm=D;                                                                             //load parameters from input, might be replaced by wrapper class
  Dwm=Dw;
  complex_type_plh ****innercontainer, ****outercontainer;                         //container arrays to significantly reduce computational effort by storing intermediate results
  if((position==1) && (direction==-1)){                                            // b_-1  can only take one value
    Dwm=1;
  }
  if((position==(pars.L-2)) && (direction==1)){                                    //as does b_L-1
    Dw=1;
  }
  icrit=pars.L/2;
  for(int i=0;i<=pars.L;i++){
    if(D<pow(sitedim,i+1)){
      icrit=i;
      break;
    }
    //TODO: Add exception throw
  }
  D=locDimR(sitedim,D,pars.L,position+direction,icrit);
  Dm=locDimL(sitedim,D,pars.L,position+direction,icrit);
  create4D(sitedim,Dwm,Dm,D,&innercontainer);
  create4D(sitedim,Dw,D,Dm,&outercontainer);
  for(int sip=0;sip<sitedim;sip++){                                                //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
    for(int bim=0;bim<Dwm;bim++){
      for(int aim=0;aim<Dm;aim++){
	for(int aip=0;aip<D;aip++){
	  innercontainer[sip][bim][aim][aip]=0;
	  for(int aimp=0;aimp<Dm;aimp++){
	    innercontainer[sip][bim][aim][aip]+=Pctr[position+direction][aim][bim][aimp]*cstate[position+direction][sip][aip][aimp]; 
	  }
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<sitedim;si++){ //calculation of outer container seems to be the computationally most effortive task (fortunately, it is easily optimized)
    for(int bi=0;bi<Dw;bi++){
      for(int aip=0;aip<D;aip++){
	for(int aim=0;aim<Dm;aim++){
	  outercontainer[si][bi][aip][aim]=0;
	  for(int sip=0;sip<sitedim;sip++){
	    for(int bim=0;bim<Dwm;bim++){
	      if(hamiltonian[position+direction][si][sip][bi][bim]!=0){            //Make use of sparse hamiltonian (very helpful) <- maybe use more sophisticated sparse format
		outercontainer[si][bi][aip][aim]+=hamiltonian[position+direction][si][sip][bi][bim]*innercontainer[sip][bim][aim][aip];
	      }
	    }
	  }
	}
      }
    }
  }
  cout<<"Completed calculation of outer container"<<endl;
  delete4D(&innercontainer);
  for(int ai=0;ai<D;ai++){
    for(int bi=0;bi<Dw;bi++){
      for(int aip=0;aip<D;aip++){
	Pctr[position][ai][bi][aip]=0;
	for(int si=0;si<sitedim;si++){
	  for(int aim=0;aim<Dm;aim++){
	    Pctr[position][ai][bi][aip]+=conj(cstate[position+direction][si][ai][aim])*outercontainer[si][bi][aip][aim];
	  }
	}
      }
    }
  }
  delete4D(&outercontainer);
}
  
int calcCtrFull(complex_type_plh ****cstate, complex_type_plh *****hamiltonian, complex_type_plh ****Pctr, const int direction, const parameters pars){
  int L=pars.L;
  if(direction==1){
    Pctr[L-1][0][0][0]=lapack_make_complex_double(1,0);
    for(int i=L-2;i>=0;i--){
      calcCtrIter(cstate,hamiltonian,Pctr,direction,i,pars);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Pctr[0][0][0][0]=lapack_make_complex_double(1,0);
      for(int i=1;i<L;i++){
	calcCtrIter(cstate,hamiltonian,Pctr,direction,i,pars);
      }
      return 0;
    }
    else{
      return -1;
    }
  }
}

//--------------------------------------------------------------------------------------------//
// The following functions are left/right normalizing the matrices of a site after 
// optimization and multiplying the remainder to the matrices of the site to the left/right
// The first two are for truncated sites, the latter two for critical sites, i.e. the site at
// which the truncation is made first. Note that there is no left-normalization of site L-1 
// and no right-normalization of site 1, these are aquired via normalization of the whole state.
//--------------------------------------------------------------------------------------------// 

int leftNormalizeState(int d, int D1, int D2, int D3, int i, complex_type_plh ****cstate){
  lapack_int info;
  complex_type_plh *Rcontainer, *Qcontainer;
  const complex_type_plh zone=lapack_make_complex_double(1.0,0);
  Qcontainer=new complex_type_plh[d*D1];
  Rcontainer=new complex_type_plh[D2*D2];
  for(int si=0;si<d;si++){
    transp(D1,D2,cstate[i][si][0]);                                      //Enables use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  }
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,d*D1,D2,cstate[i][0][0],D2,Qcontainer);  //Use thin QR decomposition
  upperdiag(D2,D2,cstate[i][0][0],Rcontainer);                               
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,d*D1,D2,D2,cstate[i][0][0],D2,Qcontainer);//Only first D columns are used -> thin QR
  for(int si=0;si<d;si++){
    transp(D2,D1,cstate[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,cstate[i+1][si][0],D2); //REMARK: Use CblasTrans because Rcontainer is in row_major while cstate[i+1][si][0] is in column_major order - this is a normal matrix multiplication
  }                                                //POSSIBLE TESTS: TEST FOR Q*R
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

int rightNormalizeState(int d, int D1, int D2, int D3, int i, complex_type_plh ****cstate){
  lapack_int info;
  complex_type_plh *Rcontainer, *Qcontainer;
  const complex_type_plh zone=lapack_make_complex_double(1.0,0);
  Qcontainer=new complex_type_plh[d*D1];
  Rcontainer=new complex_type_plh[D2*D2];
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,d*D1,cstate[i][0][0],D2,Qcontainer);
  lowerdiag(D2,D2,cstate[i][d-1][0],Rcontainer);
  info=LAPACKE_zungrq(LAPACK_COL_MAJOR,D2,d*D1,D2,cstate[i][0][0],D2,Qcontainer);
  for(int si=0;si<d;si++){
    cblas_ztrmm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,D3,D2,&zone,Rcontainer,D2,cstate[i-1][si][0],D2);
  }                                                //POSSIBLE TESTS: TEST FOR R*Q
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//-------------------------------------------------------------------------------------------//
// DEBUG FUNCTIONS
//-------------------------------------------------------------------------------------------//

void testNormalization(){ //CORRECT RESULT: (-1/sqrt(2)) and sqrt(2) as matrix entries in (0) and (1) matrices for leftnormalize and unity for rightnormalize
  complex_type_plh *****array, ****state, ****contr;
  parameters pars(2,100,2,5,10);
  int lD, rD;
  int icrit=pars.L/2;                     
  createStateArray(pars.d,pars.D,pars.L,&state);
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
	matrixprint(rD,lD,state[i][s][0]);
    }
  }
  cout<<"NORMALIZING\n";
  int d=pars.d;
  int D=pars.D;
  int L=pars.L;
  for(int i=L-1;i>0;i--){
      rightNormalizeState(d,locDimR(d,D,L,i,icrit),locDimL(d,D,L,i,icrit),locDimL(d,D,L,i-1,icrit),i,state);
  }
  complex_type_plh zone(1,0);
  complex_type_plh *un=new complex_type_plh[4];
  un[0]=-1;
  un[1]=0;
  un[2]=0;
  un[3]=-1;
  for(int i=0;i<L;i++){
    for(int si=0;si<d;si++){
    matrixprint(locDimL(d,D,L,i,icrit),locDimR(d,D,L,i,icrit),state[i][si][0]);
    }
  }
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,2,2,2,&zone,state[1][0][0],2,state[1][0][0],2,&zone,un,2);
  matrixprint(2,2,un);
  deleteStateArray(&state);
}

void testLR(){ //CORRECT RESULT: 4
  complex_type_plh *****array, ****state, ****contr;
  parameters pars(2,100,20,5,10);
  int lD, rD;
  int icrit=pars.L/2;
  for(int i=0;i<pars.L;i++){
    if(pow(pars.d,i+1)>pars.D){
      icrit=i;
      break;
    }
  }
  create5D(pars.L,pars.d,pars.d,pars.Dw,pars.Dw,&array);
  create4D(pars.L,pars.D,pars.Dw,pars.D,&contr);
  createStateArray(pars.d,pars.D,pars.L,&state);
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
  calcCtrFull(state,array,contr,1,pars);
  cout<<contr[pars.L-2][0][0][0]<<endl;
  delete4D(&contr);
  deleteStateArray(&state);
  delete5D(&array);
}
