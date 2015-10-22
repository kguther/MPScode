#include <cstdlib>
#include <complex>
#include <iostream>
#include <cblas.h>
#include <arcomp.h>
#include <arscomp.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"

using namespace std;

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//NAMING CONVENTION: i is always the site index. ALWAYS.

//---------------------------------------------------------------------------------------------------//
// 'Main' file containing the MPS ansatz solver itself (this is not the actual main.cpp, just the
// most important file)
//---------------------------------------------------------------------------------------------------//

network::network(parameters inputpars){
  pars=inputpars;                                                                 //TODO: Switch to initialization list (didnt work the first time here, dunno why)
  d=inputpars.d;
  D=inputpars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  N=inputpars.N;
  icrit=L/2;                                                                      //In case chain is too short, this is the correct value (trust me)
  for(int j=0;j<L/2;j++){
    if(pow(d,j+1)>D){
      icrit=j;
      break;
    }
  }
  create5D(L,d,d,Dw,Dw,&networkH);                                                //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied
  createStateArray(d,D,L,&networkState);                                          //Allocation of MPS - memory of the matrices has to be allocated exactly matching the dimension to use them with lapack
}

network::~network(){
  delete5D(&networkH);
  deleteStateArray(&networkState);
}

//---------------------------------------------------------------------------------------------------//
// This is the main solver, using the networkState as initial and final state and networkH as MPO
// representation of the Hamiltonian. Computed are eigenstate (initial state is overwritten) and 
// corresponding eigenvalue (output).
//---------------------------------------------------------------------------------------------------//

double network::solve(){
  lapack_complex_double normalization;
  double lambda;
  int errRet;
  for(int i=L-1;i>0;i--){
    rightNormalizeState(i);
  }
  normalizeFinal(1);
  create4D(L,D,Dw,D,&Lctr);                                                        //Allocates a 4D array of dimension LxDxDwxD which is contigous in the last index
  create4D(L,D,Dw,D,&Rctr);                                                        //for the partial contractions Lctr and Rctr
  Lctr[0][0][0][0]=1;                                                              //initialization neccessary for iterative building of Lctr
  calcCtrFull(Rctr,1);                                                             //Preparing for the first sweep 
  for(int nsweep=0;nsweep<N;nsweep++){
    cout<<"Starting rightsweep\n";
    for(int i=0;i<(L-1);i++){
      errRet=optimize(i,&lambda);                                                  //Eigenvalue solver to actualize cstate
      leftNormalizeState(i);
      calcCtrIterLeft(Lctr,i+1);
    }
    normalizeFinal(0);
    cout<<"Starting leftsweep\n";
    for(int i=L-1;i>0;i--){
      errRet=optimize(i,&lambda);
      rightNormalizeState(i);
      calcCtrIterRight(Rctr,i-1);
    }
    normalizeFinal(1);
  }
  delete4D(&Lctr);                                                                 //Frees the memory previously allocated for partial contractions
  delete4D(&Rctr);
  return lambda;
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int const i, double *iolambda){
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  arcomplex<double> *currentM;
  int nconv;
  optHMatrix HMat(Rctr[i][0][0],Lctr[i][0][0],networkH[i][0][0][0],pars,i);
  plambda=&lambda;
  currentM=networkState[i][0][0];
  ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,&optHMatrix::MultMv,"SR",0,0,1000,currentM);
  nconv=eigProblem.EigenValVectors(currentM,plambda);
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver\n";
    return 1;
  }
  else{
    cout<<"Current energy: "<<real(lambda)<<endl;
    *iolambda=real(lambda);
    return 0;
  }
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state. This is lots of numerical effort and probably the best 
// candidate for optimization
//---------------------------------------------------------------------------------------------------//

void network::calcCtrIterLeft(lapack_complex_double ****Pctr, const int i){ //iteratively builds L expression
  int DR, DL, DwR, DwL;
  int threshold=1e-30;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double ****innercontainer, ****outercontainer;                     //container arrays to significantly reduce computational effort by storing intermediate results
  if(i==1){                                                                         // b_-1  can only take one value
    DwL=1;
  }
  DR=locDimR(i-1);
  DL=locDimL(i-1);
  create4D(d,DwL,DL,DR,&innercontainer);
  create4D(d,DwR,DR,DL,&outercontainer);
  for(int sip=0;sip<d;sip++){                                                       //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
    for(int bim=0;bim<DwL;bim++){
      for(int aim=0;aim<DL;aim++){
	for(int aip=0;aip<DR;aip++){
	  simpleContainer=0;
	  for(int aimp=0;aimp<DL;aimp++){
	    simpleContainer+=Pctr[i-1][aim][bim][aimp]*networkState[i-1][sip][aip][aimp]; 
	  }
	  innercontainer[sip][bim][aim][aip]=simpleContainer;
	}
      }
    }
  }
  //cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<d;si++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	for(int aim=0;aim<DL;aim++){
	  simpleContainer=0;
	  for(int sip=0;sip<d;sip++){
	    for(int bim=0;bim<DwL;bim++){
	      if(abs(networkH[i-1][si][sip][bi][bim])>threshold){               //<- maybe use more sophisticated sparse format
	        simpleContainer+=networkH[i-1][si][sip][bi][bim]*innercontainer[sip][bim][aim][aip];
	      }
	    }
	  }
	  outercontainer[si][bi][aip][aim]=simpleContainer;
	}
      }
    }
  }
  //cout<<"Completed calculation of outer container"<<endl;
  delete4D(&innercontainer);
  for(int ai=0;ai<DR;ai++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	simpleContainer=0;
	for(int si=0;si<d;si++){
	  for(int aim=0;aim<DL;aim++){
	    simpleContainer+=conj(networkState[i-1][si][ai][aim])*outercontainer[si][bi][aip][aim];
	  }
	}
	Pctr[i][ai][bi][aip]=simpleContainer;
      }
    }
  }
  delete4D(&outercontainer);
}

//---------------------------------------------------------------------------------------------------//

void network::calcCtrIterRight(lapack_complex_double ****Pctr, const int i){ //iteratively builds R expression
  int DR, DL, DwR, DwL;
  int threshold=1e-30;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double ****innercontainer, ****outercontainer;
  if(i==L-1){                                                                       // b_L-1  can only take one value
    DwR=1;
  }
  DR=locDimR(i+1);
  DL=locDimL(i+1);
  create4D(d,DwR,DR,DL,&innercontainer);
  create4D(d,DwL,DL,DR,&outercontainer);
  for(int sip=0;sip<d;sip++){                                                       
    for(int bi=0;bi<DwR;bi++){
      for(int ai=0;ai<DR;ai++){
	for(int aimp=0;aimp<DL;aimp++){
	  simpleContainer=0;
	  for(int aip=0;aip<DR;aip++){
	    simpleContainer+=Pctr[i+1][ai][bi][aip]*networkState[i+1][sip][aip][aimp]; 
	  }
	  innercontainer[sip][bi][ai][aimp]=simpleContainer;
	}
      }
    }
  }
  //cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<d;si++){
    for(int bim=0;bim<DwL;bim++){
      for(int aimp=0;aimp<DL;aimp++){
	for(int ai=0;ai<DR;ai++){
	  simpleContainer=0;
	  for(int sip=0;sip<d;sip++){
	    for(int bi=0;bi<DwR;bi++){
	      if(abs(networkH[i+1][si][sip][bi][bim])>threshold){ 
	        simpleContainer+=networkH[i+1][si][sip][bi][bim]*innercontainer[sip][bi][ai][aimp];
	      }
	    }
	  }
	  outercontainer[si][bim][aimp][ai]=simpleContainer;
	}
      }
    }
  }
  //cout<<"Completed calculation of outer container"<<endl;
  delete4D(&innercontainer);
  for(int aim=0;aim<DL;aim++){
    for(int bim=0;bim<DwL;bim++){
      for(int aimp=0;aimp<DL;aimp++){
	simpleContainer=0;
	for(int si=0;si<d;si++){
	  for(int ai=0;ai<DR;ai++){
	    simpleContainer+=conj(networkState[i+1][si][ai][aim])*outercontainer[si][bim][aimp][ai];
	  }
	}
	Pctr[i][aim][bim][aimp]=simpleContainer;
      }
    }
  }
  delete4D(&outercontainer);
}

//---------------------------------------------------------------------------------------------------//
  
int network::calcCtrFull(lapack_complex_double ****Pctr, const int direction){
  int L=pars.L;
  if(direction==1){                                                                //This is just some ordinary iterative computation of the partial contraction Pctr (P=R,L)
    Pctr[L-1][0][0][0]=lapack_make_complex_double(1,0);
    for(int i=L-2;i>=0;i--){
      calcCtrIterRight(Pctr,i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Pctr[0][0][0][0]=lapack_make_complex_double(1,0);
      for(int i=1;i<L;i++){
	calcCtrIterLeft(Pctr,i);
      }
      return 0;
    }
    else{
      cout<<"CRITICAL ERROR: Invalid sweep direction identifier in calculation of partial contractions\n";
      return -1;
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// The following functions are left/right normalizing the matrices of a site after 
// optimization and multiplying the remainder to the matrices of the site to the left/right
// The first two are for truncated sites, the latter two for critical sites, i.e. the site at
// which the truncation is made first. Note that there is no left-normalization of site L-1 
// and no right-normalization of site 1, these are aquired via normalization of the whole state.
//---------------------------------------------------------------------------------------------------//

int network::leftNormalizeState(int const i){
  lapack_int info;
  int D1, D2, D3;                                                                     //Yep, the local dimensions are just named D1, D2, D3 - the decomposition is of a d*D1xD2 matrix
  D1=locDimL(i);
  D2=locDimR(i);
  D3=locDimR(i+1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0.0);
  Qcontainer=new lapack_complex_double[D2];
  Rcontainer=new lapack_complex_double[D2*D2];
  for(int si=0;si<d;si++){
    transp(D2,D1,networkState[i][si][0]);                                              //Enables use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  }
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,d*D1,D2,networkState[i][0][0],D2,Qcontainer);   //Use thin QR decomposition
  upperdiag(D2,D2,networkState[i][0][0],Rcontainer);                               
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,d*D1,D2,D2,networkState[i][0][0],D2,Qcontainer);//Only first D columns are used -> thin QR (below icrit, this is equivalent to a full QR)
  for(int si=0;si<d;si++){
    transp(D1,D2,networkState[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,networkState[i+1][si][0],D2); //REMARK: Use CblasTrans because Rcontainer is in row_major while networkState[i+1][si][0] is in column_major order - this is a normal matrix multiplication - here, R is packed into the matrices of the next site
  }                                                //POSSIBLE TESTS: TEST FOR Q*R - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

int network::rightNormalizeState(int const i){
  lapack_int info;
  int D1,D2,D3;                                                                         //This time, we decompose a D2xd*D1 matrix and multiply the D2xD2 matrix R with a D3xD2 matrix
  D1=locDimR(i);
  D2=locDimL(i);
  D3=locDimL(i-1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0);
  Qcontainer=new lapack_complex_double[d*D1];                                            //Used for storage of lapack-internal matrices
  Rcontainer=new lapack_complex_double[D2*D2];                                           //Used for storage of R from RQ decomposition
  /*cout<<"Normalizing\n";
    matrixprint(D2,d*D1,networkState[i][0][0]);*/
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,d*D1,networkState[i][0][0],D2,Qcontainer);     //Thats how zgerqf works: the last D2 columns contain the upper trigonal matrix R, to adress them, move D2 from the end
  lowerdiag(D2,D2,networkState[i][0][0]+D2*(d*D1-D2),Rcontainer);                        //lowerdiag does get an upper trigonal matrix in column major ordering, dont get confused
  info=LAPACKE_zungrq(LAPACK_COL_MAJOR,D2,d*D1,D2,networkState[i][0][0],D2,Qcontainer);
  /*cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,D2,d*D1,&zone,Rcontainer,D2,networkState[i][0][0],D2);
    matrixprint(D2,d*D1,networkState[i][0][0]);
    cout<<"Testing finished\n";*/
  for(int si=0;si<d;si++){
    cblas_ztrmm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,D3,D2,&zone,Rcontainer,D2,networkState[i-1][si][0],D3);
  }                                                //POSSIBLE TESTS: TEST FOR R*Q - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

void network::normalizeFinal(int const i){
  lapack_complex_double normalization;
  int site, lcD;
  if(i){
    site=0;
    lcD=locDimR(0);
  }
  else{
    site=L-1;
    lcD=locDimL(L-1);
  }
  normalization=cblas_dznrm2(d*lcD,networkState[site][0][0],1);                         //Normalize last matrices to maintain normalization of state
  normalization=1.0/normalization;
  cblas_zscal(d*lcD,&normalization,networkState[site][0][0],1);
}

//---------------------------------------------------------------------------------------------------//
// These functions only return the column/row dimension of the matrices of the i-th site
// Simple, but essential. Using local matrix dimensions make lots of stuff easier.
// NAMING CONVENTION: Any variable-name starting with l indicates a local dimension
// NAMING CONVENTION: Any variable-name of a local dimension ending with R is a local row dimension,
// any ending with L is a local column dimension
//---------------------------------------------------------------------------------------------------//

int network::locDimL(int const i){
  if(i<=icrit){
    return pow(d,i);
  }
  if(i<=L-icrit-1){
    return D;
  }
  return pow(d,L-i);
}

//---------------------------------------------------------------------------------------------------//

int network::locDimR(int const i){
  if(i<icrit){
    return pow(d,i+1);
  }
  if(i<=L-icrit-2){
    return D;
  }
  return pow(d,L-i-1);
}

//---------------------------------------------------------------------------------------------------//
// These functions compute the normalization of the network to the left of some site i. This is for
// testing purpose only since the correct algorithm ensures the results to be trivial
//---------------------------------------------------------------------------------------------------//

void network::leftNormalizationMatrixFull(){
  lapack_complex_double ***psi;
  create3D(L,D,D,&psi);
  for(int i=0;i<L;i++){
    for(int ai=0;ai<D;ai++){
      for(int aip=0;aip<D;aip++){
	psi[i][ai][aip]=0.0/0.0;
      }
    }
  }
  psi[0][0][0]=1;
  cout<<"Printing Psi expressions\n";
  for(int i=1;i<L;i++){
    leftNormalizationMatrixIter(i,psi[0][0]);
    matrixprint(D,D,psi[i][0]);
  }
  cout<<"That's it\n";
  delete3D(&psi);
}

//---------------------------------------------------------------------------------------------------//

void network::leftNormalizationMatrixIter(int i, lapack_complex_double *psi){
  lapack_complex_double psiContainer;
  for(int aim=0;aim<locDimL(i);aim++){
    for(int aimp=0;aimp<locDimL(i);aimp++){
      psiContainer=0;
      for(int si=0;si<d;si++){
	for(int aimm=0;aimm<locDimL(i-1);aimm++){
	  for(int aimmp=0;aimmp<locDimL(i-1);aimmp++){
	    psiContainer+=networkState[i-1][si][aim][aimm]*conj(networkState[i-1][si][aimp][aimm])*psi[(i-1)*D*D+aimm*D+aimmp];
	  }
	}
      }
      psi[i*D*D+aim*D+aimp]=psiContainer;
    }
  }	
}
