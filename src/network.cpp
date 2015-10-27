#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cblas.h>
#include <arcomp.h>
#include <arscomp.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include <time.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "tmpContainer.h"
#include "pContraction.h"
#include "siteoptimizer.h"
#include "mpo.h"

using namespace std;

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//NAMING CONVENTION: i is always the site index. ALWAYS.

//---------------------------------------------------------------------------------------------------//
// 'Main' file containing the MPS ansatz solver itself (this is not the actual main.cpp, just the
// most important file)
//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
// Constructors and similar stuff for the basic network class
//---------------------------------------------------------------------------------------------------//

network::network(){
  //Empty networks dont do much, use the generate function to fill them - this allows for separation of declaration and assignment
  //Set all internal pointers to NULL to allow for destruction of an empty network
  createStateArray(1,1,1,&networkState);
}
network::network(parameters inputpars){
  initialize(inputpars);
}

network::~network(){
  deleteStateArray(&networkState);
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary methods for construction and initialization of network objects
//---------------------------------------------------------------------------------------------------//

void network::generate(parameters inputpars){
  deleteStateArray(&networkState);
  initialize(inputpars);
}

void network::initialize(parameters inputpars){
  int lDR, lDL;
  pars=inputpars;
  d=inputpars.d;
  D=inputpars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  nSweeps=inputpars.nSweeps;
  icrit=L/2;//In case chain is too short, this is the correct value (trust me)
  for(int j=0;j<L/2;j++){
    if(pow(d,j+1)>D){
      icrit=j;
      break;
    }
  }
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access
  networkH.initialize(d,Dw,L);
  //Allocation of MPS - memory of the matrices has to be allocated exactly matching the dimension to use them with lapack
  createStateArray(d,D,L,&networkState);
  for(int i=0;i<pars.L;i++){
    lDL=locDimL(i);
    lDR=locDimR(i);
    for(int s=0;s<pars.d;s++){
      for(int ai=0;ai<lDR;ai++){
	for(int aim=0;aim<lDL;aim++){
	  if(ai==aim){
	    networkState[i][s][ai][aim]=1;
	  }
	  else{
	    networkState[i][s][ai][aim]=0;
	  }
	}
      }
    }
  }
  Lctr.initialize(L,D,Dw);
  Rctr.initialize(L,D,Dw);
}

//---------------------------------------------------------------------------------------------------//
// These functions can be employed to alter the algorithm parameters N and D during lifetime of a 
// network object. This allows for iteratively increasing of D.
//---------------------------------------------------------------------------------------------------//

void network::setParameterNSweeps(int Nnew){
  nSweeps=Nnew;
}

int network::setParameterD(int Dnew){
  if(Dnew<D){
    return -1;
  }
  lapack_complex_double ****newNetworkState;
  createStateArray(d,Dnew,L,&newNetworkState);
  for(int i=0;i<L;i++){
    for(int si=0;si<d;si++){
      for(int ai=0;ai<locDimR(i);ai++){
	for(int aim=0;aim<locDimL(i);aim++){
	  newNetworkState[i][si][ai][aim]=networkState[i][si][ai][aim];
	}
      }
    }
  }
  deleteStateArray(&networkState);
  networkState=newNetworkState;
  D=Dnew;
  pars.D=Dnew;
  for(int j=0;j<L/2;j++){
    if(pow(d,j+1)>D){
      icrit=j;
      break;
    }
  }
  Lctr.initialize(L,D,Dw);
  Rctr.initialize(L,D,Dw);
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// This is the main solver, using the networkState as initial and final state and networkH as MPO
// representation of the Hamiltonian. Computed are eigenstate (initial state is overwritten) and 
// corresponding eigenvalue (output).
//---------------------------------------------------------------------------------------------------//

double network::solve(){  //IMPORTANT TODO: ENHANCE STARTING POINT -> HUGE SPEEDUP
  lapack_complex_double normalization;
  double lambda;
  int errRet;
  clock_t curtime;
  curtime=clock();
  Lctr.initialize(L,D,Dw);
  Rctr.initialize(L,D,Dw);
  for(int i=L-1;i>0;i--){
    rightNormalizeState(i);
  }
  normalizeFinal(1);
  //Allocate a 4D array of dimension LxDxDwxD which is contigous in the last index for the partial contractions Lctr and Rctr
  Lctr.global_access(0,0,0,0)=1;
  calcCtrFull(1);
  for(int iSweep=0;iSweep<nSweeps;iSweep++){
    cout<<"Starting rightsweep\n";
    for(int i=0;i<(L-1);i++){
      cout<<"Optimizing site matrix\n";
      curtime=clock();
      errRet=optimize(i,lambda); 
      curtime=clock()-curtime;
      cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      leftNormalizeState(i);
      calcCtrIterLeft(i+1);
    }
    normalizeFinal(0);
    cout<<"Starting leftsweep\n";
    for(int i=L-1;i>0;i--){
      cout<<"Optimizing site matrix\n";      
      curtime=clock();
      errRet=optimize(i,lambda); 
      curtime=clock()-curtime;
      cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      rightNormalizeState(i);
      calcCtrIterRight(i-1);
    }
    normalizeFinal(1);
  }
  //Free the memory previously allocated for partial contractions
  //delete4D(&Lctr);
  //delete4D(&Rctr);
  cout<<"Determined ground state energy of: "<<lambda<<" using parameters: D="<<D<<" nSweeps="<<nSweeps<<endl;
  return lambda;
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int const i, double &iolambda){
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  arcomplex<double> *currentM;
  arcomplex<double> *RTerm, *LTerm, *HTerm;
  Lctr.subContractionStart(i,&LTerm);
  Rctr.subContractionStart(i,&RTerm);
  networkH.subMatrixStart(i,&HTerm);
  optHMatrix HMat(RTerm,LTerm,HTerm,pars,i);
  plambda=&lambda;
  currentM=networkState[i][0][0];
  ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,&optHMatrix::MultMv,"SR",0,1e-3,400,currentM);
  //One should avoid to hit the maximum number of iterations since this can lead into a suboptimal site matrix, increasing the current energy (although usually not by a lot)
  int nconv;
  nconv=eigProblem.EigenValVectors(currentM,plambda);
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver\n";
    return 1;
  }
  else{
    cout<<setprecision(21)<<"Current energy: "<<real(lambda)<<endl;
    iolambda=real(lambda);
    return 0;
  }
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state. This is lots of numerical effort and probably the best 
// candidate for optimization
//---------------------------------------------------------------------------------------------------//

void network::calcCtrIterLeft(const int i){ //iteratively builds L expression
  int DR, DL, DwR, DwL;
  int threshold=1e-30;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr, *targetPctr;
  Lctr.subContractionStart(i-1,&sourcePctr);
  Lctr.subContractionStart(i,&targetPctr);
  //container arrays to significantly reduce computational effort by storing intermediate results
  if(i==1){
    // b_-1  can only take one value
    DwL=1;
  }
  DR=locDimR(i-1);
  DL=locDimL(i-1);
  tmpContainer<lapack_complex_double> innercontainer(d,DwL,DL,DR);
  tmpContainer<lapack_complex_double> outercontainer(d,DwR,DR,DL);
  //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
  for(int sip=0;sip<d;sip++){
    for(int bim=0;bim<DwL;bim++){
      for(int aim=0;aim<DL;aim++){
	for(int aip=0;aip<DR;aip++){
	  simpleContainer=0;
	  for(int aimp=0;aimp<DL;aimp++){
	    simpleContainer+=sourcePctr[pctrIndex(aim,bim,aimp)]*networkState[i-1][sip][aip][aimp]; 
	  }
	  innercontainer.global_access(sip,bim,aim,aip)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<d;si++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	for(int aim=0;aim<DL;aim++){
	  simpleContainer=0;
	  for(int sip=0;sip<d;sip++){
	    for(int bim=0;bim<DwL;bim++){
	      simpleContainer+=networkH.global_access(i-1,si,sip,bi,bim)*innercontainer.global_access(sip,bim,aim,aip);
	    }
	  }
	  outercontainer.global_access(si,bi,aip,aim)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of outer container"<<endl;
  for(int ai=0;ai<DR;ai++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	simpleContainer=0;
	for(int si=0;si<d;si++){
	  for(int aim=0;aim<DL;aim++){
	    simpleContainer+=conj(networkState[i-1][si][ai][aim])*outercontainer.global_access(si,bi,aip,aim);
	  }
	}
	targetPctr[pctrIndex(ai,bi,aip)]=simpleContainer;
      }
    }
  }
  cout<<"Completed calculation of partial contraction"<<endl;
}

//---------------------------------------------------------------------------------------------------//

void network::calcCtrIterRight(const int i){ //iteratively builds R expression
  int DR, DL, DwR, DwL;
  int threshold=1e-30;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr, *targetPctr;
  Rctr.subContractionStart(i+1,&sourcePctr);
  Rctr.subContractionStart(i,&targetPctr);
  if(i==L-1){
    // b_L-1  can only take one value
    DwR=1;
  }
  DR=locDimR(i+1);
  DL=locDimL(i+1);
  tmpContainer<lapack_complex_double> innercontainer(d,DwR,DR,DL);
  tmpContainer<lapack_complex_double> outercontainer(d,DwL,DL,DR);
  for(int sip=0;sip<d;sip++){                                                       
    for(int bi=0;bi<DwR;bi++){
      for(int ai=0;ai<DR;ai++){
	for(int aimp=0;aimp<DL;aimp++){
	  simpleContainer=0;
	  for(int aip=0;aip<DR;aip++){
	    simpleContainer+=sourcePctr[pctrIndex(ai,bi,aip)]*networkState[i+1][sip][aip][aimp]; 
	  }
	  innercontainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<d;si++){
    for(int bim=0;bim<DwL;bim++){
      for(int aimp=0;aimp<DL;aimp++){
	for(int ai=0;ai<DR;ai++){
	  simpleContainer=0;
	  for(int sip=0;sip<d;sip++){
	    for(int bi=0;bi<DwR;bi++){
	      simpleContainer+=networkH.global_access(i+1,si,sip,bi,bim)*innercontainer.global_access(sip,bi,ai,aimp);
	    }
	  }
	  outercontainer.global_access(si,bim,aimp,ai)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of outer container"<<endl;
  //delete4D(&innercontainer);
  for(int aim=0;aim<DL;aim++){
    for(int bim=0;bim<DwL;bim++){
      for(int aimp=0;aimp<DL;aimp++){
	simpleContainer=0;
	for(int si=0;si<d;si++){
	  for(int ai=0;ai<DR;ai++){
	    simpleContainer+=conj(networkState[i+1][si][ai][aim])*outercontainer.global_access(si,bim,aimp,ai);
	  }
	}
	targetPctr[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
  cout<<"Completed calculation of partial contraction"<<endl;
}

//---------------------------------------------------------------------------------------------------//
  
int network::calcCtrFull(const int direction){
  int L=pars.L;
  //This is just some ordinary iterative computation of the partial contraction Pctr (P=R,L)
  if(direction==1){
    Rctr.global_access(L-1,0,0,0)=lapack_make_complex_double(1,0);
    for(int i=L-2;i>=0;i--){
      calcCtrIterRight(i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Lctr.global_access(0,0,0,0)=lapack_make_complex_double(1,0);
      for(int i=1;i<L;i++){
	calcCtrIterLeft(i);
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
  int D1, D2, D3;
  //Yep, the local dimensions are just named D1, D2, D3 - the decomposition is of a d*D1xD2 matrix
  D1=locDimL(i);
  D2=locDimR(i);
  D3=locDimR(i+1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0.0);
  Qcontainer=new lapack_complex_double[D2];//Used for storage of lapack-internal matrices
  Rcontainer=new lapack_complex_double[D2*D2];//Used for storage of R from RQ decomposition
  //Enable use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  for(int si=0;si<d;si++){
    transp(D2,D1,networkState[i][si][0]);
  }
  //Use thin QR decomposition
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,d*D1,D2,networkState[i][0][0],D2,Qcontainer);
  upperdiag(D2,D2,networkState[i][0][0],Rcontainer);                               
  //Only first D2 columns are used -> thin QR (below icrit, this is equivalent to a full QR)
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,d*D1,D2,D2,networkState[i][0][0],D2,Qcontainer);
  for(int si=0;si<d;si++){
    transp(D1,D2,networkState[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,networkState[i+1][si][0],D2); 
    //REMARK: Use CblasTrans because Rcontainer is in row_major while networkState[i+1][si][0] is in column_major order - this is a normal matrix multiplication - here, R is packed into the matrices of the next site
  }                                                //POSSIBLE TESTS: TEST FOR Q*R - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

int network::rightNormalizeState(int const i){
  lapack_int info;
  //This time, we decompose a D2xd*D1 matrix and multiply the D2xD2 matrix R with a D3xD2 matrix
  int D1,D2,D3;
  D1=locDimR(i);
  D2=locDimL(i);
  D3=locDimL(i-1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0);
  Qcontainer=new lapack_complex_double[d*D1];
  Rcontainer=new lapack_complex_double[D2*D2];
  //Thats how zgerqf works: the last D2 columns contain the upper trigonal matrix R, to adress them, move D2 from the end
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,d*D1,networkState[i][0][0],D2,Qcontainer);
  //lowerdiag does get an upper trigonal matrix in column major ordering, dont get confused
  lowerdiag(D2,D2,networkState[i][0][0]+D2*(d*D1-D2),Rcontainer);
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
  //Normalize last matrices to maintain normalization of state
  normalization=cblas_dznrm2(d*lcD,networkState[site][0][0],1);
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
