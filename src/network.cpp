#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <arcomp.h>
#include <arscomp.h>
#include <time.h>
#include <lapacke.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "tmpContainer.h"
#include "pContraction.h"
#include "siteoptimizer.h"
#include "mpo.h"
#include "mps.h"

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
}

//---------------------------------------------------------------------------------------------------//

network::network(parameters inputpars){
  initialize(inputpars);
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary methods for construction and initialization of network objects
//---------------------------------------------------------------------------------------------------//

void network::initialize(parameters inputpars){
  pars=inputpars;
  d=inputpars.d;
  D=inputpars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  nSweeps=inputpars.nSweeps;
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access
  networkH.initialize(d,Dw,L);
  //Allocation of MPS - memory of the matrices has to be allocated exactly matching the dimension to use them with lapack and cblas
  networkState.generate(d,D,L);
  Lctr.initialize(L,D,Dw);
  Rctr.initialize(L,D,Dw);
}

//---------------------------------------------------------------------------------------------------//
// These functions can be employed to alter the algorithm parameters N and D during lifetime of a 
// network object. This allows for iteratively increasing of D. They completely take care of all required
// updated within the network, but they should be needed only a few times.
//---------------------------------------------------------------------------------------------------//

void network::setParameterNSweeps(int Nnew){
  nSweeps=Nnew;
}

//---------------------------------------------------------------------------------------------------//

int network::setParameterD(int Dnew){
  //Decreasing D is neither required nor reasonable in any context
  if(Dnew<D){
    return -1;
  }
  networkState.setParameterD(Dnew);
  //Adapt D
  D=Dnew;
  pars.D=Dnew;
  //Adapt L and R expression
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
    networkState.rightNormalizeState(i);
  }
  networkState.normalizeFinal(1);
  //Allocate a 4D array of dimension LxDxDwxD which is contigous in the last index for the partial contractions Lctr and Rctr
  Lctr.global_access(0,0,0,0)=1;
  //In preparation of the first sweep, generate full contraction to the right (first sweeps starts at site 0)
  calcCtrFull(1);
  for(int iSweep=0;iSweep<nSweeps;iSweep++){
    cout<<"Starting rightsweep\n";
    for(int i=0;i<(L-1);i++){
      //Step of leftsweep
      cout<<"Optimizing site matrix\n";
      curtime=clock();
      errRet=optimize(i,lambda); 
      curtime=clock()-curtime;
      cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      networkState.leftNormalizeState(i);
      calcCtrIterLeft(i+1);
    }
    networkState.normalizeFinal(0);
    cout<<"Starting leftsweep\n";
    for(int i=L-1;i>0;i--){
      //Step of rightsweep
      cout<<"Optimizing site matrix\n";      
      curtime=clock();
      errRet=optimize(i,lambda); 
      curtime=clock()-curtime;
      cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      networkState.rightNormalizeState(i);
      calcCtrIterRight(i-1);
    }
    networkState.normalizeFinal(1);
  }
  cout<<"Determined ground state energy of: "<<lambda<<" using parameters: D="<<D<<" nSweeps="<<nSweeps<<endl;
  return lambda;
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int const i, double &iolambda){
  //Invokes ARPACK++ to solve the eigenvalue problem
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  arcomplex<double> *currentM;
  arcomplex<double> *RTerm, *LTerm, *HTerm;
  //Get the current partial contractions and site matrix of the Hamiltonian
  Lctr.subContractionStart(LTerm,i);
  Rctr.subContractionStart(RTerm,i);
  networkH.subMatrixStart(HTerm,i);
  //Generate matrix which is to be passed to ARPACK++
  optHMatrix HMat(RTerm,LTerm,HTerm,pars,i);
  plambda=&lambda;
  //Using the current site matrix as a starting point allows for much faster convergence as it has already been optimized in previous sweeps (except for the first sweep, this is where a good starting point has to be guessed
  networkState.subMatrixStart(currentM,i);
  ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,&optHMatrix::MultMv,"SR",0,1e-4,400,currentM);
  //One should avoid to hit the maximum number of iterations since this can lead into a suboptimal site matrix, increasing the current energy (although usually not by a lot)
  //So far it seems that the eigensolver either converges quite fast or not at all (i.e. very slow, such that the maximum number of iterations is hit) depending strongly on the tolerance
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
  int DR, DL, DwR, DwL, ld;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr, *targetPctr;
  //container arrays to significantly reduce computational effort by storing intermediate results
  Lctr.subContractionStart(sourcePctr,i-1);
  Lctr.subContractionStart(targetPctr,i);
  if(i==1){
    // b_-1  can only take one value
    DwL=1;
  }
  DR=networkState.locDimR(i-1);
  DL=networkState.locDimL(i-1);
  ld=locd(i-1);
  tmpContainer<lapack_complex_double> innercontainer(ld,DwL,DL,DR);
  tmpContainer<lapack_complex_double> outercontainer(ld,DwR,DR,DL);
  //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
  for(int sip=0;sip<ld;sip++){
    for(int bim=0;bim<DwL;bim++){
      for(int aim=0;aim<DL;aim++){
	for(int aip=0;aip<DR;aip++){
	  simpleContainer=0;
	  for(int aimp=0;aimp<DL;aimp++){
	    simpleContainer+=sourcePctr[pctrIndex(aim,bim,aimp)]*networkState.global_access(i-1,sip,aip,aimp);
	  }
	  innercontainer.global_access(sip,bim,aim,aip)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<ld;si++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	for(int aim=0;aim<DL;aim++){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;sip++){
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
	for(int si=0;si<ld;si++){
	  for(int aim=0;aim<DL;aim++){
	    simpleContainer+=conj(networkState.global_access(i-1,si,ai,aim))*outercontainer.global_access(si,bi,aip,aim);
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
  int DR, DL, DwR, DwL, ld;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double simpleContainer;
  lapack_complex_double *sourcePctr, *targetPctr;
  Rctr.subContractionStart(sourcePctr,i+1);
  Rctr.subContractionStart(targetPctr,i);
  if(i==L-1){
    // b_L-1  can only take one value
    DwR=1;
  }
  DR=networkState.locDimR(i+1);
  DL=networkState.locDimL(i+1);
  ld=locd(i+1);
  tmpContainer<lapack_complex_double> innercontainer(ld,DwR,DR,DL);
  tmpContainer<lapack_complex_double> outercontainer(ld,DwL,DL,DR);
  for(int sip=0;sip<ld;sip++){                                                       
    for(int bi=0;bi<DwR;bi++){
      for(int ai=0;ai<DR;ai++){
	for(int aimp=0;aimp<DL;aimp++){
	  simpleContainer=0;
	  for(int aip=0;aip<DR;aip++){
	    simpleContainer+=sourcePctr[pctrIndex(ai,bi,aip)]*networkState.global_access(i+1,sip,aip,aimp); 
	  }
	  innercontainer.global_access(sip,bi,ai,aimp)=simpleContainer;
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<ld;si++){
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
  for(int aim=0;aim<DL;aim++){
    for(int bim=0;bim<DwL;bim++){
      for(int aimp=0;aimp<DL;aimp++){
	simpleContainer=0;
	for(int si=0;si<ld;si++){
	  for(int ai=0;ai<DR;ai++){
	    simpleContainer+=conj(networkState.global_access(i+1,si,ai,aim))*outercontainer.global_access(si,bim,aimp,ai);
	  }
	}
	targetPctr[pctrIndex(aim,bim,aimp)]=simpleContainer;
      }
    }
  }
  cout<<"Completed calculation of partial contraction"<<endl;
}

//---------------------------------------------------------------------------------------------------//
  
int network::calcCtrFull(int const direction){
  //Full calculation of the contraction is only required once: before the first sweep
  //This is just some ordinary iterative computation of the partial contraction Pctr (P=R,L)
  lapack_complex_double *initialCtr;
  if(direction==1){
    Rctr.global_access(L-1,0,0,0)=lapack_make_complex_double(1.0,0.0);
    for(int i=L-2;i>=0;i--){
      calcCtrIterRight(i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Lctr.global_access(0,0,0,0)=lapack_make_complex_double(1.0,0.0);
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

void network::leftEnrichment(int const i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Bnew;
  int lDR,lDL, lDRR, ld, ldp, lDwR, lDwL;
  lDR=networkState.locDimR(i);
  lDL=networkState.locDimL(i);
  lDRR=networkState.locDimR(i+1);
  ld=locd(i);
  ldp=locd(i+1);
  lDwR=networkH.locDimR(i);
  lDwL=networkH.locDimL(i);
  int MNumCols=lDR*(1+lDwR);
  int MNumRows=ld*lDL;
  Mnew=new lapack_complex_double[MNumRows*MNumCols];
  Bnew=new lapack_complex_double[ldp*lDRR*MNumCols];
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+si*lDL+ai*MNumRows]=networkState.global_access(i,si,ai,aim);
      }
    }      
    for(int air=0;air<lDRR;++air){
      for(int ai=0;ai<lDR;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=networkState.global_access(i+1,si,air,ai);
      }
    }
  }
  //Add zeros and P-Expression to Mnew and Bnew
  //Singular Value Decomposition of Mnew=U*S*V
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  char uplo=(MNumRows>=MNumCols)?'U':'L';
  lapack_int info;
  double *diags=new double[containerDim];
  double *offdiags=new double[containerDim-1];
  lapack_complex_double *QContainer=new lapack_complex_double[lDL*lDL*ld*ld];
  lapack_complex_double *PContainer=new lapack_complex_double[MNumCols*MNumCols];
  lapack_complex_double *Mnewcpy=new lapack_complex_double[MNumCols*MNumCols];
  info=LAPACKE_zgebrd(LAPACK_COL_MAJOR,MNumRows,MNumCols,Mnew,MNumRows,diags,offdiags,QContainer,PContainer);
  arraycpy(MNumCols,MNumCols,Mnew,Mnewcpy);
  info=LAPACKE_zungbr(LAPACK_COL_MAJOR,'Q',MNumRows,MNumRows,MNumCols,Mnew,MNumRows,QContainer);
  info=LAPACKE_zungbr(LAPACK_COL_MAJOR,'P',MNumCols,MNumCols,MNumRows,Mnewcpy,MNumRows,PContainer);
  info=LAPACKE_zbdsqr(LAPACK_COL_MAJOR,uplo,containerDim,MNumCols,MNumRows,0,diags,offdiags,Mnewcpy,MNumCols,Mnew,MNumRows,0,1);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  for(int mi=0;mi<MNumCols;++mi){
    for(int ai=0;ai<lDR;++ai){
      Mnew[ai+lDR*mi]*=diags[ai];
    }
  }
  //From here, Mnew is to be treated as a lDR x MNumCols matrix
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL): the latter is a bit tricky because it is only necessary for i>icrit, might be achieved within lapack
  lapack_complex_double *AStart;
  lapack_complex_double *BStart;
  networkState.subMatrixStart(AStart,i);
  networkState.subMatrixStart(BStart,i+1);
  //Postprocessing: A=U, B=S*V*Bnew (use the networkState subMatrixStart for getting B)
  for(si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=Mnewcpy[aim+si*lDL+ai*MNumRows];
      }
    }
  }
  delete[] diags;
  delete[] Mnewcpy;
  delete[] QContainer;
  delete[] PContainer;
  delete[] Mnew;
  delete[] Bnew;
}


//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int network::locd(int const i){
  return d;
}

int network::locDMax(int const i){
  if(i<=L/2){
    return pow(d,i+1);
  }
  return pow(d,L-i);
}

//---------------------------------------------------------------------------------------------------//
// These functions compute the normalization of the network to the left of some site i. This is for
// testing purpose only since the correct algorithm ensures the results to be trivial.
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
  for(int aim=0;aim<networkState.locDimL(i);aim++){
    for(int aimp=0;aimp<networkState.locDimL(i);aimp++){
      psiContainer=0;
      for(int si=0;si<d;si++){
	for(int aimm=0;aimm<networkState.locDimL(i-1);aimm++){
	  for(int aimmp=0;aimmp<networkState.locDimL(i-1);aimmp++){
	    psiContainer+=networkState.global_access(i-1,si,aim,aimm)*conj(networkState.global_access(i-1,si,aimp,aimm))*psi[(i-1)*D*D+aimm*D+aimmp];
	  }
	}
      }
      psi[i*D*D+aim*D+aimp]=psiContainer;
    }
  }	
}
