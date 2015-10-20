#include <cstdlib>
#include <complex>
#include <iostream>
#include <cblas.h>
#include <arcomp.h>
#include <arlnsmat.h>
#include <arlscomp.h>
namespace LPK{          //Horrible hack to prevent conflict with superLU (why does this even work?)
#include <lapacke.h>
#include <lapacke_utils.h>
}
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"

using namespace std;
using namespace LPK;

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//---------------------------------------------------------------------------------------------------//
// Main file containing the MPS ansatz solver itself
//---------------------------------------------------------------------------------------------------//

network::network(parameters inputpars){
  pars=inputpars;
  d=inputpars.d;
  D=inputpars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  N=inputpars.N;
  icrit=L/2;                                                                      //In case chain is too short, this is the correct value (trust me)
  for(int i=0;i<L;i++){
    if(pow(d,i+1)>D){
      icrit=i;
      break;
    }
  }
  create5D(L,d,d,Dw,Dw,&networkH);
  createStateArray(d,D,L,&networkState);
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
  arcomplex<double> lambda;
  for(int i=L-1;i>0;i--){
    rightNormalizeState(i);
  }
  normalizeFinal(1);
  create4D(L,D,Dw,D,&Lctr);                                                        //Allocates a 4D array of dimension LxDxDwxD which is contigous in the last index
  create4D(L,D,Dw,D,&Rctr);                                                        //for the partial contractions Lctr and Rctr
  Lctr[0][0][0][0]=1;                                                              //initialization neccessary for iterative building of Lctr
  calcCtrFull(Rctr,1);                                                             //Preparing for the first sweep 
  for(int nsweep=0;nsweep<N;nsweep++){
    for(int i=0;i<(L-1);i++){
      lambda=optimize(i);                                                          //Eigenvalue solver to actualize cstate
      leftNormalizeState(i);
      calcCtrIter(Lctr,-1,i+1);
    }
    normalizeFinal(0);
    for(int i=L-1;i>0;i--){
      lambda=optimize(i);
      rightNormalizeState(i);
      calcCtrIter(Rctr,1,i-1);
    }
    normalizeFinal(1);
  }
  cout<<Lctr[1][0][0][0]<<endl;
  delete4D(&Lctr);                                                                 //Frees the memory previously allocated for partial contractions
  delete4D(&Rctr);
  return real(lambda);
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

arcomplex<double> network::optimize(int i){
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  plambda=&lambda;
  vector<int> irow;
  vector<int> pcol;
  int nzel, nconv;
  int dimension=d*locDimR(i)*locDimL(i);
  vector<arcomplex<double> > H;
  buildMatrix(i,&irow, &pcol, &nzel, &H);
  ARluNonSymMatrix<arcomplex<double> ,double> Hopt(dimension,nzel,&H[0],&irow[0],&pcol[0]);
  ARluCompStdEig<double> OptProblem(1,Hopt,"SM",3,1e-8);
  nconv=OptProblem.EigenValVectors(networkState[i][0][0],plambda);
  if(nconv==0){
    cout<<"Failed to converge in iterative eigensolver";
    exit(-1);
  }
  cout<<"Current energy: "<<real(lambda)<<endl;
  return lambda;
}

void network::buildMatrix(int i, vector<int> *irow, vector<int> *pcol, int *nzel, vector<arcomplex<double> > *H){
  //This function generates the sparsed matrix which represents the eigenvalue problem in CSC format. H contains the nonzero entries, irow the row indices and pcol the indices of column breaks in H
  //std::vector is used because the number of nonzero elements is unknown (and can be quite large, to guess a maximum would require for allocation of memory of size O(d^2D^4) which quickly exceeds RAM)
  //Arguments are of type vector because internal use of vectors would free the allocated memory upon end of function - they can still be handed over to ARPACK++ classes via &irow[0] etc (yes, this works)
  int lDR,lDL,DwL,DwR,colIndexCounter;
  double threshold=1e-50;  //Entries with absolute value below this are omitted, might be changed to a parameter later on, no neccessity now
  arcomplex<double> container;
  DwL=Dw;
  DwR=Dw;
  if(i==0){                                          
    DwL=1;
  }
  if(i==(L-1)){                                       
    DwR=1;
  }
  cout<<"Calculating sparse matrix H"<<endl;
  lDR=locDimR(i);
  lDL=locDimL(i);
  (*pcol).push_back(0);
  colIndexCounter=0;
  for(int si=0;si<d;si++){
    for(int ai=0;ai<lDR;ai++){
      for(int aim=0;aim<lDL;aim++){
	for(int sip=0;sip<d;sip++){
	  for(int aip=0;aip<lDR;aip++){
	    for(int aimp=0;aimp<lDL;aimp++){
	      container=0;
	      for(int bi=0;bi<DwR;bi++){
		for(int bim=0;bim<DwL;bim++){
		  container+=Lctr[i][aim][bim][aimp]*networkH[i][si][sip][bi][bim]*Rctr[i][ai][bi][aip];
		}
	      }
	      if(abs(container)>threshold){
		colIndexCounter++;                            //Counts to the next column break
		(*H).push_back(container);                    //Adds the next nonzero element
		(*irow).push_back(aimp+lDL*aip+lDL*lDR*sip);  //Stores the row index of the next nonzero element (current aimp,aip,sip)
	      }
	    }
	  }
	}
	(*pcol).push_back(colIndexCounter);
      }
    }
  }
  (*nzel)=(*H).size();                                        //Number of nonzero elements is required for ARPACK++ classes
  cout<<"Finished calculation of H"<<endl;
}

//---------------------------------------------------------------------------------------------------//
// Following functions are for calculation of partial contractions of the state for calculating the 
// action of the hamiltonian on the state. This is lots of numerical effort and probably the best 
// candidate for optimization
//---------------------------------------------------------------------------------------------------//

void network::calcCtrIter(lapack_complex_double ****Pctr, const int direction, const int i){ //direction==1 builds R expression, direction==-1 builds L expression
  int DR, DL, DwR, DwL;
  int threshold=1e-50;
  DwL=Dw;
  DwR=Dw;
  lapack_complex_double ****innercontainer, ****outercontainer;                     //container arrays to significantly reduce computational effort by storing intermediate results
  if((i==1) && (direction==-1)){                                                    // b_-1  can only take one value
    DwL=1;
  }
  if((i==(L-2)) && (direction==1)){                                                  //as does b_L-1
    DwR=1;
  }
  icrit=L/2;
  for(int iloop=0;iloop<=pars.L;iloop++){
    if(D<pow(d,iloop+1)){
      icrit=iloop;
      break;
    }
    //TODO: Add exception throw
  }
  DR=locDimR(i+direction);
  DL=locDimL(i+direction);
  create4D(d,DwL,DL,DR,&innercontainer);
  create4D(d,DwR,DR,DL,&outercontainer);
  for(int sip=0;sip<d;sip++){                                                        //horrible construct to efficiently compute the partial contraction, is parallelizable, needs to be parallelized (still huge computational effort) <-- potential for optimization
    for(int bim=0;bim<DwL;bim++){
      for(int aim=0;aim<DL;aim++){
	for(int aip=0;aip<DR;aip++){
	  innercontainer[sip][bim][aim][aip]=0;
	  for(int aimp=0;aimp<DL;aimp++){
	    innercontainer[sip][bim][aim][aip]+=Pctr[i+direction][aim][bim][aimp]*networkState[i+direction][sip][aip][aimp]; 
	  }
	}
      }
    }
  }
  cout<<"Completed calculation of inner container"<<endl;
  for(int si=0;si<d;si++){ //calculation of outer container seems to be the computationally most effortive task (fortunately, it is easily optimized)
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	for(int aim=0;aim<DL;aim++){
	  outercontainer[si][bi][aip][aim]=0;
	  for(int sip=0;sip<d;sip++){
	    for(int bim=0;bim<DwL;bim++){
	      if(abs(networkH[i+direction][si][sip][bi][bim])>threshold){               //<- maybe use more sophisticated sparse format
		outercontainer[si][bi][aip][aim]+=networkH[i+direction][si][sip][bi][bim]*innercontainer[sip][bim][aim][aip];
	      }
	    }
	  }
	}
      }
    }
  }
  cout<<"Completed calculation of outer container"<<endl;
  delete4D(&innercontainer);
  for(int ai=0;ai<DR;ai++){
    for(int bi=0;bi<DwR;bi++){
      for(int aip=0;aip<DR;aip++){
	Pctr[i][ai][bi][aip]=0;
	for(int si=0;si<d;si++){
	  for(int aim=0;aim<DL;aim++){
	    Pctr[i][ai][bi][aip]+=conj(networkState[i+direction][si][ai][aim])*outercontainer[si][bi][aip][aim];
	  }
	}
      }
    }
  }
  delete4D(&outercontainer);
}
  
int network::calcCtrFull(lapack_complex_double ****Pctr, const int direction){
  int L=pars.L;
  if(direction==1){
    Pctr[L-1][0][0][0]=lapack_make_complex_double(1,0);
    for(int i=L-2;i>=0;i--){
      calcCtrIter(Pctr,direction,i);
	}
    return 0;
  }
  else{
    if(direction==-1){
      Pctr[0][0][0][0]=lapack_make_complex_double(1,0);
      for(int i=1;i<L;i++){
	calcCtrIter(Pctr,direction,i);
      }
      return 0;
    }
    else{
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

int network::leftNormalizeState(int i){
  lapack_int info;
  int D1, D2, D3;
  D1=locDimL(i);
  D2=locDimR(i);
  D3=locDimR(i+1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0.0);
  Qcontainer=new lapack_complex_double[d*D1];
  Rcontainer=new lapack_complex_double[D2*D2];
  for(int si=0;si<d;si++){
    transp(D1,D2,networkState[i][si][0]);                                              //Enables use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  }
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,d*D1,D2,networkState[i][0][0],D2,Qcontainer);   //Use thin QR decomposition
  upperdiag(D2,D2,networkState[i][0][0],Rcontainer);                               
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,d*D1,D2,D2,networkState[i][0][0],D2,Qcontainer);//Only first D columns are used -> thin QR (below icrit, this is equivalent to a full QR)
  for(int si=0;si<d;si++){
    transp(D2,D1,networkState[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,networkState[i+1][si][0],D2); //REMARK: Use CblasTrans because Rcontainer is in row_major while networkState[i+1][si][0] is in column_major order - this is a normal matrix multiplication
  }                                                //POSSIBLE TESTS: TEST FOR Q*R
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

int network::rightNormalizeState(int i){
  lapack_int info;
  int D1,D2,D3;
  D1=locDimR(i);
  D2=locDimL(i);
  D3=locDimL(i-1);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=lapack_make_complex_double(1.0,0);
  Qcontainer=new lapack_complex_double[d*D1];
  Rcontainer=new lapack_complex_double[D2*D2];
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,d*D1,networkState[i][0][0],D2,Qcontainer);     //Thats how zgerqf works: the last D2 columns contain the upper trigonal matrix R, to adress them, move D2 from the end
  lowerdiag(D2,D2,networkState[i][0][0]+D2*(d*D1-D2),Rcontainer);                        //lowerdiag does get an upper trigonal matrix in column major ordering, dont get confused
  info=LAPACKE_zungrq(LAPACK_COL_MAJOR,D2,d*D1,D2,networkState[i][0][0],D2,Qcontainer);
  for(int si=0;si<d;si++){
    cblas_ztrmm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,D3,D2,&zone,Rcontainer,D2,networkState[i-1][si][0],D3);
  }                                                //POSSIBLE TESTS: TEST FOR R*Q
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

void network::normalizeFinal(int i){
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
// Simple, but essential
//---------------------------------------------------------------------------------------------------//

int network::locDimL(int i){
  if(i<=icrit){
    return pow(d,i);
  }
  if(i<=L-icrit-1){
    return D;
  }
  return pow(d,L-i);
}

int network::locDimR(int i){
  if(i<icrit){
    return pow(d,i+1);
  }
  if(i<=L-icrit-2){
    return D;
  }
  return pow(d,L-i-1);
}
