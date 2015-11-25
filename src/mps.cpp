#include "mps.h"
#include "arrayprocessing.h"
#include "arraycreation.h"

mps::mps():stateArray()
{
}

//---------------------------------------------------------------------------------------------------//

mps::mps(int const din, int const Din, int const Lin){
  stateArray::initialize(din,Din,Lin);
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

void mps::generate(int const din, int const Din, int const Lin){
  stateArray::generate(din,Din,Lin);
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

void mps::createInitialState(){
    int lDL, lDR;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    for(int si=0;si<d;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  if(ai==aim){
	    state_array_access_structure[i][si][ai][aim]=1;
	  }
	}
      }
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

int mps::leftNormalizeState(int const i){
  lapack_int info;
  int D1, D2, D3, ld;
  //Yep, the local dimensions are just named D1, D2, D3 - the decomposition is of a d*D1xD2 matrix
  D1=locDimL(i);
  D2=locDimR(i);
  D3=locDimR(i+1);
  ld=locd(i);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=1.0;
  Qcontainer=new lapack_complex_double[D2];//Used for storage of lapack-internal matrices
  Rcontainer=new lapack_complex_double[D2*D2];//Used for storage of R from RQ decomposition
  //Enable use of LAPACK_ROW_MAJOR which is necessary here due to the applied storage scheme
  for(int si=0;si<ld;++si){
    transp(D2,D1,state_array_access_structure[i][si][0]);
  }
  //Use thin QR decomposition
  info=LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,ld*D1,D2,state_array_access_structure[i][0][0],D2,Qcontainer);
  upperdiag(D2,D2,state_array_access_structure[i][0][0],Rcontainer);                               
  //Only first D2 columns are used -> thin QR (below icrit, this is equivalent to a full QR)
  info=LAPACKE_zungqr(LAPACK_ROW_MAJOR,ld*D1,D2,D2,state_array_access_structure[i][0][0],D2,Qcontainer);
  for(int si=0;si<ld;++si){
    transp(D1,D2,state_array_access_structure[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,state_array_access_structure[i+1][si][0],D2); 
    //REMARK: Use CblasTrans because Rcontainer is in row_major while state_array_access_structure[i+1][si][0] is in column_major order - this is a normal matrix multiplication - here, R is packed into the matrices of the next site
  }                                                //POSSIBLE TESTS: TEST FOR Q*R - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

int mps::rightNormalizeState(int const i){
  lapack_int info;
  //This time, we decompose a D2xd*D1 matrix and multiply the D2xD2 matrix R with a D3xD2 matrix
  int D1,D2,D3,ld;
  D1=locDimR(i);
  D2=locDimL(i);
  D3=locDimL(i-1);
  ld=locd(i);
  lapack_complex_double *Rcontainer, *Qcontainer;
  const lapack_complex_double zone=1.0;
  Qcontainer=new lapack_complex_double[ld*D1];
  Rcontainer=new lapack_complex_double[D2*D2];
  //Thats how zgerqf works: the last D2 columns contain the upper trigonal matrix R, to adress them, move D2 from the end
  info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,D2,ld*D1,state_array_access_structure[i][0][0],D2,Qcontainer);
  //lowerdiag does get an upper trigonal matrix in column major ordering, dont get confused
  lowerdiag(D2,D2,state_array_access_structure[i][0][0]+D2*(ld*D1-D2),Rcontainer);
  info=LAPACKE_zungrq(LAPACK_COL_MAJOR,D2,ld*D1,D2,state_array_access_structure[i][0][0],D2,Qcontainer);
  for(int si=0;si<ld;++si){
    cblas_ztrmm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,D3,D2,&zone,Rcontainer,D2,state_array_access_structure[i-1][si][0],D3);
  }                                                //POSSIBLE TESTS: TEST FOR R*Q - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

void mps::normalizeFinal(int const i){
  lapack_complex_double normalization;
  int site, lcD,ld;
  if(i){
    site=0;
    lcD=locDimR(0);
  }
  else{
    site=L-1;
    lcD=locDimL(L-1);
  }
  ld=locd(site);
  //Normalize last matrices to maintain normalization of state
  normalization=cblas_dznrm2(ld*lcD,state_array_access_structure[site][0][0],1);
  normalization=1.0/normalization;
  cblas_zscal(ld*lcD,&normalization,state_array_access_structure[site][0][0],1);
}
