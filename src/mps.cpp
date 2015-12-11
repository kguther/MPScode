#include "mps.h"
#include "arrayprocessing.h"
#include "arraycreation.h"
#include <iostream>

mps::mps():stateArray(),
	   aiBlockIndices(0),
	   siaimBlockIndices(0)
{}

//---------------------------------------------------------------------------------------------------//

mps::mps(int const din, int const Din, int const Lin, std::vector<quantumNumber> *conservedQNsin){
  stateArray::initialize(din,Din,Lin);
  conservedQNs=conservedQNsin;
  createInitialState();
}

mps::~mps(){
  delete[] aiBlockIndices;
  delete[] siaimBlockIndices;
}

//---------------------------------------------------------------------------------------------------//

void mps::generate(int const din, int const Din, int const Lin, std::vector<quantumNumber> *conservedQNsin){
  stateArray::generate(din,Din,Lin);
  conservedQNs=conservedQNsin;
  nQNs=(*conservedQNs).size();
  createInitialState();
  basisQNOrderMatrix indexCalc(dimInfo,conservedQNs);
  aiBlockIndices=new std::vector<std::vector<int> >[L];
  siaimBlockIndices=new std::vector<std::vector<multInt> >[L];
  for(int i=0;i<L;++i){
    indexCalc.blockStructure(i,0,aiBlockIndices[i],siaimBlockIndices[i]);
  }
}

//---------------------------------------------------------------------------------------------------//

void mps::createInitialState(){
  int lDL, lDR, ld;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    lDL=(lDR<lDL)?lDR:lDL;
    for(int si=0;si<d;++si){
      for(int aim=0;aim<lDL;++aim){
	if(si==0 && aim==0)
	state_array_access_structure[i][si][aim][aim]=1;
      }
    }
  }
  if(nQNs){
    int qnCriteriumCheck;
    for(int i=0;i<L;++i){
      ld=locd(i);
      lDR=locDimR(i);
      lDL=locDimL(i);
      for(int si=0;si<ld;++si){
	for(int ai=0;ai<lDR;++ai){
	  for(int aim=0;aim<lDL;++aim){
	    qnCriteriumCheck=0;
	    for(int iQN=0;iQN<nQNs;++iQN){
	      qnCriteriumCheck+=(*conservedQNs)[iQN].qnCriterium(i,si,ai,aim);
	    }
	    if(qnCriteriumCheck){
	      global_access(i,si,ai,aim)=0;
	    }
	    else{
	      global_access(i,si,ai,aim)=1;
	    }
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// The following functions are left/right normalizing the matrices of a site after 
// optimization and multiplying the remainder to the matrices of the site to the left/right.
//---------------------------------------------------------------------------------------------------//

int mps::leftNormalizeState(int const i){
  if(nQNs){
    return leftNormalizeStateBlockwise(i);
  }
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
  if(nQNs){
    //return rightNormalizeStateBlockwise(i);
  }
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


//---------------------------------------------------------------------------------------------------//

int mps::leftNormalizeStateBlockwise(int const i){
  int ld, lDR, lDL, lDRR;
  int lBlockSize,rBlockSize, minBlockSize;
  int aiCurrent, siCurrent, aimCurrent, aipCurrent;
  lapack_complex_double *M, *R;
  lapack_complex_double *Rcontainer, *Qcontainer;
  lapack_complex_double *inputA;
  lapack_complex_double zone=1.0, zzero=0.0;
  lapack_int info;
  lDL=locDimL(i);
  lDR=locDimR(i);
  lDRR=locDimR(i+1);
  ld=locd(i);
  R=new lapack_complex_double[lDR*lDR];
  for(int iBlock=0;iBlock<aiBlockIndices[i].size();++iBlock){
    rBlockSize=aiBlockIndices[i][iBlock].size();
    lBlockSize=siaimBlockIndices[i][iBlock].size();
    std::cout<<lBlockSize<<"\t"<<rBlockSize<<std::endl;
    //rBlockSize is required to be smaller than or equal to lBlockSize
    minBlockSize=(lBlockSize<rBlockSize)?lBlockSize:rBlockSize;
    if(rBlockSize!=0 && lBlockSize!=0){
      M=new lapack_complex_double[lBlockSize*rBlockSize];
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  aiCurrent=aiBlockIndices[i][iBlock][j];
	  siCurrent=siaimBlockIndices[i][iBlock][k].si;
	  aimCurrent=siaimBlockIndices[i][iBlock][k].aim;
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      Rcontainer=new lapack_complex_double[rBlockSize*rBlockSize];
      Qcontainer=new lapack_complex_double[rBlockSize];
      info=LAPACKE_zgeqrf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgeqrf: "<<info<<std::endl;
	//exit(1);
      }
      lowerdiag(rBlockSize,rBlockSize,M,Rcontainer,rBlockSize);
      info=LAPACKE_zungqr(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungqr: "<<info<<std::endl;
	//exit(1);
      }
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=aiBlockIndices[i][iBlock][j];
	for(int k=0;k<lBlockSize;++k){
	  siCurrent=siaimBlockIndices[i][iBlock][k].si;
	  aimCurrent=siaimBlockIndices[i][iBlock][k].aim;
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[k+j*lBlockSize];
	}
	for(int l=0;l<rBlockSize;++l){
	  aipCurrent=aiBlockIndices[i][iBlock][l];
	  R[aipCurrent+aiCurrent*lDR]=Rcontainer[l+j*rBlockSize];
	}
      }
      delete[] Rcontainer;
      delete[] Qcontainer;
      delete[] M;
    }
  }
  inputA=new lapack_complex_double[lDR*lDRR];
  for(int si=0;si<ld;++si){
    arraycpy(lDRR,lDR,state_array_access_structure[i+1][si][0],inputA);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,lDR,&zone,R,lDR,inputA,lDR,&zzero,state_array_access_structure[i+1][si][0],lDR);
  }
  delete[] inputA;
  delete[] R;
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int mps::rightNormalizeStateBlockwise(int const i){
  int ld, lDR, lDL, lDLL;
  int lBlockSize,rBlockSize, minBlockSize;
  int aiCurrent, siCurrent, aimCurrent, aipCurrent;
  lapack_complex_double *M, *R;
  lapack_complex_double *Rcontainer, *Qcontainer;
  lapack_complex_double *inputA;
  lapack_complex_double zone=1.0, zzero=0.0;
  lapack_int info;
  lDL=locDimL(i);
  lDR=locDimR(i);
  lDLL=locDimL(i-1);
  ld=locd(i);
  R=new lapack_complex_double[lDR*lDR];
  for(int iBlock=0;iBlock<aiBlockIndices[i].size();++iBlock){
    rBlockSize=aiBlockIndices[i][iBlock].size();
    lBlockSize=siaimBlockIndices[i][iBlock].size();
    std::cout<<lBlockSize<<"\t"<<rBlockSize<<std::endl;
    //lBlockSize is required to be smaller than or equal to rBlockSize
    // THUS, ONLY lBlockSize=rBlockSize ENABLES BOTH NORMALIZATION FUNCTIONS WITH A SINGLE BLOCKSIZE
    minBlockSize=(lBlockSize<rBlockSize)?lBlockSize:rBlockSize;
    if(rBlockSize!=0 && lBlockSize!=0){
      M=new lapack_complex_double[lBlockSize*rBlockSize];
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  aiCurrent=aiBlockIndices[i][iBlock][j];
	  siCurrent=siaimBlockIndices[i][iBlock][k].si;
	  aimCurrent=siaimBlockIndices[i][iBlock][k].aim;
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      Rcontainer=new lapack_complex_double[rBlockSize*rBlockSize];
      Qcontainer=new lapack_complex_double[rBlockSize];
      info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgeqrf: "<<info<<std::endl;
	exit(1);
      }
      lowerdiag(lBlockSize,lBlockSize,M,Rcontainer);
      info=LAPACKE_zungrq(LAPACK_COL_MAJOR,minBlockSize,rBlockSize,rBlockSize,M,minBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungqr: "<<info<<std::endl;
	exit(1);
      }
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=aiBlockIndices[i][iBlock][j];
	for(int k=0;k<lBlockSize;++k){
	  siCurrent=siaimBlockIndices[i][iBlock][k].si;
	  aimCurrent=siaimBlockIndices[i][iBlock][k].aim;
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[k+j*lBlockSize];
	}
	for(int l=0;l<rBlockSize;++l){
	  aipCurrent=aiBlockIndices[i][iBlock][l];
	  R[aipCurrent+aiCurrent*lDR]=Rcontainer[l+j*rBlockSize];
	}
      }
      delete[] Rcontainer;
      delete[] Qcontainer;
      delete[] M;
    }
  }
  inputA=new lapack_complex_double[lDL*lDLL];
  for(int si=0;si<ld;++si){
    arraycpy(lDLL,lDL,state_array_access_structure[i-1][si][0],inputA);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,lDL,&zone,inputA,lDLL,R,lDL,&zzero,state_array_access_structure[i-1][si][0],lDLL);
  }
  delete[] inputA;
  delete[] R;
  return 0;
}
