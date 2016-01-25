#include "mps.h"
#include "arrayprocessing.h"
#include "arraycreation.h"
#include <iostream>

mps::mps():stateArray()
{}

//---------------------------------------------------------------------------------------------------//

mps::mps(dimensionTable &dimInfoIn, std::vector<quantumNumber> *conservedQNsin){
  stateArray::initialize(dimInfoIn);
  conservedQNs=conservedQNsin;
  if(conservedQNs){
    nQNs=(*conservedQNs).size();
  }
  else{
    nQNs=0;
  }
  if(nQNs){
    indexTable.initialize(dimInfo,conservedQNs);
    indexTable.generateQNIndexTables();
  }
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

mps::~mps(){
}

//---------------------------------------------------------------------------------------------------//

void mps::generate(dimensionTable &dimInfoIn, std::vector<quantumNumber> *conservedQNsin){
  int direction;
  stateArray::generate(dimInfoIn);
  conservedQNs=conservedQNsin;
  if(conservedQNs){
    nQNs=(*conservedQNs).size();
  }
  else{
    nQNs=0;
  }
  if(nQNs){
    indexTable.initialize(dimInfo,conservedQNs);
    indexTable.generateQNIndexTables();
  }
  createInitialState();
}

//---------------------------------------------------------------------------------------------------//

void mps::createInitialState(){
  int lDL, lDR, ld;
  if(!nQNs){
    for(int i=0;i<L;++i){
      lDL=locDimL(i);
      lDR=locDimR(i);
      lDL=(lDR<lDL)?lDR:lDL;
      ld=locd(i);
      for(int si=0;si<ld;++si){
	for(int aim=0;aim<lDL;++aim){
	  state_array_access_structure[i][si][aim][aim]=1;
	}
      }
    }
  }
  else{
    //This is the exact ground state at the critical point for fixed particle number and subchain parity. It turns out that this is a nice guess for the ground state of the perturbed system (for small perturbations).
    int numBlocks, lBlockSize, rBlockSize;
    for(int i=0;i<L;++i){
      numBlocks=indexTable.numBlocksLP(i);
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	rBlockSize=indexTable.rBlockSizeLP(i,iBlock);
	lBlockSize=indexTable.lBlockSizeLP(i,iBlock);
	for(int j=0;j<rBlockSize;++j){
	  for(int k=0;k<lBlockSize;++k){
	    state_array_access_structure[i][indexTable.siBlockIndexLP(i,iBlock,k)][indexTable.aiBlockIndexLP(i,iBlock,j)][indexTable.aimBlockIndexLP(i,iBlock,k)]=exactGroundStateEntry(i,indexTable.siBlockIndexLP(i,iBlock,k),indexTable.aiBlockIndexLP(i,iBlock,j),indexTable.aimBlockIndexLP(i,iBlock,k));
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

lapack_complex_double mps::exactGroundStateEntry(int const i, int const si, int const ai, int const aim){
  if(si==0 || si==2){
    return 1.0;
  }
  if((*conservedQNs)[1].QNLabel(i-1,aim)==1){
    return 1.0;
  }
  return -1.0;
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
  transp(D2,D2,Rcontainer);
  for(int si=0;si<ld;++si){
    transp(D1,D2,state_array_access_structure[i][si][0]);
    cblas_ztrmm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,D2,D3,&zone,Rcontainer,D2,state_array_access_structure[i+1][si][0],D2);
    //here, R is packed into the matrices of the next site
  }                                                //POSSIBLE TESTS: TEST FOR Q*R - DONE: WORKS THE WAY INTENDED
  delete[] Rcontainer;
  delete[] Qcontainer;
  return 0;  //TODO: Add exception throw
}

//---------------------------------------------------------------------------------------------------//

int mps::rightNormalizeState(int const i){
  if(nQNs){
    return rightNormalizeStateBlockwise(i);
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
  if(info){
    std::cout<<"ERROR IN LAPACKE_zungrq:"<<info<<" At site: "<<i<<" With dimensions: "<<D2<<"x"<<D1<<" and local Hilbert space dimension: "<<ld<<std::endl;
    exit(1);
  }
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
// Now, these two functions also left-/rightnormalize a site matrix, but they use the block structure
// of the MPS matrices when using conserved QNs. In particular, they only make sense when using conserved
// QNs. There, the normalization is executed on each block individually, preserving the QN constraint. 
// IMPORTANT: IN GENERAL, THE BLOCKS ARE NOT SQUARE. THE QN LABELING SCHEME MUST BE DESIGNED SUCH
// THAT EACH BLOCK CAN BE BROUGHT INTO CANONICAL FORM. THIS IS A NONTRIVIAL CONSTRAINT.
//---------------------------------------------------------------------------------------------------//

int mps::leftNormalizeStateBlockwise(int const i){
  int ld, lDR, lDL, lDRR;
  int lBlockSize,rBlockSize;
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
  for(int iBlock=0;iBlock<indexTable.numBlocksLP(i);++iBlock){
    rBlockSize=indexTable.rBlockSizeLP(i,iBlock);
    lBlockSize=indexTable.lBlockSizeLP(i,iBlock);
    //rBlockSize is required to be smaller than or equal to lBlockSize
    if(rBlockSize!=0 && lBlockSize!=0){
      //We do not normalize empty blocks
      M=new lapack_complex_double[lBlockSize*rBlockSize];
      //First, copy the content of the block to some dummy array
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesLP(i,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      Rcontainer=new lapack_complex_double[rBlockSize*rBlockSize];
      Qcontainer=new lapack_complex_double[rBlockSize];
      //Make a QR decomposition of that dummy array
      info=LAPACKE_zgeqrf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgeqrf: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
	exit(1);
      }
      lowerdiag(rBlockSize,rBlockSize,M,Rcontainer,lBlockSize);
      info=LAPACKE_zungqr(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungqr: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
	exit(1);
      }
      //Copy the resulting Q-matrix back into the block of the MPS matrix
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesLP(i,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[k+j*lBlockSize];
	}
	//And the resulting R-matrix into the corresponding elements of a uncompressed matrix R
	for(int l=0;l<rBlockSize;++l){
	  aiCurrent=indexTable.aiBlockIndexLP(i,iBlock,j);
	  aipCurrent=indexTable.aiBlockIndexLP(i,iBlock,l);
	  R[aipCurrent+aiCurrent*lDR]=Rcontainer[l+j*rBlockSize];
	}
      }
      delete[] Rcontainer;
      delete[] Qcontainer;
      delete[] M;
    }
  }
  //Finally, multiply the uncompressed matrix R into the next MPS matrix. Block structure does not need to be adressed explicitly, since they both have the required structure.
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
  int lBlockSize,rBlockSize;
  int aiCurrent, siCurrent, aimCurrent, aimpCurrent;
  lapack_complex_double *M, *R;
  lapack_complex_double *Rcontainer, *Qcontainer;
  lapack_complex_double *inputA;
  lapack_complex_double zone=1.0, zzero=0.0;
  lapack_int info;
  lDL=locDimL(i);
  lDR=locDimR(i);
  lDLL=locDimL(i-1);
  ld=locd(i);
  R=new lapack_complex_double[lDL*lDL];
  for(int iBlock=0;iBlock<indexTable.numBlocksRP(i);++iBlock){
    lBlockSize=indexTable.lBlockSizeRP(i,iBlock);
    rBlockSize=indexTable.rBlockSizeRP(i,iBlock);
    //lBlockSize is required to be smaller than or equal to rBlockSize
    if(rBlockSize!=0 && lBlockSize!=0){
      M=new lapack_complex_double[lBlockSize*rBlockSize];
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  convertIndicesRP(i,j,k,iBlock,siCurrent,aiCurrent,aimCurrent);
	  M[k+j*lBlockSize]=global_access(i,siCurrent,aiCurrent,aimCurrent);
	}
      }
      Rcontainer=new lapack_complex_double[lBlockSize*lBlockSize];
      Qcontainer=new lapack_complex_double[rBlockSize];
      info=LAPACKE_zgerqf(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zgerqf: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
	exit(1);
      }
      lowerdiag(lBlockSize,lBlockSize,M+(rBlockSize-lBlockSize)*lBlockSize,Rcontainer);
      info=LAPACKE_zungrq(LAPACK_COL_MAJOR,lBlockSize,rBlockSize,lBlockSize,M,lBlockSize,Qcontainer);
      if(info){
	std::cout<<"Error in LAPACKE_zungrq: "<<info<<" with block sizes: "<<lBlockSize<<"x"<<rBlockSize<<" at site "<<i<<std::endl;
	exit(1);
      }
      for(int j=0;j<lBlockSize;++j){
	for(int k=0;k<rBlockSize;++k){
	  convertIndicesRP(i,k,j,iBlock,siCurrent,aiCurrent,aimCurrent);
	  global_access(i,siCurrent,aiCurrent,aimCurrent)=M[j+k*lBlockSize];
	}
	for(int l=0;l<lBlockSize;++l){
	  aimCurrent=indexTable.aimBlockIndexRP(i,iBlock,j);
	  aimpCurrent=indexTable.aimBlockIndexRP(i,iBlock,l);
	  R[aimpCurrent+aimCurrent*lDL]=Rcontainer[l+j*lBlockSize];
	}
      }
      delete[] Qcontainer;
      delete[] Rcontainer;
      delete[] M;
    }
  }
  inputA=new lapack_complex_double[lDL*lDLL];
  for(int si=0;si<ld;++si){
    arraycpy(lDLL*lDL,state_array_access_structure[i-1][si][0],inputA);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,lDL,&zone,inputA,lDLL,R,lDL,&zzero,state_array_access_structure[i-1][si][0],lDLL);
  }
  delete[] inputA;
  delete[] R;
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// These functions convert block-internal indices to normal MPS bond indices.
//---------------------------------------------------------------------------------------------------//

void mps::convertIndicesLP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim){
  ai=indexTable.aiBlockIndexLP(i,iBlock,j);
  aim=indexTable.aimBlockIndexLP(i,iBlock,k);
  si=indexTable.siBlockIndexLP(i,iBlock,k);
}

void mps::convertIndicesRP(int const i, int const j, int const k, int const iBlock, int &si, int &ai, int &aim){
  aim=indexTable.aimBlockIndexRP(i,iBlock,k);
  ai=indexTable.aiBlockIndexRP(i,iBlock,j);
  si=indexTable.siBlockIndexRP(i,iBlock,j);
}

//---------------------------------------------------------------------------------------------------//
// Outdated function to enforce the QN constraint on some MPS.
//---------------------------------------------------------------------------------------------------//

void mps::restoreQN(int const i){
  int lDL, lDR, ld;
  for(int iQN=0;iQN<nQNs;++iQN){
    ld=locd(i);
    lDR=locDimR(i);
    lDL=locDimL(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  if(((*conservedQNs)[iQN].QNLabel(i,ai)-(*conservedQNs)[iQN].QNLabel(i-1,aim)-(*conservedQNs)[iQN].QNLabel(si)) && abs(global_access(i,si,ai,aim))>0.000001){
	    global_access(i,si,ai,aim)=0;
	  }
	}
      }
    }
  }
}
