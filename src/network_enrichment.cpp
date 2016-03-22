#include "mkl_complex_defined.h"
#include "arrayprocessing.h"
#include "tmpContainer.h"
#include "network.h"
#include "optHMatrix.h"
#include "blockHMatrix.h"
#include <cmath>
#include <iostream>
#include <memory>

//---------------------------------------------------------------------------------------------------//
// The implementation of the enrichment is very ugly and might still contain severe bugs which
// are irrelevant to the testing system (Heisenberg chain). It does however, improve convergence
// speed and avoids local minima (at least thats what I expect, there are no problems with local 
// minima till now).
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichment(int i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Bnew;
  lapack_complex_double *pExpression, *unity;
  int lDRR, ldp;
  getLocalDimensions(i);
  lDRR=networkState.locDimR(i+1);
  ldp=locd(i+1);
  int MNumCols=lDR*(1+lDwR);
  int MNumRows=ld*lDL;
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  //Allocate the memory needed for the output of ZUNGBR which is more than the original matrix - also prevents possible segfault when copying (size of target is given for the same reason)
  Mnew=new lapack_complex_double[maxDim*maxDim];
  Bnew=new lapack_complex_double[ldp*lDRR*MNumCols];
  pExpression=new lapack_complex_double[MNumRows*lDR*lDwR];
  getPExpressionLeft(i,pExpression);
  for(int si=0;si<ld;++si){
    //Copy A and B while rearranging to allow for expansion
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
    //Add zeros and P-Expression to Mnew and Bnew
    for(int ai=lDR;ai<MNumCols;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+si*lDL+ai*MNumRows]=alpha*pExpression[aim+si*lDL+(ai-lDR)*MNumRows];
      }
    }
    for(int air=0;air<lDRR;++air){
      for(int ai=lDR;ai<MNumCols;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=0;
      }
    }
  }
  delete[] pExpression;
  //Singular Value Decomposition of Mnew=U*S*VT
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  double *diags=new double[containerDim];
  lapack_complex_double *VT=new lapack_complex_double[maxDim*maxDim];
  lapack_complex_double *U=new lapack_complex_double[maxDim*maxDim];
  //POSSIBLE IMPROVEMENT: USE LAPACK DRIVER ROUTINE INSTEAD OF COMPUTATIONAL ROUTINES
  //lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',MNumRows,MNumCols,Mnew,MNumRows,diags,U,MNumRows,VT,MNumCols);
  //U -> A, S*VT->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int mi=0;mi<MNumCols;++mi){
    for(int ai=0;ai<lDR;++ai){
      //I think this is correct, not the original lDR as leading dimension
      VT[ai+MNumCols*mi]*=diags[ai];
    }
  }
  delete[] diags;
  //From here, Mnewcpy is to be treated as a lDR x MNumCols matrix
  //Postprocessing: A=U, B=S*V*Bnew
  lapack_complex_double *BStart;
  lapack_complex_double *networkB;
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  //Multiply S*V into the expanded B (Bnew) to create the normal-sized B (Bstart, direct access to the networkState mps)
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(networkB,i+1,si);
    BStart=Bnew+si*lDRR*MNumCols;
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,MNumCols,&zone,VT,lDR,BStart,MNumCols,&zzero,networkB,lDR);
  }
  delete[] VT;
  delete[] Bnew;
  networkState.subMatrixStart(networkB,i);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
        networkState.global_access(i,si,ai,aim)=U[aim+si*lDL+ai*lDL*ld];
      }
    }
  }
  delete[] unity;
  delete[] Mnew;
  delete[] U;
}

//---------------------------------------------------------------------------------------------------//

void network::rightEnrichment(int i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Anew;
  lapack_complex_double *pExpression;
  int lDLL, ldm;
  getLocalDimensions(i);
  lDLL=networkState.locDimL(i-1);
  ldm=locd(i-1);
  int MNumCols=ld*lDR;
  int MNumRows=lDL*(1+lDwL);
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  Mnew=new lapack_complex_double[maxDim*maxDim];
  Anew=new lapack_complex_double[ldm*lDLL*MNumRows];
  pExpression=new lapack_complex_double[ld*lDR*lDL*lDwL];
  getPExpressionRight(i,pExpression);
  //Copy A and B while rearranging to allow for expansion
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+ai*MNumRows+si*lDR*MNumRows]=networkState.global_access(i,si,ai,aim);
      }
    }
    for(int aim=0;aim<lDL;++aim){
      for(int aimm=0;aimm<lDLL;++aimm){
	Anew[aimm+aim*ld*lDLL+si*lDLL]=networkState.global_access(i-1,si,aim,aimm);
      }
    }
    //Add zeros and P-Expression to Mnew and Bnew
    for(int ai=0;ai<lDR;++ai){
      for(int aim=lDL;aim<MNumRows;++aim){
	Mnew[aim+ai*MNumRows+si*lDR*MNumRows]=alpha*pExpression[aim-lDL+si*lDR*lDwL*lDL+ai*lDwL*lDL];
      }
    }
    for(int aim=lDL;aim<MNumRows;++aim){
      for(int aimm=0;aimm<lDLL;++aimm){
        Anew[aimm+aim*ld*lDLL+si*lDLL]=0;
      }
    }
  }
  delete[] pExpression;
  //Singular Value Decomposition of Mnew=U*S*V
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  double *diags=new double[containerDim];
  lapack_complex_double *VT=new lapack_complex_double[maxDim*maxDim];
  lapack_complex_double *U=new lapack_complex_double[maxDim*maxDim];
  //lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',MNumRows,MNumCols,Mnew,MNumRows,diags,U,MNumRows,VT,MNumCols);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int aim=0;aim<lDL;++aim){
    for(int mi=0;mi<MNumRows;++mi){
      U[mi+aim*MNumRows]*=diags[aim];
    }
  }
  delete[] diags;
  //From here, Mnew is to be treated as a MNumRows x lDL matrix
  //Postprocessing: A=Anew*U*S, B=VT
  lapack_complex_double *AStart=new lapack_complex_double[lDLL*MNumRows];
  lapack_complex_double *networkA;
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  //Multiply U*S into the expanded A (AStart) to create the normal-sized A (networkA, direct access to the networkState mps)
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(networkA,i-1,si);
    for(int mi=0;mi<MNumRows;++mi){
      for(int aimm=0;aimm<lDLL;++aimm){
	AStart[aimm+mi*lDLL]=Anew[aimm+si*lDLL+mi*ld*lDLL];
      }
    }
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,MNumRows,&zone,AStart,lDLL,U,MNumRows,&zzero,networkA,lDLL);
  }
  delete[] AStart;
  delete[] U;
  delete[] Anew;
  //Store the remaining VT into B
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=VT[aim+si*lDR*MNumCols+ai*MNumCols];
      }
    }
  }
  delete[] VT;
  delete[] Mnew;
}

//---------------------------------------------------------------------------------------------------//
// If you think these functions were obscure, wait for the QN conserving variant.
// They do essentially the same thing, but expand and decompose the matrices blockwise, thus preserving
// the block structure. In principle, therefore, the expansion term P is also modified. It is unknown
// if the truncated P is as useful as the heuristic ansatz by Hubig
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichmentBlockwise(int i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Bnew, *R, *BStart, *networkB;
  lapack_complex_double *pExpression;
  lapack_complex_double *U, *VT;
  double *diags;
  int blockDimR, blockDimL, maxDim;
  int containerDim;
  int lDRR, ldp;
  int siCurrent, aiCurrent, aimCurrent, aipCurrent;
  int numBlocks, lBlockSize, rBlockSize;
  getLocalDimensions(i);
  int const MNumCols=lDR*(1+lDwR);
  int const MNumRows=lDL*ld;
  lDRR=networkState.locDimR(i+1);
  ldp=locd(i+1);
  std::auto_ptr<lapack_complex_double> Rp(new lapack_complex_double[lDR*lDR*(1+lDwR)]);
  std::auto_ptr<lapack_complex_double> pEP(new lapack_complex_double[ld*lDL*lDR*lDwR]);
  std::auto_ptr<lapack_complex_double> MP, UP, VTP;
  std::auto_ptr<double> diagsP;
#ifndef USE_MKL
  std::auto_ptr<int> iworkP;
  std::auto_ptr<double> rworkP;
  std::auto_ptr<lapack_complex_double> workP;
  int *iwork;
  double *rwork;
  lapack_complex_double *work;
#endif
  pExpression=pEP.get();
  R=Rp.get();
  getPExpressionLeft(i,pExpression);
  numBlocks=networkState.indexTable.numBlocksLP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState.indexTable.lBlockSizeLP(i,iBlock);
    rBlockSize=networkState.indexTable.rBlockSizeLP(i,iBlock);
    blockDimL=lBlockSize;
    blockDimR=rBlockSize*(1+lDwR);
    maxDim=(blockDimL>blockDimR)?blockDimL:blockDimR;
    if(lBlockSize!=0 && rBlockSize!=0){
      MP.reset(new lapack_complex_double[blockDimL*blockDimR]);
      Mnew=MP.get();
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
      	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
	  siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
	  Mnew[k+j*lBlockSize]=networkState.global_access(i,siCurrent,aiCurrent,aimCurrent);
	  for(int bi=0;bi<lDwR;++bi){
	    Mnew[k+j*lBlockSize+(bi+1)*lBlockSize*rBlockSize]=alpha*pExpression[aimCurrent+lDL*siCurrent+bi*lDL*ld*lDR+aiCurrent*lDL*ld];
	  }
	}
      }
      containerDim=(lBlockSize>(rBlockSize*(1+lDwR)))?rBlockSize*(1+lDwR):lBlockSize;
      diagsP.reset(new double[containerDim]);
      UP.reset(new lapack_complex_double[blockDimL*blockDimL]);
      VTP.reset(new lapack_complex_double[blockDimR*blockDimR]);
      diags=diagsP.get();
      U=UP.get();
      VT=VTP.get();
      //There seems to be a bug in liblapacke providing a wrong size for the work array, which can lead to a segfault. This bug is not present in the mkl implementation
#ifdef USE_MKL
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR);
#endif
#ifndef USE_MKL
      int const lwork=containerDim*containerDim+maxDim*containerDim*2;
      int const lrwork=(containerDim*5+7>2*containerDim+2*maxDim+1)?containerDim*(containerDim*5+7):containerDim*(2*maxDim+2*containerDim+1);
      iworkP.reset(new int[8*containerDim]);
      rworkP.reset(new double[lrwork*4]);
      workP.reset(new lapack_complex_double[4*lwork]);
      work=workP.get();
      rwork=rworkP.get();
      iwork=iworkP.get();
      LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,'A',blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR,work,lwork,rwork,iwork);
#endif
      for(int mi=0;mi<blockDimR;++mi){
	for(int j=0;j<rBlockSize;++j){
	  //It should be blockDimR as leading dimension, since this is the structure of VT
	  VT[j+mi*blockDimR]*=diags[j];
	}
      }
      for(int b=0;b<(1+lDwR);++b){
	for(int jp=0;jp<rBlockSize;++jp){
	  aipCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,jp);
	  for(int j=0;j<rBlockSize;++j){
	    aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	    R[aiCurrent+lDR*aipCurrent+lDR*lDR*b]=VT[j+blockDimR*jp+blockDimR*rBlockSize*b];
	  }
	}
      }
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
	  siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
	  networkState.global_access(i,siCurrent,aiCurrent,aimCurrent)=U[k+lBlockSize*j];
	}
      }
    }
  }
  std::auto_ptr<lapack_complex_double> BP(new lapack_complex_double[ldp*lDRR*lDR*(1+lDwR)]);
  Bnew=BP.get();
  for(int si=0;si<ldp;++si){
    for(int air=0;air<lDRR;++air){
      for(int ai=0;ai<lDR;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=networkState.global_access(i+1,si,air,ai);
      }
      for(int ai=lDR;ai<MNumCols;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=0;
      }
    }
  }
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  for(int si=0;si<ldp;++si){
    networkState.subMatrixStart(networkB,i+1,si);
    BStart=Bnew+si*lDRR*MNumCols;
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,lDR*(1+lDwR),&zone,R,lDR,BStart,lDR*(1+lDwR),&zzero,networkB,lDR);
  }
}

//---------------------------------------------------------------------------------------------------//

void network::rightEnrichmentBlockwise(int i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Anew, *R, *networkA;
  lapack_complex_double *pExpression;
  lapack_complex_double *U, *VT;
  double *diags;
  int blockDimR, blockDimL, maxDim;
  int containerDim;
  int lDLL, ldm;
  int siCurrent, aiCurrent, aimCurrent, aimpCurrent;
  int numBlocks, lBlockSize, rBlockSize;
  getLocalDimensions(i);
  int const MNumCols=lDR*ld;
  int const MNumRows=lDL*(1+lDwL);
  lDLL=networkState.locDimL(i-1);
  ldm=locd(i-1);

  std::auto_ptr<lapack_complex_double> Rp(new lapack_complex_double[lDL*lDL*(1+lDwL)]);
  std::auto_ptr<lapack_complex_double> pEP(new lapack_complex_double[ld*lDL*lDR*lDwL]);
  std::auto_ptr<lapack_complex_double> MP, UP, VTP;
  std::auto_ptr<double> diagsP;
#ifndef USE_MKL
  std::auto_ptr<int> iworkP;
  std::auto_ptr<double> rworkP;
  std::auto_ptr<lapack_complex_double> workP;
  int *iwork;
  double *rwork;
  lapack_complex_double *work;
#endif
  
  R=Rp.get();
  pExpression=pEP.get();
  getPExpressionRight(i,pExpression);
  numBlocks=networkState.indexTable.numBlocksRP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState.indexTable.lBlockSizeRP(i,iBlock);
    rBlockSize=networkState.indexTable.rBlockSizeRP(i,iBlock);
    blockDimL=lBlockSize*(1+lDwL);
    blockDimR=rBlockSize;
    maxDim=(blockDimL>blockDimR)?blockDimL:blockDimR;
    if(lBlockSize!=0 && rBlockSize!=0){
      MP.reset(new lapack_complex_double[lBlockSize*rBlockSize*(1+lDwL)]);
      Mnew=MP.get();
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexRP(i,iBlock,j);
	siCurrent=networkState.indexTable.siBlockIndexRP(i,iBlock,j);
      	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexRP(i,iBlock,k);
	  Mnew[k+j*blockDimL]=networkState.global_access(i,siCurrent,aiCurrent,aimCurrent);
	  for(int bim=0;bim<lDwL;++bim){
	    Mnew[k+j*blockDimL+(bim+1)*lBlockSize]=alpha*pExpression[aimCurrent*lDwL+lDR*lDwL*lDL*siCurrent+bim+aiCurrent*lDL*lDwL];
	  }
	}
      }
      containerDim=(blockDimL>blockDimR)?blockDimR:blockDimL;
      diagsP.reset(new double[containerDim]);
      UP.reset(new lapack_complex_double[blockDimL*blockDimL]);
      VTP.reset(new lapack_complex_double[blockDimR*blockDimR]);
      diags=diagsP.get();
      U=UP.get();
      VT=VTP.get();
#ifdef USE_MKL
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR);
#endif
#ifndef USE_MKL
      int const lwork=containerDim*containerDim+maxDim*containerDim*2;
      int const lrwork=(containerDim*5+7>2*containerDim+2*maxDim+1)?containerDim*(containerDim*5+7):containerDim*(2*maxDim+2*containerDim+1);
      iworkP.reset(new int[8*containerDim]);
      rworkP.reset(new double[lrwork*4]);
      workP.reset(new lapack_complex_double[4*lwork]);
      work=workP.get();
      rwork=rworkP.get();
      iwork=iworkP.get();
      LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,'A',blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR,work,lwork,rwork,iwork);
#endif
      for(int mi=0;mi<blockDimL;++mi){
	for(int k=0;k<lBlockSize;++k){
	  //It should be blockDimL as leading dimension, since this is the structure of U
	  U[mi+k*blockDimL]*=diags[k];
	}
      }
      for(int b=0;b<(1+lDwL);++b){
	for(int kp=0;kp<lBlockSize;++kp){
	  aimpCurrent=networkState.indexTable.aimBlockIndexRP(i,iBlock,kp);
	  for(int k=0;k<lBlockSize;++k){
	    aimCurrent=networkState.indexTable.aimBlockIndexRP(i,iBlock,k);
	    R[aimCurrent+lDL*(1+lDwL)*aimpCurrent+lDL*b]=U[k+blockDimL*kp+lBlockSize*b];
	  }
	}
      }
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexRP(i,iBlock,j);
	siCurrent=networkState.indexTable.siBlockIndexRP(i,iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexRP(i,iBlock,k);
	  networkState.global_access(i,siCurrent,aiCurrent,aimCurrent)=VT[k+blockDimR*j];
	}
      }
    }
  }
  std::auto_ptr<lapack_complex_double> AP(new lapack_complex_double[ldm*lDLL*lDL*(1+lDwL)]);
  Anew=AP.get();
  for(int si=0;si<ldm;++si){
    for(int aimm=0;aimm<lDLL;++aimm){
      for(int aim=0;aim<lDL;++aim){
	Anew[aim*lDLL*ldm+aimm+si*lDLL]=networkState.global_access(i-1,si,aim,aimm);
      }
      for(int aim=lDL;aim<MNumRows;++aim){
	Anew[aim*ldm*lDLL+aimm+si*lDLL]=0;
      }
    }
  }
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  std::auto_ptr<lapack_complex_double> ASP(new lapack_complex_double[lDLL*MNumRows]);
  lapack_complex_double *AStart=ASP.get();
  for(int si=0;si<ldm;++si){
    networkState.subMatrixStart(networkA,i-1,si);
    for(int mi=0;mi<MNumRows;++mi){
      for(int aimm=0;aimm<lDLL;++aimm){
	AStart[aimm+mi*lDLL]=Anew[aimm+si*lDLL+mi*ld*lDLL];
      }
    }
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,MNumRows,&zone,AStart,lDLL,R,MNumRows,&zzero,networkA,lDLL);
  }
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary functions to determine the heuristic expansion term for the enrichment step.
//---------------------------------------------------------------------------------------------------//

void network::getPExpressionLeft(int const i, lapack_complex_double *pExpr){
  int const numBlocks=networkState.indexTable.numBlocksLP(i);
  int lBlockSize, rBlockSize;
  int aiCurrent, siCurrent, aimCurrent;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> outerContainer(ld,lDwR,lDR,lDL);
  //Use the measurement class to calculate the p-expression, which is one of the containers used in calculation of partial contractions. 
  pCtr.calcOuterContainerLeft(i+1,outerContainer);
  for(int m=0;m<lDL*lDR*ld*lDwR;++m){
    pExpr[m]=0;
  }
  //Copying has to be done explicitly for an optimized storage scheme of pExpr
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState.indexTable.lBlockSizeLP(i,iBlock);
    rBlockSize=networkState.indexTable.rBlockSizeLP(i,iBlock);
    for(int k=0;k<lBlockSize;++k){
      aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
      siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	for(int bi=0;bi<lDwR;++bi){
	  pExpr[aimCurrent+lDL*siCurrent+bi*lDL*ld*lDR+aiCurrent*lDL*ld]=outerContainer.global_access(siCurrent,bi,aiCurrent,aimCurrent);
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void network::getPExpressionRight(int i, lapack_complex_double *pExpr){
  //This works just like the left version, except for the shape of the used storage scheme (and other functions are called to get the right expression)
  int const numBlocks=networkState.indexTable.numBlocksRP(i);
  int lBlockSize, rBlockSize;
  int aiCurrent, siCurrent, aimCurrent;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> outerContainer(lDL,lDwL,ld,lDR);
  pCtr.calcOuterContainerRight(i-1,outerContainer);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState.indexTable.lBlockSizeRP(i,iBlock);
    rBlockSize=networkState.indexTable.rBlockSizeRP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiCurrent=networkState.indexTable.aiBlockIndexRP(i,iBlock,j);
      siCurrent=networkState.indexTable.siBlockIndexRP(i,iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimCurrent=networkState.indexTable.aimBlockIndexRP(i,iBlock,k);
	for(int bim=0;bim<lDwL;++bim){
	  pExpr[aimCurrent*lDwL+lDR*lDwL*lDL*siCurrent+bim+aiCurrent*lDL*lDwL]=outerContainer.global_access(aimCurrent,bim,siCurrent,aiCurrent);
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// Another auxiliary function to determine the current energy value for estimation of a new alpha
// This is quite effortive, so we currently use a guess for alpha and just decrease it continously 
// after each step
//---------------------------------------------------------------------------------------------------//

double network::getCurrentEnergy(int i){
  lapack_complex_double simpleContainer=0;
  lapack_complex_double *siteMatrixContainer=new lapack_complex_double [ld*lDR*lDL];
  lapack_complex_double *currentM, *LTerm, *RTerm;
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  networkState.subMatrixStart(currentM,i);
  //To compute the current energy, we multiply the current state with H and contract with the current state. This is done locally, the cached partial contractions are used for the rest of the network
  if(i==0 || i==L-1 || !pars.nQNs){
    optHMatrix gather(RTerm,LTerm,&networkH,networkDimInfo,networkH.maxDim(),i,0,0,0);
    gather.MultMv(currentM,siteMatrixContainer);
  }
  else{
    blockHMatrix BMat(RTerm,LTerm,&networkH,networkDimInfo,Dw,i,&(networkState.indexTable),0,0,&conservedQNs,1);
    BMat.prepareInput(currentM);
    BMat.MultMvBlockedLP(BMat.compressedVector,BMat.compressedVector);
    BMat.readOutput(siteMatrixContainer);
  }
  //Here, we contract the result with the current state. There is a version with and one without the use of QNs
  if(pars.nQNs){
    int const numBlocks=networkState.indexTable.numBlocksLP(i);
    int lBlockSize, rBlockSize;
    int siCurrent, aiCurrent, aimCurrent;
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=networkState.indexTable.lBlockSizeLP(i,iBlock);
      rBlockSize=networkState.indexTable.rBlockSizeLP(i,iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
	siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	  simpleContainer+=conj(currentM[aimCurrent+aiCurrent*lDL+siCurrent*lDR*lDL])*siteMatrixContainer[aimCurrent+aiCurrent*lDL+siCurrent*lDL*lDR];
	}
      }
    }
  }
  else{
    getLocalDimensions(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  simpleContainer+=conj(currentM[aim+ai*lDL+si*lDR*lDL])*siteMatrixContainer[aim+ai*lDL+si*lDR*lDL];
	}
      }
    }
  }
  delete[] siteMatrixContainer;
  return real(simpleContainer);
}

//---------------------------------------------------------------------------------------------------//

void network::getNewAlpha(int i, double lambda, double prevLambda){
  /*
  double dET=getCurrentEnergy(i);
  double dE0=abs(prevLambda-lambda);
  double const pre=1.2;
  dET-=lambda;
  std::cout<<dET<<" "<<dE0<<std::endl;
  if(dET<0){
    alpha*=pre;
  }
  else{
    if(abs(dET/dE0)<0.3){
      alpha*=pre;
    }
    else{
      alpha/=pre;
    }
  }
  */
  //This is numerically somewhat simpler and works also quite nice.
  //alpha*=0.9765;
  alpha*=pow(0.5,1/static_cast<double>(2*L));
  std::cout<<"New alpha="<<alpha<<std::endl;
}
