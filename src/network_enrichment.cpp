#include "mkl_complex_defined.h"
#include "arrayprocessing.h"
#include "tmpContainer.h"
#include "network.h"
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// The implementation of the enrichment is very ugly and might still contain severe bugs which
// are irrelevant to the current testing system (Heisenberg chain). It does however, improve convergence
// speed and avoids local minima (at least thats what I expect, there are no problems with local 
// minima till now).
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichment(double const alpha, int const i){
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
  /*unity=new lapack_complex_double[MNumRows*MNumRows];
  for(int ai=0;ai<MNumRows;++ai){
    for(int aip=0;aip<MNumRows;++aip){
      unity[ai+aip*MNumRows]=(ai==aip)?-1:0;
    }
  }
  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,MNumRows,MNumRows,MNumRows,&zone,U,MNumRows,U,MNumRows,&zone,unity,MNumRows);
  std::cout<<"Verification of U: "<<cblas_dznrm2(MNumRows*MNumRows,unity,1)<<std::endl;
  delete[] unity;*/
  networkState.subMatrixStart(networkB,i);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
        networkState.global_access(i,si,ai,aim)=U[aim+si*lDL+ai*lDL*ld];
      }
    }
  }
  unity=new lapack_complex_double[lDR*lDR];
  for(int ai=0;ai<lDR;++ai){
    for(int aip=0;aip<lDR;++aip){
      unity[ai+aip*lDR]=(ai==aip)?-1:0;
    }
  }
  for(int si=0;si<ld;++si){
    transp(lDR,lDL,networkB+si*lDL*lDR);
  }
  cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,lDR,lDR,ld*lDL,&zone,networkB,lDR,networkB,lDR,&zone,unity,lDR);
  for(int si=0;si<ld;++si){
    transp(lDL,lDR,networkB+si*lDL*lDR);
  }
  std::cout<<"Verification: "<<cblas_dznrm2(lDR*lDR,unity,1)<<std::endl;
  delete[] unity;
  delete[] Mnew;
  delete[] U;
}

//---------------------------------------------------------------------------------------------------//

void network::rightEnrichment(double const alpha, int const i){
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
  //Store the remaining U into B
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=VT[aim+si*lDR*MNumCols+ai*MNumCols];
      }
    }
  }
  lapack_complex_double *unity=new lapack_complex_double[lDL*lDL];
  for(int ai=0;ai<lDL;++ai){
    for(int aip=0;aip<lDL;++aip){
      unity[ai+aip*lDL]=(ai==aip)?-1:0;
    }
  }
  networkState.subMatrixStart(networkA,i);
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,lDL,lDL,ld*lDR,&zone,networkA,lDL,networkA,lDL,&zone,unity,lDL);
  std::cout<<"Verification: "<<cblas_dznrm2(lDL*lDL,unity,1)<<std::endl;
  delete[] unity;
  delete[] VT;
  delete[] Mnew;
}

//---------------------------------------------------------------------------------------------------//

void network::leftEnrichmentBlockwise(double const alpha, int const i){
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
  R=new lapack_complex_double[lDR*lDR*(1+lDwR)];
  pExpression=new lapack_complex_double[ld*lDL*lDR*lDwR];
  getPExpressionLeft(i,pExpression);
  numBlocks=networkState.indexTable.numBlocksLP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState.indexTable.lBlockSizeLP(i,iBlock);
    rBlockSize=networkState.indexTable.rBlockSizeLP(i,iBlock);
    blockDimL=lBlockSize;
    blockDimR=rBlockSize*(1+lDwR);
    maxDim=(blockDimL>blockDimR)?blockDimL:blockDimR;
    if(lBlockSize!=0 && rBlockSize!=0){
      Mnew=new lapack_complex_double[lBlockSize*rBlockSize*(1+lDwR)];
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
      	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
	  siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
	  Mnew[k+j*lBlockSize]=networkState.global_access(i,siCurrent,aiCurrent,aimCurrent);
	  for(int bi=0;bi<lDwR;++bi){
	    Mnew[k+j*lBlockSize+(bi+1)*lBlockSize*rBlockSize]=alpha*pExpression[aimCurrent+lDL*siCurrent+bi*lDL*ld*lDR+aiCurrent*lDL*ld];
	    // This is some kind of voodoo
	    if(Mnew[k+j*lBlockSize+(bi+1)*lBlockSize*rBlockSize]==1e100){
	      //std::cout<<"Debug\n";
	    }
	  }
	}
      }
      containerDim=(lBlockSize>(rBlockSize*(1+lDwR)))?rBlockSize*(1+lDwR):lBlockSize;
      diags=new double[containerDim];
      U=new lapack_complex_double[blockDimL*blockDimL];
      VT=new lapack_complex_double[blockDimR*blockDimR];
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR);
      delete[] Mnew;
      for(int mi=0;mi<blockDimR;++mi){
	for(int j=0;j<rBlockSize;++j){
	  //It should be blockDimR as leading dimension, since this is the structure of VT
	  VT[j+mi*blockDimR]*=diags[j];
	}
      }
      delete[] diags;
      for(int b=0;b<(1+lDwR);++b){
	for(int jp=0;jp<rBlockSize;++jp){
	  aipCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,jp);
	  for(int j=0;j<rBlockSize;++j){
	    aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	    R[aiCurrent+lDR*aipCurrent+lDR*lDR*b]=VT[j+blockDimR*jp+blockDimR*rBlockSize*b];
	  }
	}
      }
      delete[] VT;
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=networkState.indexTable.aiBlockIndexLP(i,iBlock,j);
	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=networkState.indexTable.aimBlockIndexLP(i,iBlock,k);
	  siCurrent=networkState.indexTable.siBlockIndexLP(i,iBlock,k);
	  networkState.global_access(i,siCurrent,aiCurrent,aimCurrent)=U[k+lBlockSize*j];
	}
      }
      delete[] U;
    }
  }
  Bnew=new lapack_complex_double[ldp*lDRR*lDR*(1+lDwR)];
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
  delete[] Bnew;
  delete[] pExpression;
}

//---------------------------------------------------------------------------------------------------//



//---------------------------------------------------------------------------------------------------//
// Auxiliary functions to determine the heuristic expansion term for the enrichment step.
//---------------------------------------------------------------------------------------------------//

void network::getPExpressionLeft(int const i, lapack_complex_double *pExpr){
  lapack_complex_double simpleContainer;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> innerContainer(ld,lDL,lDR,lDwL);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDL;++aim){
      for(int ai=0;ai<lDR;++ai){
	for(int bim=0;bim<lDwL;++bim){
	  simpleContainer=0;
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=pCtr.Lctr.global_access(i,aim,bim,aimp)*networkState.global_access(i,si,ai,aimp);
	  }
	  innerContainer.global_access(si,aim,ai,bim)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int bi=0;bi<lDwR;++bi){
	for(int aim=0;aim<lDL;++aim){
	  //BEWARE: IM UNSURE WHETHER bi OR ai IS THE OUTERMOST INDEX - FORGET IT, ITS bi
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bim=0;bim<lDwL;++bim){
	      simpleContainer+=networkH.global_access(i,si,sip,bi,bim)*innerContainer.global_access(sip,aim,ai,bim);
	    }
	  }
	  pExpr[aim+lDL*si+bi*lDL*ld*lDR+ai*lDL*ld]=simpleContainer;    
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void network::getPExpressionRight(int const i, lapack_complex_double *pExpr){
  lapack_complex_double simpleContainer;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> innerContainer(ld,lDL,lDR,lDwR);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDL;++aim){
      for(int ai=0;ai<lDR;++ai){
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=pCtr.Rctr.global_access(i,ai,bi,aip)*networkState.global_access(i,si,aip,aim);
	  }
	  innerContainer.global_access(si,aim,ai,bi)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int bim=0;bim<lDwL;++bim){
	for(int aim=0;aim<lDL;++aim){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=networkH.global_access(i,si,sip,bi,bim)*innerContainer.global_access(sip,aim,ai,bi);
	    }
	  }
	  pExpr[aim*lDwL+bim+ai*lDL*lDwL+si*lDL*lDwL*lDR]=simpleContainer;
	}
      }
    }
  }
}

