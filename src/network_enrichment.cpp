#include "mkl_complex_defined.h"
#include "arrayprocessing.h"
#include "tmpContainer.h"
#include "network.h"

//---------------------------------------------------------------------------------------------------//
// The implementation of the enrichment is very ugly and might still contain severe bugs which
// are irrelevant to the current testing system (Heisenberg chain). It does however, improve convergence
// speed and avoids local minima (at least thats what I expect, there are no problems with local 
// minima till now).
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichment(double const alpha, int const i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Bnew;
  lapack_complex_double *pExpression;
  int lDRR, ldp;
  getLocalDimensions(i);
  lDRR=networkState.locDimR(i+1);
  ldp=locd(i+1);
  int MNumCols=lDR*(1+lDwL);
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
  //Singular Value Decomposition of Mnew=U*S*V
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  double *diags=new double[containerDim];
  lapack_complex_double *Mnewcpy=new lapack_complex_double[maxDim*maxDim];
  //POSSIBLE IMPROVEMENT: USE LAPACK DRIVER ROUTINE INSTEAD OF COMPUTATIONAL ROUTINES
  lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int mi=0;mi<MNumCols;++mi){
    for(int ai=0;ai<lDR;++ai){
      Mnewcpy[ai+lDR*mi]*=diags[ai];
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
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,MNumCols,&zone,Mnewcpy,lDR,BStart,MNumCols,&zzero,networkB,lDR);
  }
  delete[] Mnewcpy;
  delete[] Bnew;
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=Mnew[aim+si*lDL+ai*MNumRows];
      }
    }
  }
  delete[] Mnew;
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
  lapack_complex_double *Mnewcpy=new lapack_complex_double[maxDim*maxDim];
  lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int aim=0;aim<lDL;++aim){
    for(int mi=0;mi<MNumRows;++mi){
      Mnew[mi+aim*MNumRows]*=diags[aim];
    }
  }
  delete[] diags;
  //From here, Mnew is to be treated as a MNumRows x lDL matrix
  //Postprocessing: A=Anew*U*S, B=U
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
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,MNumRows,&zone,AStart,lDLL,Mnew,MNumRows,&zzero,networkA,lDLL);
  }
  delete[] AStart;
  delete[] Mnew;
  delete[] Anew;
  //Store the remaining U into B
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=Mnewcpy[aim+si*lDR*MNumCols+ai*MNumCols];
      }
    }
  }
  delete[] Mnewcpy;
}

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
	  pExpr[aim+lDL*si+bi*lDL*d*lDR+ai*lDL*d]=simpleContainer;    
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

