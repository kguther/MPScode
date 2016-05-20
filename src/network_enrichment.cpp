#include "mkl_complex_defined.h"
#include "arrayprocessing.h"
#include "tmpContainer.h"
#include "network.h"
#include "optHMatrix.h"
#include "blockHMatrix.h"
#include "truncation.h"
#include "verifyQN.h"
#include "exceptionClasses.h"
#include <cmath>
#include <algorithm>
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
  lapack_complex_double *pExpression, *unity, *localMatrix, *localBMatrix;
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
  networkState.subMatrixStart(localMatrix,i);
  networkState.subMatrixStart(localBMatrix,i+1);
  for(int si=0;si<ld;++si){
    //Copy A and B while rearranging to allow for expansion
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+si*lDL+ai*MNumRows]=localMatrix[stateIndex(si,ai,aim)];
      }
    }
    for(int air=0;air<lDRR;++air){
      for(int ai=0;ai<lDR;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=localBMatrix[si*lDR*lDRR+air*lDR+ai];
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
  networkState.subMatrixStart(localMatrix,i);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
        localMatrix[stateIndex(si,ai,aim)]=U[aim+si*lDL+ai*lDL*ld];
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
  lapack_complex_double *pExpression, *localMatrix, *localAMatrix;
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
  networkState.subMatrixStart(localMatrix,i);
  networkState.subMatrixStart(localAMatrix,i-1);
  //Copy A and B while rearranging to allow for expansion
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+ai*MNumRows+si*lDR*MNumRows]=localMatrix[stateIndex(si,ai,aim)];
      }
    }
    for(int aim=0;aim<lDL;++aim){
      for(int aimm=0;aimm<lDLL;++aimm){
	Anew[aimm+aim*ld*lDLL+si*lDLL]=localAMatrix[si*lDL*lDLL+aim*lDLL+aimm];
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
  networkState.subMatrixStart(localMatrix,i);
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
        localMatrix[stateIndex(si,ai,aim)]=VT[aim+si*lDR*MNumCols+ai*MNumCols];
      }
    }
  }
  delete[] VT;
  delete[] Mnew;
}

//---------------------------------------------------------------------------------------------------//
// If you think these functions were obscure, wait for the QN conserving variant.
// They do essentially the same thing, but expand and decompose the matrices blockwise, thus preserving
// the block structure.

// Also, the QN labeling scheme is dynamically updated during enrichment
// In particular, each block gets expanded on its own, this seems to be the only 
// way to guarantee a well defined labeling scheme is generated.

// This is the most complicated part of the program
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichmentBlockwise(int i){
  lapack_complex_double *Bnew, *R, *BStart, *localMatrix;
  lapack_complex_double *pExpression;
  lapack_complex_double *globalU, *globalVT;
  int blockDimR, blockDimL;
  int containerDim;
  int lDRR, ldp;
  int siCurrent, aiCurrent, aimCurrent, aipCurrent, sipCurrent;
  int lBlockSize, rBlockSize;


  //siteQNOrderMatrix const localIndexTable=networkState.indexTable().getLocalIndexTable(i);
  siteQNOrderMatrix localIndexTableFull;
  localIndexTableFull.generateFull(networkState.indexTable().getLocalIndexTable(i));

  getLocalDimensions(i);
  //pSize is the right dimension of the expansion term
  int const pSize=lDwR*lDR;
  int const numBlocks=localIndexTableFull.numBlocksLP();
  int const MNumCols=lDR+numBlocks*pSize;
  int const MNumRows=lDL*ld;
  lDRR=networkState.locDimR(i+1);
  ldp=locd(i+1);
  std::unique_ptr<lapack_complex_double[]> Rp(new lapack_complex_double[lDR*lDR]);
  std::unique_ptr<lapack_complex_double[]> pEP(new lapack_complex_double[ld*lDL*lDR*lDwR]);
  pExpression=pEP.get();
  R=Rp.get();
  
  //Prepare global containers
  std::unique_ptr<lapack_complex_double[]> globalUP(new lapack_complex_double[MNumRows*MNumRows]);
  //Formally, globalVT is of size MNumColsxMNumCols. But we never use this due to truncation
  //Save memory, use only MNumCols*lDR
  std::unique_ptr<lapack_complex_double[]> globalVTP(new lapack_complex_double[MNumCols*lDR]);
  std::vector<auxiliary::sortData> comparer(MNumRows);
  globalU=globalUP.get();
  globalVT=globalVTP.get();
  //Initialize with 0:
  for(int m=0;m<MNumRows*MNumRows;++m){
    globalU[m]=0.0;
  }
  for(int m=0;m<MNumCols*lDR;++m){
    globalVT[m]=0.0;
  }
  //Also diags, although not strictly necessary
  for(int m=0;m<comparer.size();++m){
    comparer[m].lambda=0.0;
  }

  //Buffer for indices
  //For parallelization, all rowOffsets must be known in advance
  std::unique_ptr<int[]> rowOffsetVT(new int[numBlocks]);
  rowOffsetVT[0]=0;
  for(int iBlock=0;iBlock<numBlocks-1;++iBlock){
    lBlockSize=localIndexTableFull.lBlockSizeLP(iBlock);
    rBlockSize=localIndexTableFull.rBlockSizeLP(iBlock);
    containerDim=(lBlockSize>(rBlockSize+pSize))?(rBlockSize+pSize):lBlockSize;
    rowOffsetVT[iBlock+1]=rowOffsetVT[iBlock]+containerDim;
  }


  getPExpressionLeft(i,pExpression);
  networkState.subMatrixStart(localMatrix,i);

#pragma omp parallel for private(lBlockSize, rBlockSize, blockDimL, blockDimR, containerDim, siCurrent, aiCurrent, aimCurrent, aipCurrent, sipCurrent) schedule(dynamic,1) 
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=localIndexTableFull.lBlockSizeLP(iBlock);
    rBlockSize=localIndexTableFull.rBlockSizeLP(iBlock);
    blockDimL=lBlockSize;
    blockDimR=rBlockSize+pSize;
    //rBlockSize may be zero, then, only the expansion term is taken into account
    if(lBlockSize!=0){
      //Mnew has to be private for each thread
      lapack_complex_double *Mnew;
      std::unique_ptr<lapack_complex_double[]> MP(new lapack_complex_double[blockDimL*blockDimR]);
      Mnew=MP.get();
      for(int k=0;k<lBlockSize;++k){
	aimCurrent=localIndexTableFull.aimBlockIndexLP(iBlock,k);
	siCurrent=localIndexTableFull.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiCurrent=localIndexTableFull.aiBlockIndexLP(iBlock,j);
	  Mnew[k+j*lBlockSize]=localMatrix[stateIndex(siCurrent,aiCurrent,aimCurrent)];
	}
	//The expansion term only has a charge associated with (si,aim) but not with ai
	//Each block own a copy of pExpression
	for(int ai=0;ai<lDR;++ai){
	  for(int bi=0;bi<lDwR;++bi){
	    Mnew[rBlockSize*lBlockSize+k+lBlockSize*(ai+lDR*bi)]=alpha*pExpression[aimCurrent+lDL*siCurrent+bi*lDL*ld*lDR+ai*lDL*ld];
	  }
	}
      }
      containerDim=(blockDimL>blockDimR)?blockDimR:blockDimL;
      //Allocate containers for lapack 
      std::unique_ptr<double[]> diagsP(new double[containerDim]);
      std::unique_ptr<lapack_complex_double[]> UP(new lapack_complex_double[blockDimL*blockDimL]);
      std::unique_ptr<lapack_complex_double[] >VTP(new lapack_complex_double[blockDimR*blockDimR]);
      double *diags=diagsP.get();
      lapack_complex_double *U=UP.get();
      lapack_complex_double *VT=VTP.get();

      //We only need the first containerDim rows of VT, but all cols of U 
      //-> U and VT have to be computed full in the rare case that there are more left singular vectors than right singular vectors (this CAN occur!) 
      //else, the first blockDimL right singular vectors are enough
      char jobz=(blockDimL>blockDimR)?'A':'S';
      //There seems to be a bug in liblapacke providing a wrong size for the work array, which can lead to a segfault. This bug is not present in the mkl implementation
#ifdef USE_MKL
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,jobz,blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR);
#endif
#ifndef USE_MKL
      int const maxDim=(blockDimR>blockDimL)?blockDimR:blockDimL;
      int const lwork=containerDim*containerDim+maxDim*containerDim*2;
      int const lrwork=(containerDim*5+7>2*containerDim+2*maxDim+1)?containerDim*(containerDim*5+7):containerDim*(2*maxDim+2*containerDim+1);
      std::unique_ptr<int[]> iworkP(new int[8*containerDim]);
      std::unique_ptr<double[]> rworkP(new double[lrwork*4]);
      std::unique_ptr<lapack_complex_double[]> workP(new lapack_complex_double[4*lwork]);
      lapack_complex_double *work=workP.get();
      double *rwork=rworkP.get();
      int *iwork=iworkP.get();
      LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,jobz,blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR,work,lwork,rwork,iwork);
#endif
      //Insert the block matrices into the global ones at the respective position

      //Insert VT into globalVT (b does not carry a qn label)
      std::unique_ptr<int[]> globalVTIndices(new int[containerDim]);
      for(int j=0;j<containerDim;++j){
	globalVTIndices[j]=rowOffsetVT[iBlock]+j;
      }
      //Only the first containerDim rows of VT are used (only that many SVs are nonzero)
      for(int jp=0;jp<rBlockSize;++jp){
	for(int j=0;j<containerDim;++j){
	  //Only the first lDL cols of VT are used, the rest is multiplied with 0
	  //Each block has its own designated rows in globalVT, leading to a blockstructure:
	  aiCurrent=localIndexTableFull.aiBlockIndexLP(iBlock,jp);
	  globalVT[globalVTIndices[j]+MNumCols*aiCurrent]=VT[j+blockDimR*jp];
	}
      }
      
      //Insert U into globalU (storage scheme for globalU: (sip,aip,si,ai)
      for(int kp=0;kp<blockDimL;++kp){
	//Inserting U is straighforward as the index table can be used
	aiCurrent=localIndexTableFull.aimBlockIndexLP(iBlock,kp);
	sipCurrent=localIndexTableFull.siBlockIndexLP(iBlock,kp);
	for(int k=0;k<blockDimL;++k){
	  aimCurrent=localIndexTableFull.aimBlockIndexLP(iBlock,k);
	  siCurrent=localIndexTableFull.siBlockIndexLP(iBlock,k);
	  globalU[aimCurrent+lDL*siCurrent+lDL*ld*aiCurrent+lDL*lDL*ld*sipCurrent]=U[k+blockDimL*kp];
	}
      }
     
      //Insert diags into comparer (for label-tracking truncation)
      auxiliary::sortData sortDataBuf;
      for(int j=0;j<containerDim;++j){
	//containerDim<=lBlockSize -> no conversion needed
	aimCurrent=localIndexTableFull.aimBlockIndexLP(iBlock,j);
	siCurrent=localIndexTableFull.siBlockIndexLP(iBlock,j);
	sortDataBuf.lambda=diags[j];
	//These will be the new labels. Because we need to be able to generate new labels, a full index table is required here
        sortDataBuf.QN=localIndexTableFull.qnLabelLP(iBlock);
	sortDataBuf.index=aimCurrent+lDL*siCurrent;
        sortDataBuf.indexExp=globalVTIndices[j];
	comparer[aimCurrent+lDL*siCurrent]=sortDataBuf;
      }
    }
  }

  //For testing
  quantumNumber gqn=networkState.getConservedQNs()[0];

  //Get the columns of globalU corresponding to the lDR highest SVs and copy them to A
  std::sort(comparer.begin(),comparer.end(),auxiliary::compareSortData);
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(localMatrix,i,si);
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	localMatrix[aim+lDL*ai]=globalU[aim+lDL*si+lDL*ld*comparer[ai].index];
      }
    }
  }

  //Construct the matrix R=SVT (truncated, lDRxlDR*(1+lDwR))
  for(int m=0;m<lDR;++m){
    for(int ai=0;ai<lDR;++ai){
      R[ai+lDR*m]=globalVT[comparer[ai].indexExp+MNumCols*m]*comparer[ai].lambda;

      //if(abs(R[ai+lDR*m])>1e-12 && m<lDR)
      //std::cout<<"Entry at "<<gqn.QNLabel(i,ai)<<"->"<<comparer[ai].QN<<std::endl;
    }
  }

  //Copy B for multiplication (virtually add zeros, but they do not contribute, so leave them in the container)
  std::unique_ptr<lapack_complex_double[]> BP(new lapack_complex_double[ldp*lDRR*lDR]);
  Bnew=BP.get();
  networkState.subMatrixStart(localMatrix,i+1);
  for(int si=0;si<ldp;++si){
    for(int air=0;air<lDRR;++air){
      for(int ai=0;ai<lDR;++ai){
	Bnew[ai+air*lDR+si*lDRR*lDR]=localMatrix[si*lDRR*lDR+air*lDR+ai];
      }
    }
  }
  
  //Multiply R into B to get a guess for the next matrix
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  for(int si=0;si<ldp;++si){
    networkState.subMatrixStart(localMatrix,i+1,si);
    BStart=Bnew+si*lDRR*lDR;
    //Formally, R is of size lDRxMNumCols and B of size MNumColsxlDRR but only the first lDR rows of B are nonzero
    //Therefore, only calculate the matrix-matrix product for the first lDR cols
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,lDR,&zone,R,lDR,BStart,lDR,&zzero,localMatrix,lDR);
  }
 
  //Refine QN
  std::vector<std::complex<int> > newLabels(lDR);
  for(int ai=0;ai<lDR;++ai){
    newLabels[ai]=comparer[ai].QN;
  }
  try{
    networkState.refineQNLabels(i+1,0,newLabels);
  }
  catch(empty_table &err){
    networkState.adaptLabels(err.site(),1);
  }
  std::cout<<"Applied enrichment\n";
  //Works nicely, A and B keep the QNC - testing is not required in general
  /*
  info=0;

  info=checkQNConstraint(networkState,i);
  if(info)
    std::cout<<"Matrix A violates QNC\n";
  info=checkQNConstraint(networkState,i+1);
  if(info)
    std::cout<<"Matrix B violates QNC\n";

  for(int i=0;i<L;++i){
    info+=checkQNConstraint(networkState,i);
  }
  if(info){
    std::cout<<"QNC VIOLATION DETECTED\n";
    exit(1);
  }
  */
}

//---------------------------------------------------------------------------------------------------//

void network::rightEnrichmentBlockwise(int i){
  lapack_complex_double *Anew, *R, *networkA, *localMatrix;
  lapack_complex_double *pExpression;
  lapack_complex_double *globalU, *globalVT;
  int blockDimR, blockDimL;
  int containerDim;
  int lDLL, ldm;
  int siCurrent, aiCurrent, aimCurrent, aimpCurrent, sipCurrent, aipCurrent;
  int lBlockSize, rBlockSize;
  getLocalDimensions(i);

  //Full index table, derived using the reachable indices (contains empty blocks, therefore slower and not used in the other components)
  siteQNOrderMatrix localIndexTableFull;
  localIndexTableFull.generateFull(networkState.indexTable().getLocalIndexTable(i));

  //pSize is the left dimension of the expansion term
  int const pSize=lDL*lDwL;
  int const numBlocks=localIndexTableFull.numBlocksRP();
  int const MNumCols=lDR*ld;
  int const MNumRows=lDL+numBlocks*pSize;
  lDLL=networkState.locDimL(i-1);
  ldm=locd(i-1);

  //MPS index table
  //siteQNOrderMatrix const localIndexTable=networkState.indexTable().getLocalIndexTable(i);

  //Prepare global containers
  std::unique_ptr<lapack_complex_double[]> Rp(new lapack_complex_double[lDL*lDL]);
  std::unique_ptr<lapack_complex_double[]> globalUP(new lapack_complex_double[MNumRows*lDL]);
  std::unique_ptr<lapack_complex_double[]> globalVTP(new lapack_complex_double[MNumCols*MNumCols]);
  //comparer has to have at least size min(MNumCols,MNumRows)
  std::vector<auxiliary::sortData> comparer(MNumCols);
  globalU=globalUP.get();
  globalVT=globalVTP.get();
  //Initialize with zero - this is important since not all entries of VT and U are set
  for(int m=0;m<MNumCols*MNumCols;++m){
    globalVT[m]=0.0;
  }
  for(int m=0;m<MNumRows*lDL;++m){
    globalU[m]=0.0;
  }
  for(int m=0;m<MNumCols;++m){
    comparer[m].lambda=0.0;
  }
  //Buffer for indices
  std::unique_ptr<int[]> colOffsetU(new int[numBlocks]);
  colOffsetU[0]=0;
  //For parallelization, all colOffsets have to be known in advance
  for(int iBlock=0;iBlock<numBlocks-1;++iBlock){
    lBlockSize=localIndexTableFull.lBlockSizeRP(iBlock);
    rBlockSize=localIndexTableFull.rBlockSizeRP(iBlock);
    containerDim=(rBlockSize>(lBlockSize+pSize))?(lBlockSize+pSize):rBlockSize;
    colOffsetU[iBlock+1]=colOffsetU[iBlock]+containerDim;
  }
    
  //Prepare local containers
  std::unique_ptr<lapack_complex_double[]> pEP(new lapack_complex_double[ld*lDL*lDR*lDwL]);
  R=Rp.get();
  pExpression=pEP.get();
  getPExpressionRight(i,pExpression);
  networkState.subMatrixStart(localMatrix,i);

#pragma omp parallel for private(lBlockSize, rBlockSize, blockDimL, blockDimR, containerDim, siCurrent, aiCurrent, aimCurrent, aimpCurrent, sipCurrent, aipCurrent) schedule(dynamic,1)
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=localIndexTableFull.lBlockSizeRP(iBlock);
    rBlockSize=localIndexTableFull.rBlockSizeRP(iBlock);
    blockDimL=lBlockSize+pSize;
    blockDimR=rBlockSize;
    //By construction of localIndexTableFull, rBlockSize can not be 0 (that is, such blocks are not listed in localIndexTableFull)
    if(rBlockSize!=0){
      std::unique_ptr<lapack_complex_double[]> MP(new lapack_complex_double[blockDimL*blockDimR]);
      lapack_complex_double *Mnew=MP.get();
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=localIndexTableFull.aiBlockIndexRP(iBlock,j);
	siCurrent=localIndexTableFull.siBlockIndexRP(iBlock,j);
      	for(int k=0;k<lBlockSize;++k){
	  aimCurrent=localIndexTableFull.aimBlockIndexRP(iBlock,k);
	  Mnew[k+j*blockDimL]=localMatrix[stateIndex(siCurrent,aiCurrent,aimCurrent)];
	}
	for(int aim=0;aim<lDL;++aim){
	  for(int bim=0;bim<lDwL;++bim){
	    Mnew[lBlockSize+(aim+lDL*bim)+j*blockDimL]=alpha*pExpression[aim*lDwL+lDR*lDwL*lDL*siCurrent+bim+aiCurrent*lDL*lDwL];
	  }
	}
      }
      //containerDim is the number of SVs
      containerDim=(blockDimL>blockDimR)?blockDimR:blockDimL;
      std::unique_ptr<double[]> diagsP(new double[containerDim]);
      std::unique_ptr<lapack_complex_double[]> UP(new lapack_complex_double[blockDimL*blockDimL]);
      std::unique_ptr<lapack_complex_double[]> VTP(new lapack_complex_double[blockDimR*blockDimR]);
      double *diags=diagsP.get();
      lapack_complex_double *U=UP.get();
      lapack_complex_double *VT=VTP.get();

      //Same simplification as in leftEnrichment
      char jobz=(blockDimR>blockDimL)?'A':'S';
#ifdef USE_MKL
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,jobz,blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR);
#endif
#ifndef USE_MKL
      int const maxDim=(blockDimR>blockDimL)?blockDimR:blockDimL;
      int const lwork=containerDim*containerDim+maxDim*containerDim*2;
      int const lrwork=(containerDim*5+7>2*containerDim+2*maxDim+1)?containerDim*(containerDim*5+7):containerDim*(2*maxDim+2*containerDim+1);
      std::unique_ptr<int[]> iworkP(new int[8*containerDim]);
      std::unique_ptr<double[]> rworkP(new double[lrwork*4]);
      std::unique_ptr<lapack_complex_double[]> workP(new lapack_complex_double[4*lwork]);
      lapack_complex_double *work=workP.get();
      double *rwork=rworkP.get();
      int *iwork=iworkP.get();
      LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,jobz,blockDimL,blockDimR,Mnew,blockDimL,diags,U,blockDimL,VT,blockDimR,work,lwork,rwork,iwork);
#endif

      //Insert VT into globalVT (storage scheme for VT: (sip,aip,si,ai)
      for(int j=0;j<blockDimR;++j){
	aiCurrent=localIndexTableFull.aiBlockIndexRP(iBlock,j);
	siCurrent=localIndexTableFull.siBlockIndexRP(iBlock,j);
	for(int jp=0;jp<blockDimR;++jp){
	  aipCurrent=localIndexTableFull.aiBlockIndexRP(iBlock,jp);
	  sipCurrent=localIndexTableFull.siBlockIndexRP(iBlock,jp);
	  globalVT[aipCurrent+lDR*sipCurrent+ld*lDR*(aiCurrent+lDR*siCurrent)]=VT[jp+blockDimR*j];
	}
      }
      
      //Insert U into globalU
      //The right indices of U can be given any order
      std::unique_ptr<int[]> globalUIndices(new int[containerDim]);
      for(int k=0;k<containerDim;++k){
	globalUIndices[k]=colOffsetU[iBlock]+k;
      }
      //Only the first containerDim cols of U are used
      for(int k=0;k<containerDim;++k){
	for(int kp=0;kp<lBlockSize;++kp){
	  //The left indices need to respect the order of A in the first lDL rows of U
	  aimCurrent=localIndexTableFull.aimBlockIndexRP(iBlock,kp);
	  globalU[aimCurrent+lDL*globalUIndices[k]]=U[kp+blockDimL*k];
	}
      }
      
      //Insert diags into comparer
      auxiliary::sortData sortDataBuf;
      for(int k=0;k<containerDim;++k){
	//containerDim<=rBlockSize -> no conversion needed
	aiCurrent=localIndexTableFull.aiBlockIndexRP(iBlock,k);
	siCurrent=localIndexTableFull.siBlockIndexRP(iBlock,k);
	sortDataBuf.lambda=diags[k];
	sortDataBuf.QN=localIndexTableFull.qnLabelRP(iBlock);
	sortDataBuf.index=aiCurrent+lDR*siCurrent;
	sortDataBuf.indexExp=globalUIndices[k];
	comparer[aiCurrent+lDR*siCurrent]=sortDataBuf;
      }
    }
  }
  //For testing
  quantumNumber gqn=networkState.getConservedQNs()[0];

  //Sort comparer, only the first lDL elements are used from now on
  std::sort(comparer.begin(),comparer.end(),auxiliary::compareSortData);

  //Get the rows of globalVT corresponding to the lDL highest SVs and copy them to B
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(localMatrix,i,si);
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	localMatrix[aim+lDL*ai]=globalVT[comparer[aim].index+ld*lDR*(ai+lDR*si)];

	//if(abs(localMatrix[aim+lDL*ai])>1e-10)
	//std::cout<<"Entry "<<localMatrix[aim+lDL*ai]<<" at "<<comparer[aim].QN<<"+"<<gqn.QNLabel(si)<<"="<<gqn.QNLabel(i,ai)<<std::endl;
      }
    }
  }
  
  //Construct the matrix R=US (truncated, lDL*(lDL+lDL*lDwL*numBlocks))
  //Only the first lDL rows are used, do not even allocate the rest
  for(int m=0;m<lDL;++m){
    for(int aimp=0;aimp<lDL;++aimp){
      R[m+lDL*aimp]=globalU[m+lDL*comparer[aimp].indexExp]*comparer[aimp].lambda;
      
      //if(abs(R[m+MNumRows*aimp])>1e-15)
	//std::cout<<"Entry at "<<gqn.QNLabel(i-1,m)<<"->"<<comparer[aimp].QN<<" with indices "<<m<<"->"<<aimp<<"("<<comparer[aimp].indexExp<<")"<<std::endl;
    }
  }    
  
  //Add zeros to the A-Matrix (using an external container)
  std::unique_ptr<lapack_complex_double[]> AP(new lapack_complex_double[ldm*lDLL*lDL]);
  Anew=AP.get();
  networkState.subMatrixStart(localMatrix,i-1);
  for(int si=0;si<ldm;++si){
    for(int aimm=0;aimm<lDLL;++aimm){
      for(int aim=0;aim<lDL;++aim){
	Anew[aimm+lDLL*aim+lDLL*lDL*si]=localMatrix[si*lDLL*lDL+aim*lDLL+aimm];
      }
    }
  }

  //Multiply R into A to get a guess for the next matrix
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  std::unique_ptr<lapack_complex_double[]> ASP(new lapack_complex_double[lDLL*lDL]);
  lapack_complex_double *AStart=ASP.get();
  for(int si=0;si<ldm;++si){
    networkState.subMatrixStart(networkA,i-1,si);
    //Anew is just a copy of A, take the si-submatrix
    AStart=Anew+lDLL*lDL*si;
    //And multiplied with R to give the si-component of the next guess
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,lDL,&zone,AStart,lDLL,R,lDL,&zzero,networkA,lDLL);
  }

  //Refine QN labels
  std::vector<std::complex<int> > newLabels(lDL);
  for(int aim=0;aim<lDL;++aim){
    newLabels[aim]=comparer[aim].QN;
  }
  try{
    networkState.refineQNLabels(i,0,newLabels);
  }
  catch(empty_table &err){
    networkState.adaptLabels(err.site(),-1);
  }
  std::cout<<"Applied enrichment\n";
  //Check for QNC conservation has been removed for better performance
  /*
  int info=0;
  info=checkQNConstraint(networkState,i-1);
  if(info)
    std::cout<<"Matrix A violates QNC\n";
  info=checkQNConstraint(networkState,i);
  if(info)
    std::cout<<"Matrix B violates QNC\n";

  for(int i=0;i<L;++i){
    info+=checkQNConstraint(networkState,i);
  }
  if(info){
    std::cout<<"QNC VIOLATION DETECTED\n";
    exit(1);
  }
  */
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary functions to determine the heuristic expansion term for the enrichment step.
//---------------------------------------------------------------------------------------------------//

void network::getPExpressionLeft(int const i, lapack_complex_double *pExpr){
  int lBlockSize, rBlockSize;
  int aiCurrent, siCurrent, aimCurrent;
  siteQNOrderMatrix const localIndexTable=networkState.indexTable().getLocalIndexTable(i);
  int const numBlocks=localIndexTable.numBlocksLP();
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> outerContainer(ld,lDwR,lDR,lDL);
  //Use the measurement class to calculate the p-expression, which is one of the containers used in calculation of partial contractions. 
  pCtr.calcOuterContainerLeft(i+1,outerContainer);
  for(int m=0;m<lDL*lDR*ld*lDwR;++m){
    pExpr[m]=0;
  }
  //Copying has to be done explicitly for an optimized storage scheme of pExpr
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
    rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
    for(int k=0;k<lBlockSize;++k){
      aimCurrent=localIndexTable.aimBlockIndexLP(iBlock,k);
      siCurrent=localIndexTable.siBlockIndexLP(iBlock,k);
      for(int j=0;j<rBlockSize;++j){
	aiCurrent=localIndexTable.aiBlockIndexLP(iBlock,j);
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
  siteQNOrderMatrix const localIndexTable=networkState.indexTable().getLocalIndexTable(i);
  int const numBlocks=localIndexTable.numBlocksRP();
  int lBlockSize, rBlockSize;
  int aiCurrent, siCurrent, aimCurrent;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> outerContainer(lDL,lDwL,ld,lDR);
  pCtr.calcOuterContainerRight(i-1,outerContainer);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=localIndexTable.lBlockSizeRP(iBlock);
    rBlockSize=localIndexTable.rBlockSizeRP(iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiCurrent=localIndexTable.aiBlockIndexRP(iBlock,j);
      siCurrent=localIndexTable.siBlockIndexRP(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimCurrent=localIndexTable.aimBlockIndexRP(iBlock,k);
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
  lapack_complex_double simpleContainer=0.0;
  getLocalDimensions(i);
  std::unique_ptr<lapack_complex_double> siteMatrixContainerP(new lapack_complex_double [ld*lDR*lDL]);
  lapack_complex_double *siteMatrixContainer=siteMatrixContainerP.get();
  lapack_complex_double *currentM, *LTerm, *RTerm;
  siteQNOrderMatrix const localIndexTable=networkState.indexTable().getLocalIndexTable(i);
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  networkState.subMatrixStart(currentM,i);
  //To compute the current energy, we multiply the current state with H and contract with the current state. This is done locally, the cached partial contractions are used for the rest of the network
  if(i==0 || i==L-1 || !pars.nQNs){
    optHMatrix gather(RTerm,LTerm,&networkH,networkDimInfo,networkH.maxDim(),i,0,0,0);
    gather.MultMv(currentM,siteMatrixContainer);
  }
  else{
    blockHMatrix BMat(RTerm,LTerm,&networkH,networkDimInfo,Dw,i,&(localIndexTable),0,0,&conservedQNs,1);
    BMat.prepareInput(currentM);
    BMat.MultMvBlockedLP(BMat.getCompressedVector(),BMat.getCompressedVector());
    BMat.readOutput(siteMatrixContainer);
  }
  //Here, we contract the result with the current state. There is a version with and one without the use of QNs
  if(pars.nQNs && 0){
    int const numBlocks=localIndexTable.numBlocksLP();
    int lBlockSize, rBlockSize;
    int siCurrent, aiCurrent, aimCurrent;
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=localIndexTable.lBlockSizeLP(iBlock);
      rBlockSize=localIndexTable.rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	aimCurrent=localIndexTable.aimBlockIndexLP(iBlock,k);
	siCurrent=localIndexTable.siBlockIndexLP(iBlock,k);
	for(int j=0;j<rBlockSize;++j){
	  aiCurrent=localIndexTable.aiBlockIndexLP(iBlock,j);
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
  return real(simpleContainer);
}

//---------------------------------------------------------------------------------------------------//

void network::getNewAlpha(int i, double &lambda, double prevLambda){
  //Is this worht the effort? (one additional application of H + one contraction of two tensors of rank 3
  double lambdaBuf=getCurrentEnergy(i);
  double dE0=prevLambda-lambda;

  if(dE0<0){
    //At this point, we have to make sure dE0 is positive. If it is negative, this is a problem which has to be handled elsewhere
    dE0*=-1.0;
  }
  double const pre=pow(0.1,1/static_cast<double>(L));
  double const tol=1e-12;
  double dET=lambdaBuf-(lambda-shift);
  if(dET<0){
    alpha*=pre;
  }
  else{
    if(dE0<tol){
      alpha/=pre;
    }
    else{
      //Both dET and dE0 are now positive
      if(dET/dE0<0.3){
	alpha*=pre;
      }
      else{
	alpha/=pre;
      }
    }
  }
  lambda=lambdaBuf;
  //This is numerically somewhat simpler but a more sophisticated scheme is desirable
  //alpha*=0.9765;
  //alpha*=pow(0.1,1/static_cast<double>(2*L));
  std::cout<<"New alpha="<<alpha<<std::endl;
}
