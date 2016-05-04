#include "infiniteNetwork.h"
#include "mkl_complex_defined.h"
#include <memory>
#include <vector>
#include <algorithm>

#include "verifyQN.h"
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// Order sortDatas with respect to lambdas, as it is used to truncate to the largest SVs
//---------------------------------------------------------------------------------------------------//


bool compareSortData(sortData const &a, sortData const &b){
  double const tol=1e-10;
  double buf=(a.lambda>b.lambda)?(a.lambda-b.lambda):(b.lambda-a.lambda);
  if(buf>tol)
    return a.lambda>b.lambda;
  if(a.QN.real()!=b.QN.real())
    return a.QN.real()>b.QN.real();
  return a.QN.imag()>b.QN.imag();
}

bool compareSortDataQNBased(sortData const &a, sortData const &b){
  if(a.QN.real()!=b.QN.real()){
    return a.QN.real()>b.QN.real();
  }
  return a.QN.imag()>b.QN.imag();  
  //The ordering has to be strictly deterministic, such that two arrays of sortData with the same entries of lambdas and QNs are ordered in the same way, even if some lambdas are degenerate
}

void verifyCompression(arcomplex<double> *cVector, int dim){
  double const norm=cblas_dznrm2(dim,cVector,1);
  double const tol=1e-7;
  if(norm<tol){
    std::cout<<"WARNING: STARTING VECTOR IS ZERO\n";
    //Maybe use sqrt(dim), but this will most likely not do a lot 
    arcomplex<double> defaultEntry=1.0/static_cast<double>(dim);
    for(int m=0;m<dim;++m){
      cVector[m]=defaultEntry;
    }
  }
  std::cout<<"Checked vector norm\n";
}

//---------------------------------------------------------------------------------------------------//

int infiniteNetwork::statePrediction(arcomplex<double> *target){
  //Submatrices of already existing matrices
  arcomplex<double> *subMatrix, *aMatrix, *bMatrix;
  //Buffer matrices for intermediate results
  arcomplex<double> *leftBuf, *rightBuf;
  arcomplex<double> const zzero=0.0;
  arcomplex<double> const zone=1.0;
  int const lDL=dimInfo.locDimL(i-1);
  int const lDRR=lDL;
  int const lDR=dimInfo.locDimR(i-1);
  std::unique_ptr<arcomplex<double>[] > leftBufP(new arcomplex<double> [lDL*lDR]);
  std::unique_ptr<arcomplex<double>[] > rightBufP(new arcomplex<double> [lDL*lDR]);
  leftBuf=leftBufP.get();
  rightBuf=rightBufP.get();
  
  int info=0;
  /*
  int info=checkQNConstraint(*networkState,i-1);
  if(info){
    std::cout<<"QN constraint violation at A-matrix"<<std::endl;
    return -1;
  }
  info=checkQNConstraint(*networkState,i+2);
  if(info){
    std::cout<<"QN constraint violation at B-matrix"<<std::endl;
    return 1;
  }
  */
  int const tol=1e-12;
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int sip=0;sip<dimInfo.locd(i+1);++sip){
      subMatrix=target+sip*lDL*lDRR+si*dimInfo.locd(i+1)*lDL*lDRR;
      bBuf.getPtr(bMatrix,si);
      aBuf.getPtr(aMatrix,sip);
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  leftBuf[ai+lDR*aim]=bMatrix[ai+lDR*aim]*diags[ai];
	  if(diagsm[aim]>tol){
	    //SVs are always positive
	    //ADD EXCEPTION
	    rightBuf[aim+lDL*ai]=aMatrix[aim+lDL*ai]*diags[ai]/diagsm[aim];
	  }
	  else{
	    rightBuf[aim+lDL*ai]=0;
	  }
	}
      }
      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDR,lDL,&zone,leftBuf,lDR,rightBuf,lDL,&zzero,subMatrix,lDR);
    }
  }

  /*
  for(int m=0;m<dimInfo.locd(i)*dimInfo.locd(i+1)*dimInfo.locDimL(i)*dimInfo.locDimL(i);++m){
    target[m]=0;
  }

  qnEnforcedPrediction(target);
  */

  std::cout<<"Guess norm: "<<cblas_dznrm2(lDL*lDL*dimInfo.locd(i)*dimInfo.locd(i+1),target,1)<<std::endl;
  
  //info=twositeCheck(*networkState,target);
  if(info){
    std::cout<<"Twosite constraint violation at site "<<networkState->currentSite()<<std::endl;
    return 2;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::qnEnforcedPrediction(arcomplex<double> *target){
  int const numBlocks=networkState->centralIndexTable().numBlocks();
  int const ld=dimInfo.locd(i);
  int const lDL=dimInfo.locDimL(i);
  int lBlockSize, rBlockSize;
  int aimB, airB, siB, sipB;
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=networkState->centralIndexTable().lBlockSize(iBlock);
    rBlockSize=networkState->centralIndexTable().rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      airB=networkState->centralIndexTable().airBlockIndex(iBlock,j);
      sipB=networkState->centralIndexTable().sipBlockIndex(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimB=networkState->centralIndexTable().aimBlockIndex(iBlock,k);
	siB=networkState->centralIndexTable().siBlockIndex(iBlock,k);
	target[aimB+lDL*airB+siB*lDL*lDL+sipB*ld*lDL*lDL]=1;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::updateMPS(arcomplex<double> *source){
  diagsm=diags;
  int const ld=dimInfo.locd(i);
  int const ldp=dimInfo.locd(i+1);
  int const lDL=dimInfo.locDimL(i);
  int const lDR=dimInfo.locDimR(i);
  int const lDRR=dimInfo.locDimR(i+1);
  int lBlockSize, rBlockSize;
  int maxTargetBlockSize, maxBlockSize;
  int aimB, siB, airB, sipB;
  arcomplex<double> *aMatrix, *bMatrix;
  std::unique_ptr<arcomplex<double>[] > sourceBlockP, aBlockP, bBlockP;
  std::unique_ptr<arcomplex<double>[] > aFullP(new arcomplex<double> [ld*lDL*ld*lDL]);
  std::unique_ptr<arcomplex<double>[] > bFullP(new arcomplex<double> [ldp*lDRR*ldp*lDRR]);
  std::unique_ptr<double[]> diagsFullLP(new double [ld*lDL]);
  std::unique_ptr<double[]> diagsFullRP(new double [ldp*lDRR]);
  arcomplex<double> *sourceBlock, *aBlock, *bBlock, *aFull, *bFull;
  std::unique_ptr<double[]> diagsBlockP;
  double *diagsFullL, *diagsFullR, *diagsBlock;
  aFull=aFullP.get();
  bFull=bFullP.get();
  diagsFullL=diagsFullLP.get();
  diagsFullR=diagsFullRP.get();
  networkState->subMatrixStart(aMatrix,i);
  networkState->subMatrixStart(bMatrix,i+1);

  //Initialize arrays with 0 where necessary
  for(int m=0;m<ld*lDL*ld*lDL;++m){
    aFull[m]=0.0;
  }
  for(int m=0;m<ldp*lDRR*ldp*lDRR;++m){
    bFull[m]=0.0;
  }
  //Note that lDL==lDRR is both required and guaranteed by the algorithm. For better readability, we distinguish
  for(int m=0;m<ld*lDL;++m){
    diagsFullL[m]=0.0;
  }
  for(int m=0;m<ldp*lDRR;++m){
    diagsFullR[m]=0.0;
  }
  
  //Invalidate block QNs
  int const nQNs=networkState->centralIndexTable().nQNs();
  optLocalQNsL=std::vector<std::complex<int> >(ld*lDL,std::complex<int>(-100,2));
  optLocalQNsR=optLocalQNsL;

  int const numBlocks=networkState->centralIndexTable().numBlocks();
  int minBlockSize;
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    rBlockSize=networkState->centralIndexTable().rBlockSize(iBlock);
    lBlockSize=networkState->centralIndexTable().lBlockSize(iBlock);
    minBlockSize=(lBlockSize>rBlockSize)?rBlockSize:lBlockSize;
    if(lBlockSize!=0 && rBlockSize!=0){
      //Setup temporary container
      sourceBlockP.reset(new arcomplex<double>[rBlockSize*lBlockSize]);
      aBlockP.reset(new arcomplex<double>[lBlockSize*lBlockSize]);
      bBlockP.reset(new arcomplex<double>[rBlockSize*rBlockSize]);
      diagsBlockP.reset(new double[minBlockSize]);
      sourceBlock=sourceBlockP.get();
      aBlock=aBlockP.get();
      bBlock=bBlockP.get();
      diagsBlock=diagsBlockP.get();
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  //Read the current block
	  sourceBlock[k+lBlockSize*j]=source[explicitIndex(iBlock,j,k)];
	}
      }
      //SVD of the current block
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',lBlockSize,rBlockSize,sourceBlock,lBlockSize,diagsBlock,aBlock,lBlockSize,bBlock,rBlockSize);

      for(int k=0;k<lBlockSize;++k){
	aimB=networkState->centralIndexTable().aimBlockIndex(iBlock,k);
	siB=networkState->centralIndexTable().siBlockIndex(iBlock,k);
	for(int kp=0;kp<lBlockSize;++kp){
	  airB=networkState->centralIndexTable().aimBlockIndex(iBlock,kp);
	  sipB=networkState->centralIndexTable().siBlockIndex(iBlock,kp);
	  //Extract the current results into the total, untruncated A
	  aFull[airB+sipB*lDL+aimB*lDL*ld+siB*lDL*lDL*ld]=aBlock[kp+lBlockSize*k];
	}

	//Get the current block QN and SV for truncation
	if(k<rBlockSize){
	  diagsFullL[aimB+lDL*siB]=diagsBlock[k];
	}
	else{
	  diagsFullL[aimB+lDL*siB]=0.0;
	}

	//This label is only bound to the left indices
	optLocalQNsL[aimB+lDL*siB]=networkState->centralIndexTable().blockQN(0,iBlock);
      }
      for(int j=0;j<rBlockSize;++j){
	aimB=networkState->centralIndexTable().airBlockIndex(iBlock,j);
	siB=networkState->centralIndexTable().sipBlockIndex(iBlock,j);
	for(int jp=0;jp<rBlockSize;++jp){
	  airB=networkState->centralIndexTable().airBlockIndex(iBlock,jp);
	  sipB=networkState->centralIndexTable().sipBlockIndex(iBlock,jp);
	  //Extract the current results into the total, untruncated B
	  bFull[airB+sipB*lDRR+aimB*lDRR*ld+siB*lDRR*lDRR*ld]=bBlock[jp+rBlockSize*j];
	}

	if(j<lBlockSize){
	  diagsFullR[aimB+lDRR*siB]=diagsBlock[j];
	}
	else{
	  diagsFullR[aimB+lDRR*siB]=0.0;
	}
	
	//This label is now only bound to the right indices
	optLocalQNsR[aimB+lDRR*siB]=networkState->centralIndexTable().blockQN(0,iBlock);

      }	
    }
  }

  //BLOCKWISE TRUNCATION REQUIRED FOR INITIAL STATE SEARCH
  std::vector<sortData> comparerL, comparerR;
  comparerL.resize(ld*lDL);
  comparerR.resize(ld*lDL);
  for(int ai=0;ai<ld*lDL;++ai){
    comparerL[ai].index=ai;
    comparerL[ai].QN=optLocalQNsL[ai];
    comparerL[ai].lambda=diagsFullL[ai];

    comparerR[ai].index=ai;
    comparerR[ai].QN=optLocalQNsR[ai];
    comparerR[ai].lambda=diagsFullR[ai];
    //std::cout<<"Available QN: "<<optLocalQNs[ai]<<std::endl;
  }

  std::sort(comparerL.begin(),comparerL.end(),compareSortData);
  std::sort(comparerR.begin(),comparerR.end(),compareSortData);

  /*
  for(int m=0;m<lDR;++m){
    std::cout<<"Left-bound SV: "<<comparerL[m].lambda<<" Right-bound SV: "<<comparerR[m].lambda<<" Left-bound QN: "<<comparerL[m].QN<<" Right-bound QN: "<<comparerR[m].QN<<std::endl;
  }
  */

  //We will need to declare D labels, even in the rare case less than D labels are used
  int const D=dimInfo.D();
  optLocalQNsL.resize(D);
  optLocalQNsR.resize(D);
  
  //In this rare case, the additional labels are never called, and do not need to be set
  
  //Truncate to lDR largest SVs
  diags.resize(lDR);

  for(int ai=0;ai<lDR;++ai){
    //Copy the remaining diagonal matrix into diags
    optLocalQNsL[ai]=comparerL[ai].QN;
    optLocalQNsR[ai]=comparerR[ai].QN;
    //Both comparers contain the same lambda-values, only bound to different QN labels
    diags[ai]=comparerL[ai].lambda;
  }

  //Copy the corresponding columns(A)/rows(B) into the MPS
  for(int ai=0;ai<lDR;++ai){
    for(int si=0;si<ld;++si){
      for(int aim=0;aim<lDL;++aim){
	aMatrix[aim+ai*lDL+si*lDL*lDR]=aFull[aim+si*lDL+ld*lDL*comparerL[ai].index];
      }
      for(int air=0;air<lDRR;++air){
	bMatrix[ai+air*lDR+si*lDR*lDRR]=bFull[comparerR[ai].index+ld*lDRR*air+ld*lDRR*lDRR*si];

	//std::cout<<"Matrix element: "<<bMatrix[ai+air*lDR+si*lDR*lDRR]<<" with labels: "<<comparerL[ai].QN<<"+"<<gqn.QNLabel(si)<<"="<<gqn.QNLabel(i+1,air)<<std::endl;
      }
    }
  }

  //Problem: TRUNCATE THE RESULTING MATRICES WITHOUT DESTROYING THE QN LABELING SCHEME
  //Truncation will be carried out by deleting individual columns/rows from a/b. This leads to an outdated QN labeling. Possible solution: Create new QN labeling here and make local QN labelings storable. Then: CHANGE HOW QNs ARE UPDATED WHEN INSERTING A SITE. This will currently destroy any stored data.
}

//---------------------------------------------------------------------------------------------------//


int infiniteNetwork::explicitIndex(int iBlock, int j, int k){
  return networkState->centralIndexTable().aimBlockIndex(iBlock,k)+dimInfo.locDimL(i)*dimInfo.locDimR(i+1)*networkState->centralIndexTable().siBlockIndex(iBlock,k)+networkState->centralIndexTable().airBlockIndex(iBlock,j)*dimInfo.locDimL(i)+networkState->centralIndexTable().sipBlockIndex(iBlock,j)*dimInfo.locDimL(i)*dimInfo.locDimR(i+1)*dimInfo.locd(i);
}
