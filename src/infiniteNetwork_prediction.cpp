#include "infiniteNetwork.h"
#include "mkl_complex_defined.h"
#include "arrayprocessing.h"
#include <memory>
#include <vector>
#include <algorithm>

//---------------------------------------------------------------------------------------------------//

bool compareSortData(sortData const &a, sortData const &b){
  return a.lambda>b.lambda;
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::statePrediction(arcomplex<double> *target){
  //Submatrices of already existing matrices
  arcomplex<double> *subMatrix, *aMatrix, *bMatrix;
  //Buffer matrices for intermediate results
  arcomplex<double> *leftBuf, *rightBuf;
  arcomplex<double> const zzero=0.0;
  arcomplex<double> const zone=1.0;
  int const lDL=dimInfo.locDimL(i-1);
  int const lDRR=lDL;
  int const lDR=dimInfo.locDimR(i-1);
  std::unique_ptr<arcomplex<double> > leftBufP(new arcomplex<double> [lDL*lDR]);
  std::unique_ptr<arcomplex<double> > rightBufP(new arcomplex<double> [lDL*lDR]);
  leftBuf=leftBufP.get();
  rightBuf=rightBufP.get();
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int sip=0;sip<dimInfo.locd(i+1);++sip){
      subMatrix=target+sip*lDL*lDRR+si*dimInfo.locd(i+1)*lDL*lDRR;
      networkState->subMatrixStart(bMatrix,i+2,si);
      networkState->subMatrixStart(aMatrix,i-1,sip);
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  leftBuf[ai+lDR*aim]=bMatrix[ai+lDR*aim]*diags[ai];
	  rightBuf[aim+lDL*ai]=aMatrix[aim+lDL*ai]*diags[ai]/diagsm[aim];
	}
      }
      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDR,lDL,&zone,leftBuf,lDR,rightBuf,lDL,&zzero,subMatrix,lDR);
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
  std::unique_ptr<arcomplex<double> > sourceBlockP, aBlockP, bBlockP;
  std::unique_ptr<arcomplex<double> > aFullP(new arcomplex<double> [ld*lDL*ld*lDL]);
  std::unique_ptr<arcomplex<double> > bFullP(new arcomplex<double> [ldp*lDRR*ldp*lDRR]);
  std::unique_ptr<double > diagsFullP(new double [ld*lDL]);
  arcomplex<double> *sourceBlock, *aBlock, *bBlock, *aFull, *bFull;
  std::unique_ptr<double> diagsBlockP;
  double *diagsFull, *diagsBlock;
  aFull=aFullP.get();
  bFull=bFullP.get();
  diagsFull=diagsFullP.get();
  networkState->subMatrixStart(aMatrix,i);
  networkState->subMatrixStart(bMatrix,i+1);

  
  int const nQNs=networkState->centralIndexTable().nQNs();
  optLocalQNs.resize(ld*lDL);
  for(int m=0;m<ld*lDL;++m){
    optLocalQNs[m]=std::complex<int>(-100,2);
  }

  int const numBlocks=networkState->centralIndexTable().numBlocks();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    rBlockSize=networkState->centralIndexTable().rBlockSize(iBlock);
    lBlockSize=networkState->centralIndexTable().lBlockSize(iBlock);
    if(lBlockSize!=0 && rBlockSize!=0){
      sourceBlockP.reset(new arcomplex<double>[rBlockSize*lBlockSize]);
      aBlockP.reset(new arcomplex<double>[lBlockSize*lBlockSize]);
      bBlockP.reset(new arcomplex<double>[rBlockSize*rBlockSize]);
      diagsBlockP.reset(new double[lBlockSize]);
      sourceBlock=sourceBlockP.get();
      aBlock=aBlockP.get();
      bBlock=bBlockP.get();
      diagsBlock=diagsBlockP.get();
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  sourceBlock[k+lBlockSize*j]=source[explicitIndex(iBlock,j,k)];
	}
      }
      LAPACKE_zgesdd(LAPACK_COL_MAJOR,'A',lBlockSize,rBlockSize,sourceBlock,lBlockSize,diagsBlock,aBlock,lBlockSize,bBlock,rBlockSize);
      for(int k=0;k<lBlockSize;++k){
	aimB=networkState->centralIndexTable().aimBlockIndex(iBlock,k);
	siB=networkState->centralIndexTable().siBlockIndex(iBlock,k);
	for(int kp=0;kp<lBlockSize;++kp){
	  airB=networkState->centralIndexTable().aimBlockIndex(iBlock,kp);
	  sipB=networkState->centralIndexTable().siBlockIndex(iBlock,kp);
	  aFull[airB+sipB*lDL+aimB*lDL*ld+siB*lDL*lDL*ld]=aBlock[kp+lBlockSize*k];
	}
	if(k<rBlockSize){
	  diagsFull[aimB+lDL*siB]=diagsBlock[k];
	  optLocalQNs[aimB+lDL*siB]=networkState->centralIndexTable().blockQN(0,iBlock);
	}
      }
      for(int j=0;j<rBlockSize;++j){
	aimB=networkState->centralIndexTable().airBlockIndex(iBlock,j);
	siB=networkState->centralIndexTable().sipBlockIndex(iBlock,j);
	for(int jp=0;jp<rBlockSize;++jp){
	  airB=networkState->centralIndexTable().airBlockIndex(iBlock,jp);
	  sipB=networkState->centralIndexTable().sipBlockIndex(iBlock,jp);
	  bFull[airB+sipB*lDRR+aimB*lDRR*ld+siB*lDRR*lDRR*ld]=bBlock[jp+rBlockSize*j];
	}
      }	
    }
  }

  //BLOCKWISE TRUNCATION REQUIRED FOR INITIAL STATE SEARCH
  std::vector<sortData> comparer;
  comparer.resize(ld*lDL);
  for(int ai=0;ai<ld*lDL;++ai){
    comparer[ai].index=ai;
    comparer[ai].QN=optLocalQNs[ai];
    comparer[ai].lambda=diagsFull[ai];
    //std::cout<<"Available QN: "<<optLocalQNs[ai]<<std::endl;
  }

  sort(comparer.begin(),comparer.end(),compareSortData);

  //Truncate to lDR largest SVs
  optLocalQNs.resize(lDR);
  diags.resize(lDR);
  for(int ai=0;ai<lDR;++ai){
    //Copy the remaining diagonal matrix into diags
    optLocalQNs[ai]=comparer[ai].QN;
    diags[ai]=comparer[ai].lambda;
  }
   
  for(int ai=0;ai<lDR;++ai){
    //std::cout<<"Kept QN label "<<optLocalQNs[ai]<<" with SV "<<diags[ai]<<std::endl;
  }

  //Copy the corresponding columns(A)/rows(B) into the MPS
  for(int ai=0;ai<lDR;++ai){
    for(int si=0;si<ld;++si){
      for(int aim=0;aim<lDL;++aim){
	aMatrix[aim+ai*lDL+si*lDL*lDR]=aFull[aim+si*lDL+ld*lDL*comparer[ai].index];
      }
      for(int air=0;air<lDRR;++air){
	bMatrix[ai+air*lDR+si*lDR*lDRR]=bFull[comparer[ai].index+ld*lDRR*air+ld*lDRR*lDRR*si];
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

//---------------------------------------------------------------------------------------------------//
