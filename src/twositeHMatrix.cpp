#include "twositeHMatrix.h"
#include "tmpContainer.h"

twositeHMatrix::twositeHMatrix(arcomplex<double> *R, arcomplex<double> *L, mpo<arcomplex<double> > *Hin, dimensionTable const &dimInfoIn, twositeQNOrderMatrix *indexTableIn, projector *excitedStateP):
  HMPO(Hin),
  dimInfo(dimInfoIn),
  Lctr(L),
  Rctr(R),
  indexTable(indexTableIn)
{
  int cBlockSize;
  int const numBlocks=indexTable->numBlocks();
  blockOffset.resize(numBlocks);
  blockOffset[0]=0;
  dimension=0;
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    cBlockSize=indexTable->lBlockSize(iBlock)*indexTable->rBlockSize(iBlock);
    if(iBlock<numBlocks-1){
      blockOffset[iBlock+1]=blockOffset[iBlock]+cBlockSize;
    }
    dimension+=cBlockSize;
  }
  
  std::cout<<"Eigenvalue problem dimension: "<<dimension<<std::endl;

  compressedVector=new arcomplex<double>[dimension];
  i=indexTable->getSite();
  ld=dimInfo.locd(i);
  ldp=dimInfo.locd(i+1);
  lDL=dimInfo.locDimL(i);
  lDRR=dimInfo.locDimR(i+1);
  lDwL=HMPO->locDimL(i);
  lDwRR=HMPO->locDimR(i+1);
  Dw=HMPO->maxDim();
  D=dimInfo.D();
  int const lDwR=HMPO->locDimR(i);
  arcomplex<double> simpleContainer;
  W=new arcomplex<double>[ld*ld*ldp*ldp*lDwRR*lDwL];
  for(int si=0;si<ld;++si){
    for(int sip=0;sip<ld;++sip){
      for(int sit=0;sit<ldp;++sit){
	for(int sitp=0;sitp<ldp;++sitp){
	  for(int bir=0;bir<lDwRR;++bir){
	    for(int bim=0;bim<lDwL;++bim){
	      simpleContainer=0;
	      for(int bi=0;bi<lDwR;++bi){
		simpleContainer+=HMPO->global_access(i,si,sip,bi,bim)*HMPO->global_access(i+1,sit,sitp,bir,bi);
	      }
	    }
	  }
	}
      }
    }
  }
}

twositeHMatrix::~twositeHMatrix(){
  delete[] compressedVector;
  delete[] W;
}

//---------------------------------------------------------------------------------------------------//

void twositeHMatrix::MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w){
  int lBlockSize, rBlockSize;
  int const numBlocks=indexTable->numBlocks();
  dynamic5DContainer<arcomplex<double> > innerContainer, outerContainer;
  innerContainer.generate(ld,ldp,lDL,lDRR,lDwL);
  arcomplex<double> *dummy;
  innerContainer.getPtr(dummy);
  int aimB,siB,sipB,airB;
  for(int m=0;m<ld*ldp*lDL*lDRR*lDwL;++m){
    dummy[m]=0;
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSize(iBlock);
    rBlockSize=indexTable->rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      airB=indexTable->airBlockIndex(iBlock,j);
      sipB=indexTable->sipBlockIndex(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimB=indexTable->aimBlockIndex(iBlock,k);
	siB=indexTable->siBlockIndex(iBlock,k);
	for(int aim=0;aim<lDL;++aim){
	  for(int bim=0;bim<lDwL;++bim){
	    innerContainer.global_access(siB,sipB,aim,airB,bim)+=v[vecIndex(sipB,airB,siB,aimB)]*Lctr[ctrIndex(aim,bim,aimB)];
	  }
	}
      }
    }
  }
  outerContainer.generate(ld,ldp,lDL,lDwRR,lDRR);
  outerContainer.getPtr(dummy);
  for(int m=0;m<ld*ldp*lDL*lDwRR*lDRR;++m){
    dummy[m]=0;
  }
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSize(iBlock);
    rBlockSize=indexTable->rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      airB=indexTable->airBlockIndex(iBlock,j);
      sipB=indexTable->sipBlockIndex(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	siB=indexTable->siBlockIndex(iBlock,k);
	for(int aim=0;aim<lDL;++aim){
	  for(int sip=0;sip<ld;++sip){
	    for(int sitp=0;sitp<ldp;++sitp){
	      for(int bir=0;bir<lDwRR;++bir){
		for(int bim=0;bim<lDwL;++bim){
		  outerContainer.global_access(sip,sitp,aim,bir,airB)+=innerContainer.global_access(siB,sipB,aim,airB,bim)*W[hIndex(siB,sip,sipB,sitp,bim,bir)];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  arcomplex<double> simpleContainer;
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSize(iBlock);
    rBlockSize=indexTable->rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      airB=indexTable->airBlockIndex(iBlock,j);
      sipB=indexTable->sipBlockIndex(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	aimB=indexTable->aimBlockIndex(iBlock,k);
	siB=indexTable->siBlockIndex(iBlock,k);
	simpleContainer=0;
	for(int bir=0;bir<lDwRR;++bir){
	  for(int air=0;air<lDRR;++air){
	    simpleContainer+=outerContainer.global_access(siB,sipB,aimB,bir,air)*Rctr[ctrIndex(airB,bir,air)];
	  }
	}
      }
    }
  }	    
}

//---------------------------------------------------------------------------------------------------//

void twositeHMatrix::storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed){
  int lBlockSize, rBlockSize;
  int const numBlocks=indexTable->numBlocks();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSize(iBlock);
    rBlockSize=indexTable->rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vCompressed[vecBlockIndex(iBlock,j,k)]=v[vecIndex(indexTable->sipBlockIndex(iBlock,j),indexTable->airBlockIndex(iBlock,j),indexTable->siBlockIndex(iBlock,k),indexTable->aimBlockIndex(iBlock,k))];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void twositeHMatrix::storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded){
  int lBlockSize, rBlockSize;
  int const numBlocks=indexTable->numBlocks();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSize(iBlock);
    rBlockSize=indexTable->rBlockSize(iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vExpanded[vecIndex(indexTable->sipBlockIndex(iBlock,j),indexTable->airBlockIndex(iBlock,j),indexTable->siBlockIndex(iBlock,k),indexTable->aimBlockIndex(iBlock,k))]=v[vecBlockIndex(iBlock,j,k)];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void twositeHMatrix::readOutput(arcomplex<double> *outputVector){
  storageExpand(compressedVector,outputVector);
}

//---------------------------------------------------------------------------------------------------//

void twositeHMatrix::prepareInput(arcomplex<double> *inputVector){
  storageCompress(inputVector,compressedVector);
}

