#include "blockHMatrix.h"
#include "tmpContainer.h"

blockHMatrix::blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, basisQNOrderMatrix *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein),
  conservedQNsB(conservedQNsin)
{
  int cBlockSize;
  dimension=0;
  blockOffset[0]=0;
  for(int iBlock=0;iBlock<(indexTable->numBlocksLP(i));++iBlock){
    cBlockSize=(indexTable->lBlockSizeLP(i,iBlock))*(indexTable->rBlockSizeLP(i,iBlock));
    dimension+=cBlockSize;
    if(iBlock>0){
      blockOffset[iBlock]=blockOffset[iBlock-1]+cBlockSize;
    }
  }
}

void blockHMatrix::MultMv(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innerContainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outerContainer(d,lDwL,lDR,lDL);
  arcomplex<double> simpleContainer;
  int numBlocks=indexTable->numBlocksLP(i);
  int lBlockSize, rBlockSize, siBlockSize, aimBlockSize;
  excitedStateProject(v,i);
  for(int bi=0;bi<lDwR;++bi){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
      for(int k=0;k<lBlockSize;++k){
	for(int jp=0;jp<rBlockSize;++jp){
	  simpleContainer=0;
	  for(int j=0;j<rBlockSize;++j){
	    simpleContainer+=Rctr[ctrIndex(indexTable->aiBlockIndexLP(i,iBlock,jp),bi,indexTable->aiBlockIndexLP(i,iBlock,j))]*v[vecBlockIndex(iBlock,j,k)];
	  }
	  innerContainer.global_access(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aimBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,jp),bi)=simpleContainer;
	}
      }
    }
  }
  for(int bim=0;bim<lDwL;++bim){
    for(int si=0;si<d;++si){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
	rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
	siBlockSize=indexTable->siBlockSizeSplit(i,iBlock);
	aimBlockSize=indexTable->aimBlockSizeSplit(i,iBlock);
	for(int j=0;j<rBlockSize;++j){
	  for(int kp=0;kp<aimBlockSize;++kp){
	    simpleContainer=0;
	    for(int k=0;k<siBlockSize;++k){
	      for(int bi=0;bi<lDwR;++bi){
		simpleContainer+=H[hIndex(si,indexTable->siBlockIndexLP(i,iBlock,k),bi,bim)]*innerContainer.global_access(indexTable->siBlockIndexSplit(i,iBlock,k),indexTable->aimBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),bi);
	      }
	    }
	    outerContainer.global_access(si,bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexSplit(i,iBlock,kp))=simpleContainer;
	  }
	}
      }
    }
  }	 
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
    aimBlockSize=indexTable->aimBlockSizeSplit(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	simpleContainer=0;
	for(int bim=0;bim<lDwL;++bim){
	  for(int kp=0;kp<aimBlockSize;++kp){
	    simpleContainer+=Lctr[ctrIndex(indexTable->aimBlockIndexLP(i,iBlock,k),bim,indexTable->aimBlockIndexSplit(i,iBlock,kp))]*outerContainer.global_access(indexTable->siBlockIndexLP(i,iBlock,k),bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexSplit(i,iBlock,kp));
	  }
	}
	w[vecBlockIndex(iBlock,j,k)]=simpleContainer;
      }
    }
  }
  excitedStateProject(w,i);
}

void blockHMatrix::excitedStateProject(arcomplex<double> *v, int const i){
  arcomplex<double> *vExpanded=new arcomplex<double>[d*lDR*lDL];
  storageExpand(v,vExpanded);
  P->project(v,i);
  storageCompress(vExpanded,v);
  delete[] vExpanded;
}

void blockHMatrix::storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded){
  int rBlockSize, lBlockSize;
  int numBlocks=indexTable->numBlocksLP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vExpanded[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))]=v[vecBlockIndex(iBlock,j,k)];
      }
    }
  }
}

void blockHMatrix::storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed){
  int rBlockSize, lBlockSize;
  int numBlocks=indexTable->numBlocksLP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vCompressed[vecBlockIndex(iBlock,j,k)]=v[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))];
      }
    }
  }
}
