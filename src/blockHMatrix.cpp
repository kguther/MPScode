#include "blockHMatrix.h"
#include "tmpContainer.h"
#include "arrayprocessing.h"
#include <iostream>
#include <time.h>


blockHMatrix::blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, int sweepDirectionIn, basisQNOrderMatrix *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein),
  conservedQNsB(conservedQNsin),
  sweepDirection(sweepDirectionIn)
{
  int cBlockSize;
  int const numBlocks=indexTable->numBlocksLP(i);
  dimension=0;
  blockOffset.clear();
  blockOffset.push_back(0);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    cBlockSize=(indexTable->lBlockSizeLP(i,iBlock))*(indexTable->rBlockSizeLP(i,iBlock));
    dimension+=cBlockSize;
    if(iBlock<numBlocks-1){
      blockOffset.push_back(blockOffset[iBlock]+cBlockSize);
    }
  }
  compressedVector=new arcomplex<double>[dimension];
}

//---------------------------------------------------------------------------------------------------//

blockHMatrix::~blockHMatrix(){
  delete[] compressedVector;
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlockedLP(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innerContainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outerContainer(d,lDwL,lDR,lDL);
  arcomplex<double> simpleContainer;
  int const numBlocks=indexTable->numBlocksLP(i);
  int lBlockSize, rBlockSize, siBlockSize, aimBlockSize;
  clock_t curtime;
  curtime=clock();
  //excitedStateProject(v,i);
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
	/*
	aimBlockSize=indexTable->aimBlockSizeSplit(i,iBlock);
   	for(int kp=0;kp<aimBlockSize;++kp){
	  siBlockSize=indexTable->siBlockSizeSplitFixedaim(i,iBlock,kp);
	  for(int j=0;j<rBlockSize;++j){
	    simpleContainer=0;
	    for(int k=0;k<siBlockSize;++k){
	      for(int bi=0;bi<lDwR;++bi){
		simpleContainer+=H[hIndex(si,indexTable->siBlockIndexSplitFixedaim(i,iBlock,kp,k),bi,bim)]*innerContainer.global_access(indexTable->siBlockIndexSplitFixedaim(i,iBlock,kp,k),indexTable->aimBlockIndexSplit(i,iBlock,kp),indexTable->aiBlockIndexLP(i,iBlock,j),bi);
	      }
	    }
	    outerContainer.global_access(si,bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexSplit(i,iBlock,kp))=simpleContainer;
	  }
	}
	*/
	for(int j=0;j<rBlockSize;++j){
	  for(int k=0;k<lBlockSize;++k){
	    outerContainer.global_access(si,bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))=0;
	  }
	  for(int k=0;k<lBlockSize;++k){
	    for(int bi=0;bi<lDwR;++bi){
	      outerContainer.global_access(si,bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))+=H[hIndex(si,indexTable->siBlockIndexLP(i,iBlock,k),bi,bim)]*innerContainer.global_access(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aimBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),bi);
	    }
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
	w[vecBlockIndex(iBlock,j,k)]=simpleContainer+shift*v[vecBlockIndex(iBlock,j,k)];
      }
    }
  }
  //excitedStateProject(w,i);
  if(0){
  curtime=clock()-curtime;
  std::cout<<"Matrix multiplication took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  exit(1);
  }
}

//---------------------------------------------------------------------------------------------------//
// Function used as an interface to the projector class used in the computation of excited states.
//---------------------------------------------------------------------------------------------------//

void blockHMatrix::excitedStateProject(arcomplex<double> *v, int const i){
  arcomplex<double> *vExpanded=new arcomplex<double>[d*lDR*lDL];
  storageExpand(v,vExpanded);
  P->project(v,i);
  storageCompress(vExpanded,v);
  delete[] vExpanded;
}

//---------------------------------------------------------------------------------------------------//
// These funcions can swap between the full storage scheme as used in the mps and the block-storage
// scheme as used in this class.
//---------------------------------------------------------------------------------------------------//

void blockHMatrix::storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded){
  int rBlockSize, lBlockSize;
  int const numBlocks=indexTable->numBlocksLP(i);
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

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed){
  int rBlockSize, lBlockSize;
  int const numBlocks=indexTable->numBlocksLP(i);
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

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::prepareInput(arcomplex<double> *inputVector){
  storageCompress(inputVector,compressedVector);
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::readOutput(arcomplex<double> *outputVector){
  storageExpand(compressedVector,outputVector);
}
