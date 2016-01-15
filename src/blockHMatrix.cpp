#include "blockHMatrix.h"
#include "tmpContainer.h"
#include "arrayprocessing.h"
#include <iostream>
#include <time.h>


blockHMatrix::blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, int sweepDirectionIn, basisQNOrderMatrix *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein),
  conservedQNsB(conservedQNsin),
  sweepDirection(1)
{
  int cBlockSize;
  int numBlocks;
  if(sweepDirection){
    numBlocks=indexTable->numBlocksLP(i);
  }
  else{
    numBlocks=indexTable->numBlocksRP(i);
  }
  dimension=0;
  blockOffset.clear();
  blockOffset.push_back(0);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    if(sweepDirection){
      cBlockSize=(indexTable->lBlockSizeLP(i,iBlock))*(indexTable->rBlockSizeLP(i,iBlock));
    }
    else{
      cBlockSize=(indexTable->lBlockSizeRP(i,iBlock))*(indexTable->rBlockSizeRP(i,iBlock));
    }
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

void blockHMatrix::MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w){
  if(sweepDirection){
    MultMvBlockedLP(v,w);
  }
  else{
    MultMvBlockedRP(v,w);
  }
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlockedLP(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innerContainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outerContainer(d,lDwL,lDR,lDL);
  arcomplex<double> simpleContainer;
  int const numBlocks=indexTable->numBlocksLP(i);
  int const aimBlockSize=indexTable->aimBlockSizeSplit(i,0);
  int lBlockSize, rBlockSize, siBlockSize, rBlockSizep;
  clock_t curtime;
  curtime=clock();
  //excitedStateProject(v,i);
  for(int bi=0;bi<lDwR;++bi){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      for(int iBlockp=0;iBlockp<numBlocks;++iBlockp){
	lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
	rBlockSizep=indexTable->rBlockSizeLP(i,iBlockp);
	rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
	for(int k=0;k<lBlockSize;++k){
	  for(int jp=0;jp<rBlockSizep;++jp){
	    simpleContainer=0;
	    for(int j=0;j<rBlockSize;++j){
	      simpleContainer+=Rctr[ctrIndex(indexTable->aiBlockIndexLP(i,iBlockp,jp),bi,indexTable->aiBlockIndexLP(i,iBlock,j))]*v[vecBlockIndexLP(iBlock,j,k)];
	    }
	    innerContainer.global_access(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aimBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlockp,jp),bi)=simpleContainer;
	  }
	}
      }
    }
  }
  for(int bim=0;bim<lDwL;++bim){
    for(int si=0;si<d;++si){
      for(int iBlock=0;iBlock<numBlocks;++iBlock){
	for(int iBlockp=0;iBlockp<numBlocks;++iBlockp){
	  lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
	  rBlockSize=indexTable->rBlockSizeLP(i,iBlockp);
	  for(int kp=0;kp<aimBlockSize;++kp){
	    siBlockSize=indexTable->siBlockSizeSplitFixedaim(i,iBlock,kp);
	    for(int jp=0;jp<rBlockSize;++jp){
	      simpleContainer=0;
	      for(int k=0;k<siBlockSize;++k){
		for(int bi=0;bi<lDwR;++bi){
		  simpleContainer+=H[hIndex(si,indexTable->siBlockIndexSplitFixedaim(i,iBlock,kp,k),bi,bim)]*innerContainer.global_access(indexTable->siBlockIndexSplitFixedaim(i,iBlock,kp,k),indexTable->aimBlockIndexSplit(i,iBlock,kp),indexTable->aiBlockIndexLP(i,iBlockp,jp),bi);
		}
	      }
	      outerContainer.global_access(si,bim,indexTable->aiBlockIndexLP(i,iBlockp,jp),indexTable->aimBlockIndexSplit(i,iBlock,kp))=simpleContainer;
	    }
	  }
	}
	/*
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
	*/
      }
    }
  }	 
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	simpleContainer=0;
	for(int bim=0;bim<lDwL;++bim){
	  for(int kp=0;kp<aimBlockSize;++kp){
	    simpleContainer+=Lctr[ctrIndex(indexTable->aimBlockIndexLP(i,iBlock,k),bim,indexTable->aimBlockIndexSplit(i,iBlock,kp))]*outerContainer.global_access(indexTable->siBlockIndexLP(i,iBlock,k),bim,indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexSplit(i,iBlock,kp));
	  }
	}
	w[vecBlockIndexLP(iBlock,j,k)]=simpleContainer+shift*v[vecBlockIndexLP(iBlock,j,k)];
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
  if(sweepDirection){
    int const numBlocks=indexTable->numBlocksLP(i);
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  vExpanded[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))]=v[vecBlockIndexLP(iBlock,j,k)];
	}
      }
    }
  }
  else{
    int const numBlocks=indexTable->numBlocksRP(i);
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeRP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  vExpanded[vecIndex(indexTable->siBlockIndexRP(i,iBlock,j),indexTable->aiBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,k))]=v[vecBlockIndexRP(iBlock,j,k)];
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed){
  int rBlockSize, lBlockSize;
  if(sweepDirection){
    int const numBlocks=indexTable->numBlocksLP(i);
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  vCompressed[vecBlockIndexLP(iBlock,j,k)]=v[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))];
	}
      }
    }
  }
  else{
    int const numBlocks=indexTable->numBlocksRP(i);
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeRP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  vCompressed[vecBlockIndexRP(iBlock,j,k)]=v[vecIndex(indexTable->siBlockIndexRP(i,iBlock,j),indexTable->aiBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,k))];
	}
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

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlockedRP(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innerContainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outerContainer(d,lDwL,lDR,lDL);
  arcomplex<double> simpleContainer;
  int const numBlocks=indexTable->numBlocksRP(i);
  int lBlockSize, rBlockSize, siBlockSize, aiBlockSize;
  clock_t curtime;
  curtime=clock();
  //excitedStateProject(v,i);
  for(int bi=0;bi<lDwR;++bi){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeRP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(i,iBlock);
      aiBlockSize=indexTable->aiBlockSizeSplit(i,iBlock);
      for(int k=0;k<lBlockSize;++k){
	for(int jp=0;jp<aiBlockSize;++jp){
	  for(int j=0;j<rBlockSize;++j){
	    innerContainer.global_access(indexTable->siBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,k),indexTable->aiBlockIndexSplit(i,iBlock,jp),bi)=0;
	  }
	  for(int j=0;j<rBlockSize;++j){
    	     innerContainer.global_access(indexTable->siBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,k),indexTable->aiBlockIndexSplit(i,iBlock,jp),bi)+=Rctr[ctrIndex(indexTable->aiBlockIndexSplit(i,iBlock,jp),bi,indexTable->aiBlockIndexRP(i,iBlock,j))]*v[vecBlockIndexRP(iBlock,j,k)];
	  }
	}
      }
    }
  }
  for(int bim=0;bim<lDwL;++bim){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeRP(i,iBlock);
      rBlockSize=indexTable->rBlockSizeRP(i,iBlock);
      siBlockSize=indexTable->siBlockSizeSplit(i,iBlock);
      for(int j=0;j<rBlockSize;++j){
	for(int k=0;k<lBlockSize;++k){
	  simpleContainer=0;
	  for(int bi=0;bi<lDwR;++bi){
	    for(int jp=0;jp<siBlockSize;++jp){
	      simpleContainer+=H[hIndex(indexTable->siBlockIndexRP(i,iBlock,j),indexTable->siBlockIndexSplit(i,iBlock,jp),bi,bim)]*innerContainer.global_access(indexTable->siBlockIndexSplit(i,iBlock,jp),indexTable->aimBlockIndexRP(i,iBlock,k),indexTable->aiBlockIndexRP(i,iBlock,j),bi);
	    }
	  }
	  outerContainer.global_access(indexTable->siBlockIndexRP(i,iBlock,j),bim,indexTable->aiBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,k))=simpleContainer;
	}
      }
    }
  }	 
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeRP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeRP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	simpleContainer=0;
	for(int bim=0;bim<lDwL;++bim){
	  for(int kp=0;kp<lBlockSize;++kp){
	    simpleContainer+=Lctr[ctrIndex(indexTable->aimBlockIndexRP(i,iBlock,k),bim,indexTable->aimBlockIndexRP(i,iBlock,kp))]*outerContainer.global_access(indexTable->siBlockIndexRP(i,iBlock,j),bim,indexTable->aiBlockIndexRP(i,iBlock,j),indexTable->aimBlockIndexRP(i,iBlock,kp));
	  }
	}
	w[vecBlockIndexRP(iBlock,j,k)]=simpleContainer+shift*v[vecBlockIndexRP(iBlock,j,k)];
      }
    }
  }
  //excitedStateProject(w,i);
}
