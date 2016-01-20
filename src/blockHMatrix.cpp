#include "blockHMatrix.h"
#include "tmpContainer.h"
#include "arrayprocessing.h"
#include <iostream>
#include <time.h>
#include "mkl_complex_defined.h"


blockHMatrix::blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, arcomplex<double> *Hin, dimensionTable &dimInfo, int Dwin, int iIn, basisQNOrderMatrix *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein),
  conservedQNsB(conservedQNsin)
{
  int cBlockSize;
  int numBlocks;
  numBlocks=indexTable->numBlocksLP(i);
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
  std::cout<<"Current eigenvalue problem dimension: "<<dimension<<std::endl;
  buildSparseHBlocked();
}

//---------------------------------------------------------------------------------------------------//

blockHMatrix::~blockHMatrix(){
  delete[] compressedVector;
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w){
  lapack_complex_double *proxy=new lapack_complex_double[dimension];
  arraycpy(dimension,v,proxy);
  arcomplex<double> simpleContainer;
  for(int m=0;m<dimension;++m){
    simpleContainer=0;
    for(int spIndex=rowPtr[m];spIndex<rowPtr[m+1];++spIndex){
      simpleContainer+=sparseMatrix[spIndex]*proxy[colIndices[spIndex]];
    }
    w[m]=simpleContainer+shift*proxy[m];
  }
  delete[] proxy;
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

void blockHMatrix::buildSparseHBlocked(){
  clock_t curtime;
  curtime=clock();
  double treshold=1e-12;
  arcomplex<double> currentEntry;
  int lBlockSize, rBlockSize, lBlockSizep, rBlockSizep;
  sparseMatrix.clear();
  rowPtr.clear();
  colIndices.clear();
  int const numBlocks=indexTable->numBlocksLP(i);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(i,iBlock);
    rBlockSize=indexTable->rBlockSizeLP(i,iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	rowPtr.push_back(sparseMatrix.size());
	for(int iBlockp=0;iBlockp<numBlocks;++iBlockp){
	  lBlockSizep=indexTable->lBlockSizeLP(i,iBlockp);
	  rBlockSizep=indexTable->rBlockSizeLP(i,iBlockp);
	  for(int jp=0;jp<rBlockSizep;++jp){
	    for(int kp=0;kp<lBlockSizep;++kp){
	      currentEntry=HEffEntry(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aimBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->siBlockIndexLP(i,iBlockp,kp),indexTable->aimBlockIndexLP(i,iBlockp,kp),indexTable->aiBlockIndexLP(i,iBlockp,jp));
	      if(abs(currentEntry)>treshold){
		sparseMatrix.push_back(currentEntry);
		colIndices.push_back(vecBlockIndexLP(iBlockp,jp,kp));
	      }
	    }
	  }
	}
      }
    }
  }
  rowPtr.push_back(sparseMatrix.size());
  /*
  curtime=clock()-curtime;
  std::cout<<"Matrix construction took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
  exit(1);
  */
}

//---------------------------------------------------------------------------------------------------//

arcomplex<double> blockHMatrix::HEffEntry(int const si, int const aim, int const ai, int const sip, int const aimp, int const aip){
  arcomplex<double> simpleContainer=0.0;
  for(int bi=0;bi<lDwR;++bi){
    if(bi==0){
      for(int bim=0;bim<lDwL;++bim){
	simpleContainer+=Lctr[ctrIndex(aim,bim,aimp)]*Rctr[ctrIndex(ai,bi,aip)]*H[hIndex(si,sip,bi,bim)];
      }
    }
    else{
      simpleContainer+=Lctr[ctrIndex(aim,lDwL-1,aimp)]*Rctr[ctrIndex(ai,bi,aip)]*H[hIndex(si,sip,bi,lDwL-1)];
    }
  }
  return simpleContainer;
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
	vExpanded[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))]=v[vecBlockIndexLP(iBlock,j,k)];
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
	vCompressed[vecBlockIndexLP(iBlock,j,k)]=v[vecIndex(indexTable->siBlockIndexLP(i,iBlock,k),indexTable->aiBlockIndexLP(i,iBlock,j),indexTable->aimBlockIndexLP(i,iBlock,k))];
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
