#include "blockHMatrix.h"
#include "tmpContainer.h"
#include "arrayprocessing.h"
#include <iostream>
#include <memory>

blockHMatrix::blockHMatrix(arcomplex<double> *R, arcomplex<double> *L, mpo<arcomplex<double> > const *Hin, dimensionTable &dimInfo, int Dwin, int iIn, siteQNOrderMatrix const *indexTablein, projector *excitedStateP, double shift, std::vector<quantumNumber> *conservedQNsin, int const cached):
  optHMatrix(R,L,Hin,dimInfo,Dwin,iIn,excitedStateP,shift,conservedQNsin),
  indexTable(indexTablein)
{
  int cBlockSize;
  int numBlocks;
  numBlocks=indexTable->numBlocksLP();
  dimension=0;
  blockOffset.clear();
  blockOffset.push_back(0);
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    cBlockSize=(indexTable->lBlockSizeLP(iBlock))*(indexTable->rBlockSizeLP(iBlock));
    dimension+=cBlockSize;
    if(iBlock<numBlocks-1){
      blockOffset.push_back(blockOffset[iBlock]+cBlockSize);
    }
  }
  compressedVector.resize(dimension);
  std::cout<<"Current eigenvalue problem dimension: "<<dimension<<std::endl;
  if(lDR<350 && lDL<350 && !cached){
    explicitMv=1;
  }
  else{
    explicitMv=0;
  }
  if(explicitMv){
    buildSparseHBlocked();
  }
}

//---------------------------------------------------------------------------------------------------//
// In contrast to the case without QN conservation, when using blocked MPS matrices, it can be much 
// faster to use a sparse matrix-vector multiplication instead of the cached approach. We use our own 
// matrix-vector multiplication to implement the projection on the excited state, which is carried
// out in the unblocked matrix (currently).
//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlocked(arcomplex<double> *v, arcomplex<double> *w){
  if(explicitMv){
    std::unique_ptr<arcomplex<double>[]> proxyP(new arcomplex<double>[dimension]);
    arcomplex<double> *proxy=proxyP.get();
    excitedStateProject(v);
    auxiliary::arraycpy(dimension,v,proxy);
    arcomplex<double> simpleContainer;
    for(int m=0;m<dimension;++m){
      simpleContainer=0;
      for(int spIndex=rowPtr[m];spIndex<rowPtr[m+1];++spIndex){
	simpleContainer+=sparseMatrix[spIndex]*proxy[colIndices[spIndex]];
      }
      w[m]=simpleContainer+shift*proxy[m];
    }
    excitedStateProject(w);
  }
  else{
    MultMvBlockedLP(v,w);
  }
}

//---------------------------------------------------------------------------------------------------//
// The previous function is very efficient for low bond dimension, but scales really poorly. Therefore,
// for higher bond dimensions (as they are required for excited state search), we return to the cached
// version of matrix vector multiplication which scales with D^3.
//---------------------------------------------------------------------------------------------------//

void blockHMatrix::MultMvBlockedLP(arcomplex<double> *v, arcomplex<double> *w){
  tmpContainer<arcomplex<double> > innerContainer(d,lDL,lDR,lDwR);
  tmpContainer<arcomplex<double> > outerContainer(d,lDwL,lDR,lDL);
  arcomplex<double> simpleContainer;
  int const numBlocks=indexTable->numBlocksLP();
  int const sparseSize=HMPO->numEls(i);
  int const nThreads=20;
  int lBlockSize, rBlockSize, siBlockSize, rBlockSizep;
  int const *biIndices, *siIndices, *bimIndices, *sipIndices;
  HMPO->biSubIndexArrayStart(biIndices,i);
  HMPO->siSubIndexArrayStart(siIndices,i);
  HMPO->bimSubIndexArrayStart(bimIndices,i);
  HMPO->sipSubIndexArrayStart(sipIndices,i);
  arcomplex<double> const *H;
  HMPO->sparseSubMatrixStart(H,i);
  int siB, aiB, aimB, sipS;
  excitedStateProject(v);
#pragma omp parallel for private(simpleContainer,siB,aimB,lBlockSize,rBlockSize) schedule(dynamic,1)
  for(int aip=0;aip<lDR;++aip){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(iBlock);
      rBlockSize=indexTable->rBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	siB=indexTable->siBlockIndexLP(iBlock,k);
	aimB=indexTable->aimBlockIndexLP(iBlock,k);
      	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0.0;
	  for(int j=0;j<rBlockSize;++j){
	    simpleContainer+=Rctr[ctrIndex(aip,bi,indexTable->aiBlockIndexLP(iBlock,j))]*v[vecBlockIndexLP(iBlock,j,k)];
	  }
	  innerContainer.global_access(siB,aimB,aip,bi)=simpleContainer;
	}
      }
    }
  }

#pragma omp parallel for schedule(static,1)
  for(int si=0;si<d;++si){
    for(int bim=0;bim<lDwL;++bim){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  outerContainer.global_access(si,bim,ai,aim)=0.0;
	}
      }
    }
  }

#pragma omp parallel for private(siB,aimB,lBlockSize,rBlockSize,sipS) schedule(dynamic,1)
  for(int ai=0;ai<lDR;++ai){
    for(int iBlock=0;iBlock<numBlocks;++iBlock){
      lBlockSize=indexTable->lBlockSizeLP(iBlock);
      for(int k=0;k<lBlockSize;++k){
	siB=indexTable->siBlockIndexLP(iBlock,k);
	aimB=indexTable->aimBlockIndexLP(iBlock,k);
	for(int nSparse=0;nSparse<sparseSize;++nSparse){
	  sipS=sipIndices[nSparse];
	  if(sipS==siB){
	    outerContainer.global_access(siIndices[nSparse],bimIndices[nSparse],ai,aimB)+=H[nSparse]*innerContainer.global_access(sipS,aimB,ai,biIndices[nSparse]);
	  }
	}
      }
    }	 
  }

#pragma omp parallel for private(simpleContainer,aiB,siB,aimB,lBlockSize,rBlockSize) schedule(dynamic,1)
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(iBlock);
    rBlockSize=indexTable->rBlockSizeLP(iBlock);
    for(int k=0;k<lBlockSize;++k){
      siB=indexTable->siBlockIndexLP(iBlock,k);
      aimB=indexTable->aimBlockIndexLP(iBlock,k);
      for(int j=0;j<rBlockSize;++j){
	aiB=indexTable->aiBlockIndexLP(iBlock,j);
	simpleContainer=0.0;
	for(int aim=0;aim<lDL;++aim){
	  for(int bim=0;bim<lDwL;++bim){
	    simpleContainer+=Lctr[ctrIndex(aimB,bim,aim)]*outerContainer.global_access(siB,bim,aiB,aim);
	  }
	}
	w[vecBlockIndexLP(iBlock,j,k)]=simpleContainer+shift*v[vecBlockIndexLP(iBlock,j,k)];
      }
    }
  }
  excitedStateProject(w);
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::buildSparseHBlocked(){
  int siB, aiB, aimB, aimpB, sipB;
  double treshold=1e-12;
  arcomplex<double> currentEntry;
  int lBlockSize, rBlockSize, lBlockSizep, rBlockSizep;
  sparseMatrix.clear();
  rowPtr.clear();
  colIndices.clear();
  int const numBlocks=indexTable->numBlocksLP();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(iBlock);
    rBlockSize=indexTable->rBlockSizeLP(iBlock);
    for(int j=0;j<rBlockSize;++j){
      aiB=indexTable->aiBlockIndexLP(iBlock,j);
      for(int k=0;k<lBlockSize;++k){
	siB=indexTable->siBlockIndexLP(iBlock,k);
	aimB=indexTable->aimBlockIndexLP(iBlock,k);
	rowPtr.push_back(sparseMatrix.size());
	for(int iBlockp=0;iBlockp<numBlocks;++iBlockp){
	  lBlockSizep=indexTable->lBlockSizeLP(iBlockp);
	  rBlockSizep=indexTable->rBlockSizeLP(iBlockp);
	  for(int kp=0;kp<lBlockSizep;++kp){
	    sipB=indexTable->siBlockIndexLP(iBlockp,kp);
	    aimpB=indexTable->aimBlockIndexLP(iBlockp,kp);
	    for(int jp=0;jp<rBlockSizep;++jp){
	      currentEntry=HEffEntry(siB,aimB,aiB,sipB,aimpB,indexTable->aiBlockIndexLP(iBlockp,jp));
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
}

//---------------------------------------------------------------------------------------------------//
// This function computes the entry of the local effective Hamiltonian. While this is only used implicitly
// in the normal algorithm, in the presence of conserved QNs, it is more efficient to explicitly 
// calculate the matrix elements.
//---------------------------------------------------------------------------------------------------//

arcomplex<double> blockHMatrix::HEffEntry(int const si, int const aim, int const ai, int const sip, int const aimp, int const aip){
  arcomplex<double> simpleContainer=0.0;
  arcomplex<double> const *H;
  HMPO->subMatrixStart(H,i);
  for(int bi=0;bi<lDwR;++bi){
    //Use the canonical form of a local MPO. This only works for nearest-neighbour terms, for longer ranged interactions, a general sparse storage has to be used.
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

void blockHMatrix::excitedStateProject(arcomplex<double> *v){
  if(P){
    if(P->nEigen()){
      arcomplex<double> *vExpanded=new arcomplex<double>[d*lDR*lDL];
      for(int m=0;m<d*lDR*lDL;++m){
	vExpanded[m]=0;
      }
      storageExpand(v,vExpanded);
      P->project(vExpanded,i);
      storageCompress(vExpanded,v);
      delete[] vExpanded;
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// These funcions can swap between the full storage scheme as used in the mps and the block-storage
// scheme as used in this class.
//---------------------------------------------------------------------------------------------------//

void blockHMatrix::storageExpand(arcomplex<double> *v, arcomplex<double> *vExpanded){
  int rBlockSize, lBlockSize;
  int const numBlocks=indexTable->numBlocksLP();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(iBlock);
    rBlockSize=indexTable->rBlockSizeLP(iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vExpanded[vecIndex(indexTable->siBlockIndexLP(iBlock,k),indexTable->aiBlockIndexLP(iBlock,j),indexTable->aimBlockIndexLP(iBlock,k))]=v[vecBlockIndexLP(iBlock,j,k)];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::storageCompress(arcomplex<double> *v, arcomplex<double> *vCompressed){
  int rBlockSize, lBlockSize;
  int const numBlocks=indexTable->numBlocksLP();
  for(int iBlock=0;iBlock<numBlocks;++iBlock){
    lBlockSize=indexTable->lBlockSizeLP(iBlock);
    rBlockSize=indexTable->rBlockSizeLP(iBlock);
    for(int j=0;j<rBlockSize;++j){
      for(int k=0;k<lBlockSize;++k){
	vCompressed[vecBlockIndexLP(iBlock,j,k)]=v[vecIndex(indexTable->siBlockIndexLP(iBlock,k),indexTable->aiBlockIndexLP(iBlock,j),indexTable->aimBlockIndexLP(iBlock,k))];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::prepareInput(arcomplex<double> *inputVector){
  storageCompress(inputVector,getCompressedVector());
}

//---------------------------------------------------------------------------------------------------//

void blockHMatrix::readOutput(arcomplex<double> *outputVector){
  storageExpand(getCompressedVector(),outputVector);
}
