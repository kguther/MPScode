#include <float.h>
#include "projector.h"
#include "siteArray.h"
#include "mps.h"
#include "overlap.h"
#include "arrayprocessing.h"

//---------------------------------------------------------------------------------------------------//
// The projector class is still experimental, it seems to yield correct results, but I am not sure
// whether everything is correct.
//---------------------------------------------------------------------------------------------------//

projector::projector(){
  orthoStates=0;
  scalarProducts=0;
  nCurrentEigen=0;
}

//---------------------------------------------------------------------------------------------------//

projector::~projector(){
  delete[] orthoStates;
  delete[] scalarProducts;
}

//---------------------------------------------------------------------------------------------------//

void projector::initialize(int const nEigsin){
  nEigs=nEigsin;
  orthoStates=new mps[nEigs];
  scalarProducts=new overlap[nEigs-1];
}

//---------------------------------------------------------------------------------------------------//
// When adjusting the simulation parameter D, all stored states are updated
//---------------------------------------------------------------------------------------------------//

void projector::setParameterD(int const Dnew){
  for(int iEigen=0;iEigen<nEigs;++iEigen){
    orthoStates[iEigen].setParameterD(Dnew);
  }
}

//---------------------------------------------------------------------------------------------------//
// This functions loads the state which shall be orthogonalized, and uses its index to determine to
// which states it shall be orthogonalized.
//---------------------------------------------------------------------------------------------------//

void projector::loadScalarProducts(mps *const variationalState,int const iEigen){
  //offset marks the beginning of the scalar products with the current state
  for(int k=0;k<iEigen;++k){
    scalarProducts[k].loadMPS(variationalState,&orthoStates[k]);
  }
  nCurrentEigen=iEigen;
}

//---------------------------------------------------------------------------------------------------//

void projector::updateScalarProducts(int const i, int const direction){
  if(nCurrentEigen>0){
    for(int k=0;k<nCurrentEigen;++k){
      if(direction==1){
	scalarProducts[k].stepRight(i);
      }
      if(direction==-1){
	scalarProducts[k].stepLeft(i);
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// Besides featuring the obvious functionality of projecting out already obtained states, the projector
// class also serves as a storage for the MPS representations of the states obtained when searching
// excited states (since the network only has a single state). These functions grant access to the
// stored states in various ways.
//---------------------------------------------------------------------------------------------------//


void projector::storeOrthoState(mps const &source, int const iEigen){
  orthoStates[iEigen].mpsCpy(source);
}

//---------------------------------------------------------------------------------------------------//

void projector::storeCurrentState(mps &source){
  orthoStates[nCurrentEigen].mpsCpy(source);
}

//---------------------------------------------------------------------------------------------------//

void projector::getStoredState(mps *target, int const iEigen){
  target=&(orthoStates[iEigen]);
}

//---------------------------------------------------------------------------------------------------//

int projector::loadNextState(mps &target, int const iEigen){
  target.mpsCpy(orthoStates[iEigen]);
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int projector::loadNextState(mps &target){
  if(nCurrentEigen<(nEigs-1)){
    target.mpsCpy(orthoStates[nCurrentEigen+1]);
    return 0;
  }
  return 1;
}		

//---------------------------------------------------------------------------------------------------//

void projector::getLocalDimensions(int const i){
  //Note that it is required for all stored states to have the same bond dimension
  if(nEigs>0){
    ld=orthoStates[0].locd(i);
    lDL=orthoStates[0].locDimL(i);
    lDR=orthoStates[0].locDimR(i);
  }
}

//---------------------------------------------------------------------------------------------------//

lapack_complex_double projector::fullOverlap(int const k){
  return scalarProducts[k].fullOverlap();
}

//---------------------------------------------------------------------------------------------------//
// Although it is self-explanatory what project does, keep in mind that one has to use getProjector
// for the current site before using project to ensure the auxiliary matrix is set.
//---------------------------------------------------------------------------------------------------//

void projector::project(lapack_complex_double *vec, int const i){
  //vec is required to be of dimension lDL x ld*lDR
  if(nCurrentEigen>0){
    lapack_complex_double *trContainer;
    lapack_complex_double *G;
    lapack_complex_double *vecContainer;
    lapack_complex_double simpleContainer;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    getLocalDimensions(i);
    vecContainer=new lapack_complex_double[ld*lDR*lDL];
    trContainer=new lapack_complex_double[ld*lDR*ld*lDR];
    //Initialization is required (i.e. it is as fast as any other way)
    for(int mi=0;mi<ld*lDR*lDL;++mi){
      vecContainer[mi]=0;
    }
    for(int mu=0;mu<nRelevantEigens;++mu){
      auxiliaryMatrix.subMatrixStart(G,mu);
      //Compute Tr(G^dagger *vec) which gives the projector onto the space spanned by lower lying states
      cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,ld*lDR,ld*lDR,lDL,&zone,G,lDL,vec,lDL,&zzero,trContainer,ld*lDR);
      simpleContainer=0;
      for(int mi=0;mi<ld*lDR;++mi){
	simpleContainer+=trContainer[mi+ld*lDR*mi];
      }
      cblas_zaxpy(ld*lDR*lDL,&simpleContainer,G,1,vecContainer,1);
    }
    //Use the projector to the orthogonal complement of the space spanned by lower lying states
    for(int mi=0;mi<ld*lDR*lDL;++mi){
      vec[mi]-=vecContainer[mi];
    }
    delete[] trContainer;
    delete[] vecContainer;
  }
}

//---------------------------------------------------------------------------------------------------//
// Compute the auxiliary matrices of the projector from the moore-penrose pseudoinverse of the 
// gram matrix
//---------------------------------------------------------------------------------------------------//

int projector::getProjector(int const i){
  //Only apply getProjector on the site next (in direction of sweep) to the last updated
  //The gram matrix is used in constructing the projector onto the space orthogonal to the lower lying states (if any). This allows for computation of excited states.
  if(nCurrentEigen>0){
    lapack_complex_double *gram;
    lapack_int nGramEigens;
    lapack_int *suppZ;
    lapack_complex_double *gramEigenvecs;
    double *gramEigens;
    lapack_complex_double *workingMatrix, *Fki;
    lapack_complex_double zFactor;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    lapack_int info;
    gram=new lapack_complex_double[nCurrentEigen*nCurrentEigen];
    gramEigenvecs=new lapack_complex_double[nCurrentEigen*nCurrentEigen];
    suppZ=new lapack_int[2*nCurrentEigen];
    gramEigens=new double[nCurrentEigen];
    getGramMatrix(gram,i);
    lapack_int gramDim=nCurrentEigen;
    info=LAPACKE_zheevr(LAPACK_COL_MAJOR,'V','A','U',gramDim,gram,gramDim,0.0,0.0,0,0,1e-5,&nGramEigens,gramEigens,gramEigenvecs,gramDim,suppZ);
    double const minRelevantEigens=LDBL_EPSILON*nCurrentEigen*gramEigens[nGramEigens-1];
    nRelevantEigens=0;
    for(int iEigen=0;iEigen<nGramEigens;++iEigen){
      if(gramEigens[iEigen]>minRelevantEigens){
	++nRelevantEigens;
      }
    }
    getLocalDimensions(i);
    //stateArray::initialize automatically intializes all elements with zero
    auxiliaryMatrix.initialize(nRelevantEigens,lDL,lDR*ld);
    //delete[] projectionMatrix;
    //projectionMatrix=new lapack_complex_double[lDL*lDL];
    for(int mu=0;mu<nRelevantEigens;++mu){
      auxiliaryMatrix.subMatrixStart(workingMatrix,mu);
      //Eigenvectors are returned in ascending order
      for(int k=0;k<nCurrentEigen;++k){
	scalarProducts[k].F.subMatrixStart(Fki,i);
	//F^dagger * F is proportional to identity (?), but a normation is required to ensure that the projector squares to itself
	zFactor=gramEigenvecs[k+(nGramEigens-1-mu)*nCurrentEigen]*1.0/sqrt(gramEigens[nGramEigens-1-mu]);
	cblas_zaxpy(ld*lDR*lDL,&zFactor,Fki,1,workingMatrix,1);
      }
    }
    delete[] gram;
    delete[] suppZ;
    delete[] gramEigenvecs;
    delete[] gramEigens;
    if(info){
      return 1;
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// Construct the Gram matrix gram=Tr(F^dagger * F)
//---------------------------------------------------------------------------------------------------//


void projector::getGramMatrix(lapack_complex_double *gram, int const i){
  //Input has to be a nCurrentEigen x nCurrentEigen array
  if(nCurrentEigen>0){
    getLocalDimensions(i);
    lapack_complex_double simpleContainer;
    lapack_complex_double *matrixContainer, *Fki, *Fkpi;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    matrixContainer=new lapack_complex_double[ld*lDR*ld*lDR];
    for(int kp=0;kp<nCurrentEigen;++kp){
      scalarProducts[kp].F.subMatrixStart(Fkpi,i);
      for(int k=0;k<nCurrentEigen;++k){
	scalarProducts[k].F.subMatrixStart(Fki,i);
	cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,ld*lDR,ld*lDR,lDL,&zone,Fki,lDL,Fkpi,lDL,&zzero,matrixContainer,ld*lDR);
	simpleContainer=0;
	for(int mi=0;mi<ld*lDR;++mi){
	  simpleContainer+=matrixContainer[mi+ld*lDR*mi];
	}
	gram[k+kp*nCurrentEigen]=simpleContainer;
      }
    }
    delete[] matrixContainer;
  }
}
