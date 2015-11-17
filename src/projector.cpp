#include <cblas.h>
#include <lapacke.h>
#include <float.h>
#include "projector.h"
#include "siteArray.h"
#include "mps.h"
#include "overlap.h"

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
  scalarProducts=new overlap[((nEigs-1)*nEigs)/2];
}

//---------------------------------------------------------------------------------------------------//

void projector::setParameterD(int const Dnew){
  for(int iEigen=0;iEigen<nEigs;++iEigen){
    orthoStates[iEigen].setParameterD(Dnew);
  }
}

//---------------------------------------------------------------------------------------------------//

void projector::loadScalarProducts(mps *variationalState,int const iEigen){
  //offset marks the beginning of the scalar products with the current state
  int const offset=(iEigen*(iEigen-1))/2;
  for(int k=0;k<iEigen;++k){
    scalarProducts[k+offset].loadMPS(&orthoStates[k],variationalState);
  }
}

//---------------------------------------------------------------------------------------------------//

void projector::updateScalarProducts(int const i, int const direction){
  if(nCurrentEigen>0){
    int const offset=(nCurrentEigen*(nCurrentEigen-1))/2;
    for(int k=0;k<nCurrentEigen;++k){
      if(direction==1){
	scalarProducts[k+offset].stepRight(i);
      }
      else{
	scalarProducts[k+offset].stepLeft(i);
      }
    }
  }
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
// Although it is self-explanatory what project does, keep in mind that one has to use getProjector
// for the current site before using project to ensure the auxiliary matrix is set.
//---------------------------------------------------------------------------------------------------//

void projector::project(int const i, void *vec){
  lapack_complex_double *trContainer;
  lapack_complex_double *G;
  lapack_complex_double *vecContainer;
  lapack_complex_double simpleContainer;
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  getLocalDimensions(currentSite);
  trContainer=new lapack_complex_double[ld*lDR*ld*lDR];
  vecContainer=new lapack_complex_double[ld*lDR*lDL];
  for(int mi=0;mi<ld*lDR*lDL;++mi){
    vecContainer[mi]=0;
  }
  for(int mu=0;mu<nRelevantEigens;++mu){
    auxiliaryMatrix.subMatrixStart(G,mu);
    cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,ld*lDR,ld*lDR,lDL,&zone,G,lDL,vec,lDL,&zzero,trContainer,ld*lDR);
    simpleContainer=0;
    for(int mi=0;mi<ld*lDR;++mi){
      simpleContainer+=trContainer[mi+ld*lDR*mi];
    }
    simpleContainer*=-1.0;
    cblas_zaxpy(ld*lDR*lDL,&simpleContainer,G,1,vecContainer,1);
  }
  for(int mi=0;mi<ld*lDR*lDL;++mi){
    vec[mi]+=vecContainer[mi];
  }
  delete[] vecContainer;
  delete[] trContainer;
}

//---------------------------------------------------------------------------------------------------//

void projector::getProjector(int const i){
  //The gram matrix is used in constructing the projector onto the space orthogonal to the lower lying states (if any). This allows for computation of excited states.
  lapack_complex_double *gram;
  gram=new lapack_complex_double[nCurrentEigen*nCurrentEigen];
  getGramMatrix(gram,i);
  lapack_int nGramEigens;
  lapack_int *suppZ;
  lapack_complex_double *gramEigenvecs;
  double *gramEigens;
  gramEigenvecs=new lapack_complex_double[nCurrentEigen*nCurrentEigen];
  suppZ=new lapack_int[2*nCurrentEigen];
  gramEigens=new double[nCurrentEigen];
  LAPACKE_zheevr(LAPACK_COL_MAJOR,'V','A','U',nCurrentEigen,gram,nCurrentEigen,0.0,0.0,0,0,1e-5,&nGramEigens,gramEigens,gramEigenvecs,nCurrentEigen,suppZ);
  int const minRelevantEigens=LDBL_EPSILON*nCurrentEigen*gramEigens[nGramEigens-1];
  nRelevantEigens=0;
  for(int iEigen=0;iEigen<nGramEigens;++iEigen){
    if(gramEigens[iEigen]>minRelevantEigens){
      ++nRelevantEigens;
    }
  }
  getLocalDimensions(i);
  auxiliaryMatrix.initialize(nRelevantEigens,lDL,ld*lDR);
  lapack_complex_double *workingMatrix, *Fki;
  lapack_complex_double zFactor;
  int const offset=(nCurrentEigen*(nCurrentEigen-1))/2;
  for(int mu=0;mu<nRelevantEigens;++mu){
    auxiliaryMatrix.subMatrixStart(workingMatrix,mu);
    //Eigenvectors are returned in ascending order
    for(int k=0;k<nCurrentEigen;++k){
      scalarProducts[k+offset].F.subMatrixStart(Fki,i);
      zFactor=1.0/sqrt(gramEigens[nGramEigens-1-mu])*gramEigenvecs[k+mu*nCurrentEigen];
      cblas_zaxpy(ld*lDR*lDL,&zFactor,Fki,1,workingMatrix,1);
    }
  } 
  delete[] gram;
  delete[] suppZ;
  delete[] gramEigenvecs;
  delete[] gramEigens;
}

//---------------------------------------------------------------------------------------------------//

void projector::getGramMatrix(lapack_complex_double *gram, int const i){
  //Input has to be a nCurrentEigen x nCurrentEigen array
  if(nCurrentEigen>0){
    int const offset=(nCurrentEigen*(nCurrentEigen-1))/2;
    getLocalDimensions(i);
    lapack_complex_double simpleContainer;
    lapack_complex_double *matrixContainer, *Fki, *Fkpi;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    matrixContainer=new lapack_complex_double[ld*lDR*ld*lDR];
    for(int kp=0;kp<nCurrentEigen;++kp){
      scalarProducts[kp+offset].F.subMatrixStart(Fkpi,i);
      for(int k=0;k<nCurrentEigen;++k){
	scalarProducts[k+offset].F.subMatrixStart(Fki,i);
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
