#include <float.h>
#include <memory>
#include "projector.h"
#include "mps.h"
#include "overlap.h"
#include "arrayprocessing.h"
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// The projector class is still experimental, it seems to yield correct results, but I am not sure
// whether everything is correct.
//---------------------------------------------------------------------------------------------------//

projector::projector(int nEigsin):
  nEigs(nEigsin),
  nCurrentEigen(0)
{
  orthoStates=new mps[nEigs];
  scalarProducts=new overlap[nEigs-1];
}

//---------------------------------------------------------------------------------------------------//

projector::projector(projector const &source){
  pCpy(source);
}

//---------------------------------------------------------------------------------------------------//

projector::~projector(){
  delete[] orthoStates;
  delete[] scalarProducts;
}

//---------------------------------------------------------------------------------------------------//

projector& projector::operator=(projector const &source){
  pCpy(source);
  return *this;
}

//---------------------------------------------------------------------------------------------------//

void projector::pCpy(projector const &source){
  if(source.nEigs!=nEigs){
    delete[] orthoStates;
    delete[] scalarProducts;
    nEigs=source.nEigs;
    orthoStates=new mps[nEigs];
    scalarProducts=new overlap[nEigs-1];
  }
  for(int iEigen=0;iEigen<nEigs-1;++iEigen){
    orthoStates[iEigen]=source.orthoStates[iEigen];
    scalarProducts[iEigen]=source.scalarProducts[iEigen];
  }
  orthoStates[nEigs-1]=source.orthoStates[nEigs-1];
}

//---------------------------------------------------------------------------------------------------//
// When adjusting the simulation parameter D, all stored states are updated
//---------------------------------------------------------------------------------------------------//

void projector::setParameterD(int Dnew){
  for(int iEigen=0;iEigen<nEigs;++iEigen){
    orthoStates[iEigen].setParameterD(Dnew);
  }
}

//---------------------------------------------------------------------------------------------------//

void projector::setParameterL(int Lnew){
  for(int iEigen=0;iEigen<nEigs;++iEigen){
    orthoStates[iEigen].setParameterL(Lnew);
  }
}

//---------------------------------------------------------------------------------------------------//
// This functions loads the state which shall be orthogonalized, and uses its index to determine to
// which states it shall be orthogonalized.
//---------------------------------------------------------------------------------------------------//

void projector::loadScalarProducts(mps const &variationalState, int iEigen){
  for(int k=0;k<iEigen;++k){
    scalarProducts[k].loadMPS(&variationalState,&orthoStates[k]);
  }
  nCurrentEigen=iEigen;
}

//---------------------------------------------------------------------------------------------------//

void projector::updateScalarProducts(int i, int direction){
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


void projector::storeOrthoState(mps const &source, int iEigen){
  orthoStates[iEigen]=source;
}

//---------------------------------------------------------------------------------------------------//

void projector::storeCurrentState(mps const &source){
  orthoStates[nCurrentEigen]=source;
}

//---------------------------------------------------------------------------------------------------//

void projector::getStoredState(mps *&target, int iEigen){
  target=&(orthoStates[iEigen]);
}

//---------------------------------------------------------------------------------------------------//

int projector::loadNextState(mps &target, int iEigen){
  target=orthoStates[iEigen];
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int projector::loadNextState(mps &target){
  if(nCurrentEigen<(nEigs-1)){
    target=orthoStates[nCurrentEigen+1];
    return 0;
  }
  return 1;
}		

//---------------------------------------------------------------------------------------------------//

void projector::getLocalDimensions(int i){
  //Note that it is required for all stored states to have the same bond dimension
  if(nEigs>0){
    ld=orthoStates[0].locd(i);
    lDL=orthoStates[0].locDimL(i);
    lDR=orthoStates[0].locDimR(i);
  }
}

//---------------------------------------------------------------------------------------------------//

lapack_complex_double projector::fullOverlap(int k){
  return scalarProducts[k].fullOverlap();
}

//---------------------------------------------------------------------------------------------------//
// Although it is self-explanatory what project does, keep in mind that one has to use getProjector
// for the current site before using project to ensure the auxiliary matrix is set.
//---------------------------------------------------------------------------------------------------//

void projector::project(lapack_complex_double *vec, int i){
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
      auxiliaryMatrix.getPtr(G,mu);
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
    delete[] vecContainer;
    delete[] trContainer;
  }
}

//---------------------------------------------------------------------------------------------------//
// Compute the auxiliary matrices of the projector from the moore-penrose pseudoinverse of the 
// gram matrix
//---------------------------------------------------------------------------------------------------//

int projector::getProjector(int i){
  //Only apply getProjector on the site next (in direction of sweep) to the last updated
  //The gram matrix is used in constructing the projector onto the space orthogonal to the lower lying states (if any). This allows for computation of excited states.
  if(nCurrentEigen>0){
    lapack_complex_double *gram;
    lapack_int nGramEigens;
    lapack_complex_double *gramEigenvecs;
    double *gramEigens;
    lapack_complex_double *workingMatrix, *Fki;
    lapack_complex_double zFactor;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    lapack_int info;
    std::unique_ptr<lapack_complex_double[]> gramP(new lapack_complex_double[nCurrentEigen*nCurrentEigen]);
    std::unique_ptr<lapack_complex_double[]> gramEigenvecsP(new lapack_complex_double[nCurrentEigen*nCurrentEigen]);
    std::unique_ptr<double> gramEigensP(new double[nCurrentEigen]);
    gram=gramP.get();
    gramEigenvecs=gramEigenvecsP.get();
    gramEigens=gramEigensP.get();
    getGramMatrix(gram,i);
    lapack_int const gramDim=nCurrentEigen;
    if(nCurrentEigen>1){
      lapack_int *suppZ;
      std::unique_ptr<lapack_int> suppZP(new lapack_int[2*nCurrentEigen]);
      suppZ=suppZP.get();
      info=LAPACKE_zheevr(LAPACK_COL_MAJOR,'V','A','U',gramDim,gram,gramDim,0.0,0.0,0,0,1e-5,&nGramEigens,gramEigens,gramEigenvecs,gramDim,suppZ);
    }
    else{
      gramEigenvecs[0]=1;
      gramEigens[0]=real(gram[0]);
      info=0;
      nGramEigens=1;
    }
    if(info){
      std::cout<<"Error in LAPACKE_zheevr: "<<info<<" at site "<<i<<std::endl;
      return 1;
    }
    double const minRelevantEigens=LDBL_EPSILON*nCurrentEigen*gramEigens[nGramEigens-1];
    nRelevantEigens=0;
    for(int iEigen=0;iEigen<nGramEigens;++iEigen){
      if(gramEigens[iEigen]>minRelevantEigens){
	++nRelevantEigens;
      }
    }
    getLocalDimensions(i);
    //stateArray::initialize automatically intializes all elements with zero
    std::vector<int> dimensions(3,0);
    dimensions[0]=nRelevantEigens;
    dimensions[1]=lDL;
    dimensions[2]=lDR*ld;
    auxiliaryMatrix=baseTensor<lapack_complex_double>(dimensions);
    //delete[] projectionMatrix;
    //projectionMatrix=new lapack_complex_double[lDL*lDL];
    for(int mu=0;mu<nRelevantEigens;++mu){
      auxiliaryMatrix.getPtr(workingMatrix,mu);
      //Eigenvectors are returned in ascending order
      for(int k=0;k<nCurrentEigen;++k){
	scalarProducts[k].F.subMatrixStart(Fki,i);
	//F^dagger * F is proportional to identity (?), but a normation is required to ensure that the projector squares to itself
	zFactor=gramEigenvecs[k+(nGramEigens-1-mu)*nCurrentEigen]*1.0/sqrt(gramEigens[nGramEigens-1-mu]);
	cblas_zaxpy(ld*lDR*lDL,&zFactor,Fki,1,workingMatrix,1);
      }
    }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// Construct the Gram matrix gram=Tr(F^dagger * F)
//---------------------------------------------------------------------------------------------------//


void projector::getGramMatrix(lapack_complex_double *gram, int i){
  //Input has to be a nCurrentEigen x nCurrentEigen array
  if(nCurrentEigen>0){
    getLocalDimensions(i);
    lapack_complex_double simpleContainer;
    lapack_complex_double *matrixContainer, *Fki, *Fkpi, *FkiT;
    lapack_complex_double zone=1.0;
    lapack_complex_double zzero=0.0;
    std::unique_ptr<lapack_complex_double[]> FkiTP(new lapack_complex_double[ld*lDR*lDL]);
    FkiT=FkiTP.get();
    std::unique_ptr<lapack_complex_double[]> matrixContainerP(new lapack_complex_double[ld*lDR*ld*lDR]);
    matrixContainer=matrixContainerP.get();
    for(int kp=0;kp<nCurrentEigen;++kp){
      scalarProducts[kp].F.subMatrixStart(Fkpi,i);
      for(int k=0;k<nCurrentEigen;++k){
	scalarProducts[k].F.subMatrixStart(Fki,i);
	auxiliary::arraycpy(ld*lDR*lDL,Fki,FkiT);
	auxiliary::transp(ld*lDR,lDL,FkiT);
	for(int m=0;m<ld*lDR*lDL;++m){
	  FkiT[m]=conj(FkiT[m]);
	}
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ld*lDR,ld*lDR,lDL,&zone,FkiT,ld*lDR,Fkpi,lDL,&zzero,matrixContainer,ld*lDR);
	simpleContainer=0;
	for(int mi=0;mi<ld*lDR;++mi){
	  simpleContainer+=matrixContainer[mi+ld*lDR*mi];
	}
	gram[k+kp*nCurrentEigen]=simpleContainer;
      }
    }
  }
}
