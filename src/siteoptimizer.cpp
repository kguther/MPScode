#include <arcomp.h>
#include <iostream>
#include <arlnsmat.h>
#include <arlscomp.h>
#include "siteoptimizer.h"

siteoptimizer::siteoptimizer(int dimensionin, int nzelin, arcomplex<double> *Hin, int *irowin, int *pcolin){
  H=Hin;
  irow=irowin;
  pcol=pcolin;
  dimension=dimensionin;
  nzel=nzelin;
}

//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int siteoptimizer::solveEigen(arcomplex<double> *plambda, arcomplex<double> *currentM){
  int nconv;
  ARluNonSymMatrix<arcomplex<double> ,double> Hopt(dimension,nzel,H,irow,pcol);
  ARluCompStdEig<double> OptProblem(1,Hopt,"SR",3,1e-8);
  nconv=OptProblem.EigenValVectors(currentM,plambda);
  if(nconv==0){
    std::cout<<"Failed to converge in iterative eigensolver";
    return 1;
  }
  return 0;
}
