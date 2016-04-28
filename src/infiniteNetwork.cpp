#include "infiniteNetwork.h"
#include "twositeHMatrix.h"
#include "uncachedMeasurement.h"
#include <arcomp.h>
#include <arscomp.h>
#include <memory>

#include "mkl_complex_defined.h"

#include "verifyQN.h"
#include <iostream>

//CHANGE: The iDMRG algorithm (infiniteNetwork class) shall be seperated from state management (new class)
//THEN: the infiniteNetwork shall contain an imps* networkState which is given at creation. This can then be either a imps or a timps. The state preparation class can now hand over an imps of size 10 or so. The uncached measurement class has to be updated to initialize its left/right container according to the size of the imps. Initialize diags with locD ones.

infiniteNetwork::infiniteNetwork(problemParameters const &parsIn, simulationParameters const &simParsIn, imps *MPState):
  pars(parsIn),
  simPars(simParsIn),
  firstStep(1),
  networkState(MPState)
{
  dimInfo=networkState->getDimInfo();
  networkH=mpo<lapack_complex_double>(pars.d.maxd(),pars.Dw,pars.L);
  pCtr=uncachedMeasurement(&networkH,networkState);
  //Setup of the L/R-terms 
  pCtr.getContractions(networkState->currentSite());
  int const lDR=dimInfo.locDimR(networkState->currentSite());
  diags=std::vector<double>(lDR,1.0);
  diagsm=diags;
}

//---------------------------------------------------------------------------------------------------//
// Interface functions for in/export of states.
//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::addDiags(){
  int const lDR=dimInfo.locDimR(i);
  int const lDL=dimInfo.locDimL(i);
  arcomplex<double> *aMatrix;
  networkState->subMatrixStart(aMatrix,i);
  //Multiply the center matrix into an adjacent matrix (destroying normalization) to get the full MPS into the networkState. This way, it can be copied to an mps.
  for(int si=0;si<dimInfo.locd(i);++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	aMatrix[aim+lDL*ai+si*lDL*lDR]*=diags[aim];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::iDMRGStep(){
  
  std::cout<<"Adding new sites\n";

  //In the first step, the sites are already there, this is way easier to implement since adding sites requires results for the SVs from the last step (for choice of QN Labels)
  if(firstStep==0){
    addSite();
  }
  
  arcomplex<double> *bufferMatrix;
  i=(dimInfo.L()-1)/2;
  std::unique_ptr<arcomplex<double> > bufferMatrixP(new arcomplex<double>[dimInfo.locd(i)*dimInfo.locd(i+1)*dimInfo.locDimL(i)*dimInfo.locDimR(i+1)]);
  bufferMatrix=bufferMatrixP.get();

  //State prediction needs at least two optimized matrices, not available for a two-site system
  if(i!=0){
    statePrediction(bufferMatrix);
  }
  else{
    for(int m=0;m<dimInfo.locd(i)*dimInfo.locd(i+1);++m){
      bufferMatrix[m]=1;
    }
  }
  optimize(bufferMatrix);
  updateMPS(bufferMatrix);

  firstStep=0;
}

//---------------------------------------------------------------------------------------------------//

int infiniteNetwork::optimize(arcomplex<double> *target){
  //Twosite optimization
  arcomplex<double> *Rterm, *Lterm;
  arcomplex<double> lambda;
  arcomplex<double> *pLambda;
  pLambda=&lambda;
  int nconv=0;
  int const maxIter=2000;
  //In the first step, the endpoints 0 and 1 of the MPO are used, in all other steps, the middle part is used. Due to translation invariance, we always use the site matrices of i==1 and i==2.
  int const HPos=(i==0)?0:1;
  pCtr.getLctr(Lterm);
  pCtr.getRctr(Rterm);
  twositeHMatrix BMat(Rterm,Lterm,&networkH,HPos,dimInfo,&(networkState->centralIndexTable()));

  //Compression and expansion have to be done manually
  BMat.prepareInput(target);
  
  if(BMat.dim()>1){
    arcomplex<double> *compressedVector=BMat.getCompressedVector();
    ARCompStdEig<double,twositeHMatrix> eigProblemTwoSite(BMat.dim(),1,&BMat,&twositeHMatrix::MultMvBlocked,"SR",0,simPars.tolInitial,maxIter,compressedVector);
    nconv=eigProblemTwoSite.EigenValVectors(compressedVector,pLambda);
  }
  else{
    nconv=1;
  }
  BMat.readOutput(target);
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver, number of Iterations taken: "<<maxIter<<" With tolerance "<<simPars.tolInitial<<std::endl;
    return 1;
  }
  else{
    std::cout<<"Current energy: "<<real(lambda)<<std::endl;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::addSite(){
  //Can only be used after optimization since it relies on the outcome of optimize() 
  pCtr.update();
  dimInfo.setParameterL(dimInfo.L()+2);

  //QNs have to be adjusted to fit a larger system
  std::vector<std::complex<int> > newQNs;
  newQNs.resize(pars.QNconserved.size());
  for(int iQN=0;iQN<pars.QNconserved.size();++iQN){
    newQNs[iQN].real(static_cast<int>(pars.filling[iQN]*dimInfo.L()));
    newQNs[iQN].imag(pars.QNconserved[iQN].imag());
  }
  std::cout<<"New global QN: "<<newQNs[0]<<std::endl;
  std::cout<<"Obtained from system length "<<dimInfo.L()<<" and filling "<<pars.filling[0]<<std::endl;
  
  //Only for a single QN
  
  int info;

  //The QNs determined in updateMPS are stored
  info=networkState->refineQN(i+1,optLocalQNsL);

  //This is where the actual MPS is grown
  networkState->addSite(dimInfo.L(),i+1,newQNs,optLocalQNsR);

}

//---------------------------------------------------------------------------------------------------//

imps* infiniteNetwork::getState(){
  return networkState;
}
