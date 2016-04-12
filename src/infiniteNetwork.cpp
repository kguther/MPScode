#include "infiniteNetwork.h"
#include "twositeHMatrix.h"
#include "uncachedMeasurement.h"
#include <arcomp.h>
#include <arscomp.h>
#include <memory>

infiniteNetwork::infiniteNetwork(problemParameters parsIn, simulationParameters simParsIn):
  pars(parsIn),
  simPars(simParsIn),
  dimInfo(dimensionTable(simParsIn.D,2,parsIn.d))
{
  networkH=mpo<lapack_complex_double>(pars.d.maxd(),pars.Dw,pars.L);
  std::complex<int> qnValue;
  for(int iQN=0;iQN<parsIn.QNconserved.size();++iQN){
    qnValue=std::complex<int>(pars.filling[iQN]*4,pars.QNconserved[iQN].imag());
    conservedQNs.push_back(quantumNumber(dimInfo,qnValue,pars.QNLocalList[iQN]));
  }
  networkState=imps(dimInfo,conservedQNs);
  pCtr=uncachedMeasurement(&networkH,&networkState);
  diags.push_back(1.0);
}

//---------------------------------------------------------------------------------------------------//
// Main iDMRG scheme for construction of initial states
//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::growSystem(){
  int const L=pars.L;
  for(int m=0;m<L/2;++m){
    iDMRGStep();
  }
  int const lDR=dimInfo.locDimR(i);
  int const lDL=dimInfo.locDimL(i);
  arcomplex<double> *aMatrix;
  networkState.subMatrixStart(aMatrix,i);
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
  
  arcomplex<double> *bufferMatrix;
  i=(dimInfo.L()-1)/2;
  std::unique_ptr<arcomplex<double> > bufferMatrixP(new arcomplex<double>[dimInfo.locd(i)*dimInfo.locd(i+1)*dimInfo.locDimL(i)*dimInfo.locDimR(i+1)]);
  bufferMatrix=bufferMatrixP.get();
  
  //State prediction needs at least two optimized matrices, not available in the first step
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
  addSite();
}

//---------------------------------------------------------------------------------------------------//

int infiniteNetwork::optimize(arcomplex<double> *target){
  arcomplex<double> *Rterm, *Lterm;
  arcomplex<double> lambda;
  arcomplex<double> *pLambda;
  pLambda=&lambda;
  int nconv=0;
  int const maxIter=2000;
  pCtr.getLctr(Lterm);
  pCtr.getRctr(Rterm);
  twositeHMatrix BMat(Rterm,Lterm,&networkH,dimInfo,&networkState.centralIndexTable);
  BMat.prepareInput(target);
  if(BMat.dim()>1){
    ARCompStdEig<double,twositeHMatrix> eigProblemTwoSite(BMat.dim(),1,&BMat,&twositeHMatrix::MultMvBlocked,"SR",0,simPars.tolInitial,maxIter,BMat.compressedVector);
    nconv=eigProblemTwoSite.EigenValVectors(BMat.compressedVector,pLambda);
  }
  else{
    nconv=1;
  }
  BMat.readOutput(target);
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver, number of Iterations taken: "<<maxIter<<" With tolerance "<<simPars.tolInitial<<std::endl;
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::addSite(){
  pCtr.update();
  std::vector<std::complex<int> > newQNs;
  newQNs.resize(pars.QNconserved.size());
  for(int iQN=0;iQN<pars.QNconserved.size();++iQN){
    newQNs[iQN].real(static_cast<int>(2*pars.filling[iQN]*pars.L));
    newQNs[iQN].imag(pars.QNconserved[iQN].imag());
  }
  dimInfo.setParameterL(dimInfo.L()+2);
  conservedQNs[0].refine(i+1,optLocalQNs);
  networkState.addSite(dimInfo.L()+2,i+1,newQNs);
}

//---------------------------------------------------------------------------------------------------//

void infiniteNetwork::exportState(mps &target){
  networkState.exportState(target);
}
