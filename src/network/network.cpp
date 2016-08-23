#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <arscomp.h>
#include <chrono>
#include "network.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "blockHMatrix.h"
#include "globalMeasurement.h"
#include "localMeasurementSeries.h"
#include "exactGroundState.h"
#include "exceptionClasses.h"

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//NAMING CONVENTION: i is always the site index. ALWAYS.

//---------------------------------------------------------------------------------------------------//
// 'Main' file containing the MPS ansatz solver itself (this is not the actual main.cpp, just the
// most important file)
//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
// Constructors and similar stuff for the basic network class
//---------------------------------------------------------------------------------------------------//

network::network()
  //Empty networks dont do much, use assignment to fill them
{}

//---------------------------------------------------------------------------------------------------//

network::network(problemParameters const &inputpars, simulationParameters const &inputsimPars):
  pars(inputpars),
  D(inputsimPars.D),
  L(inputpars.L),
  Dw(inputpars.Dw),
  simPars(inputsimPars),
  excitedStateP(projector(pars.nEigs)),
  networkDimInfo(dimensionTable(D,L,pars.d)),
  reSweep(0),
  check(0),
  checkParity(0)
{
  networkH=mpo<lapack_complex_double>(pars.d.maxd(),Dw,L);
  nConverged.resize(pars.nEigs);
  for(int iEigen=0;iEigen<nConverged.size();++iEigen){
    nConverged[iEigen]=1;
  }
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    conservedQNs.push_back(quantumNumber(networkDimInfo,pars.QNconserved[iQN],pars.QNLocalList[iQN]));
  }
  resetState();
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary methods for management of multiple runs
//---------------------------------------------------------------------------------------------------//

void network::loadNetworkState(mps const &source){
  networkState=source;
}

//---------------------------------------------------------------------------------------------------//

void network::exportNetworkState(mps &target){
  target=networkState;
}

//---------------------------------------------------------------------------------------------------//

void network::resetConvergence(){
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    nConverged[iEigen]=1;
  }
}

//---------------------------------------------------------------------------------------------------//

void network::resetState(){
  //Generation of MPS, that is basically the index tables for the QNs are written
  networkState.generate(networkDimInfo,conservedQNs);
  //All states have to be stored, for reusability in later solve() with other simulation parameters
  //Somewhat unelegant way to handle loading of the stored states in the first solve()
  for(int iEigen=1;iEigen<pars.nEigs;++iEigen){
    excitedStateP.storeOrthoState(networkState,iEigen);
  }
  if(pars.nQNs && pars.d.maxd()==4 && networkH.maxDim()==12){
    if(conservedQNs[0].QNValue().imag()){  
      exactGroundState gsLoader(conservedQNs[0].QNValue());
      gsLoader.writeExactGroundState(networkState);
    }
    //In the absence of Z2-symmetry, the ground state guess is the even exact ground state, the first excited state guess is the odd exact ground state
    else{
      //the Z2-number is set to 1 internally in exactGroundState if the input i 0
      std::complex<int> even=conservedQNs[0].QNValue();
      std::complex<int> odd=std::complex<int>(conservedQNs[0].QNValue().real(),-1);
      exactGroundState gsLoaderA(even);
      exactGroundState gsLoaderB(odd);
      gsLoaderA.writeExactGroundState(networkState);
      if(pars.nEigs>1){
	std::vector<quantumNumber> conjugateQNs;
	for(int iQN=0;iQN<pars.nQNs;++iQN){
	  conjugateQNs.push_back(quantumNumber(networkDimInfo,pars.QNconserved[iQN],pars.QNLocalList[iQN],-1));
	}
	mps firstInitialState(networkDimInfo,conjugateQNs);
	gsLoaderB.writeExactGroundState(firstInitialState);
	excitedStateP.storeOrthoState(firstInitialState,1);
      }
    }  
  }
  excitedStateP.storeOrthoState(networkState,0);
}

//---------------------------------------------------------------------------------------------------//
// These functions can be employed to alter the algorithm parameters nSweeps and D during lifetime of a 
// network object. This allows for iteratively increasing D. They completely take care of all required
// updated within the network, but they should be needed only a few times.
//---------------------------------------------------------------------------------------------------//

int network::setSimParameters(simulationParameters const &newPars){
  int info;
  simPars=newPars;
  info=setParameterD(simPars.D);
  if(info){
    return -1;
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

int network::setParameterD(int Dnew){
  //Decreasing D is neither required nor reasonable in any context
  if(Dnew<=D){
    return -1;
  }
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    conservedQNs[iQN].setParameterD(Dnew);
  }
  //this is problematic as scaling D requires a new QN labeling scheme -> translation between labels is required
  networkState.setParameterD(Dnew);
  networkDimInfo.setParameterD(Dnew);
  //All stored states have to be brought into the correct form for compatibility with current D
  excitedStateP.setParameterD(Dnew);
  //Adapt D
  D=Dnew;
  simPars.D=Dnew;
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void network::getLocalDimensions(int i){
  lDL=networkDimInfo.locDimL(i);
  lDR=networkDimInfo.locDimR(i);
  ld=networkDimInfo.locd(i);
  lDwR=networkH.locDimR(i);
  lDwL=networkH.locDimL(i);
}

//---------------------------------------------------------------------------------------------------//
// This is the main solver, using the networkState as initial and final state and networkH as MPO
// representation of the Hamiltonian. Computed are eigenstate (initial state is overwritten) and 
// corresponding eigenvalue (output).
//---------------------------------------------------------------------------------------------------//

int network::solve(std::vector<double> &lambda, std::vector<double> &deltaLambda){
  std::chrono::steady_clock::time_point t1=std::chrono::steady_clock::now();
  int maxIter=10000;
  int stepRet;
  int cshift=0;
  double tol;
  double spinCheck=0;
  double parCheck=0;
  lambda.resize(pars.nEigs);
  deltaLambda.resize(pars.nEigs);
  alpha=simPars.alpha;
  if(pars.nQNs || pars.nEigs>1){
    cshift=-100;
  }
  //Necessary if multiple instances of solve() are called
  excitedStateP.loadNextState(networkState,0);
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    pCtr.initialize(&networkH,&networkState);
    std::cout<<"Starting normalization\n";
    //For more than one QN, the static scheme does not necessarily yield a valid labeling scheme
    int const sec=(pars.nQNs>1)?1:0;
    double alphaBuf=alpha;
    alpha=0.0;
    for(int i=L-1;i>0;--i){
      normalize(i,0,sec);
    }
    alpha=alphaBuf;
    networkState.normalizeFinal(1);
#ifdef QNCHECK
    overlap test;
    test.loadMPS(&networkState,&networkState);
    std::cout<<"Norm: "<<test.getFullOverlap()<<std::endl;
      if(check){
	measure(check,spinCheck);
	std::cout<<"Current particle number (initial): "<<spinCheck<<std::endl;
      }
      if(checkParity){
	measure(checkParity,parCheck);
	std::cout<<"Current subchain parity (initial): "<<parCheck<<std::endl;
      }
#endif

    std::cout<<"Computing partial contractions\n";
    pCtr.Lctr.global_access(0,0,0,0)=1;
    //In preparation of the first sweep, generate full contraction to the right (first sweeps starts at site 0)
    pCtr.calcCtrFull(1);

    shift=cshift;
    std::cout<<"Computing state "<<iEigen<<std::endl;
    alpha=simPars.alpha;
    tol=simPars.tolInitial;
    //load all pairings with the current state and previous ones into the scalar products
    if(iEigen){
      std::cout<<"GENERATING PROJECTOR"<<std::endl;
    }
    excitedStateP.loadScalarProducts(networkState,iEigen);

    pCtr.calcCtrIterRightBase(-1,&expectationValue);
    std::cout<<"Initial energy: "<<expectationValue<<std::endl;

    lambda[iEigen]=real(expectationValue)+shift;
    for(int iSweep=0;iSweep<simPars.nSweeps;++iSweep){
      if(!nConverged[iEigen]){
	break;
      }
      //actual sweep is executed here
      try{
	sweep(maxIter,tol,lambda[iEigen]);
	//In calcCtrIterRightBase, the second argument has to be a pointer, because it usually is an array. No call-by-reference here, using a container as argument would also just make things complicated
	pCtr.calcCtrIterRightBase(-1,&expectationValue);

	if(tol>simPars.tolMin){
	  tol*=pow(simPars.tolMin/simPars.tolInitial,1.0/simPars.nSweeps);
	}
	deltaLambda[iEigen]=1;
	if(iSweep==simPars.nSweeps-1){
	  std::cout<<"Calculating quality of convergence\n";
	  deltaLambda[iEigen]=convergenceCheck();
	  if(deltaLambda[iEigen]<simPars.devAccuracy){
	    nConverged[iEigen]=0;
	  }
	  std::cout<<"Quality of convergence: "<<deltaLambda[iEigen]<<"\tRequired accuracy: "<<simPars.devAccuracy<<std::endl<<std::endl;
	  std::cout<<"Used bond dimension: "<<D<<std::endl;
	}
	//the last optimization has not been included in the projector's F-Matrix yet
	excitedStateP.updateScalarProducts(0,-1);
	for(int prev=0;prev<iEigen;++prev){
	  std::cout<<"Overlap with state "<<prev<<" is: "<<excitedStateP.fullOverlap(prev)<<std::endl;
	}
      }
      catch(svd_failure &err){
	if(reSweep){
	  throw critical_error();
	}
	resetSweep();
      }

      if(check){
	measure(check,spinCheck);
	std::cout<<"Current particle number (final): "<<spinCheck<<std::endl;
      }
      if(checkParity){
	measure(checkParity,parCheck);
	std::cout<<"Current subchain parity (final): "<<parCheck<<std::endl;
      }

    }
    stepRet=gotoNextEigen();
    if(!stepRet){
      std::cout<<"LOADED STATES. PREPARED COMPUTATION OF NEXT EIGENSTATE"<<std::endl;
    }
  }
  std::chrono::duration<double> deltaT=std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now()-t1);
  std::cout<<"Calculation took "<<deltaT.count()<<" seconds\n\n";
  lambda[0]-=cshift;
  for(int iEigen=1;iEigen<pars.nEigs;++iEigen){
    lambda[iEigen]-=shift;
  }
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    if(nConverged[iEigen]){
      return 1;
    }
  }    
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void network::sweep(double maxIter, double tol, double &lambda){
  // By default enrichment is used whenever conserved QNs are used
  int const expFlag=pars.nQNs;
  double lambdaCont;
  
  //Initializing the overlap takes some time, thus, it is disabled
  //double testE;
  //overlap test(&networkState,&networkState);
  std::cout<<"STARTING RIGHTSWEEP\n\n";
  for(int i=0;i<(L-1);++i){
    //Step of leftsweep
    std::cout<<"Optimizing site matrix"<<std::endl;
    lambdaCont=lambda;
    std::chrono::steady_clock::time_point t1=std::chrono::steady_clock::now();
    optimize(i,maxIter,tol,lambda);
    std::chrono::duration<double> deltaT=std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now()-t1);
    std::cout<<"Optimization took "<<deltaT.count()<<" seconds\n\n";
    //Execute left-sided enrichment step and update the coefficient of the expansion term
    normalize(i,1,expFlag);
    //Here, the scalar products with lower lying states are updated
    excitedStateP.updateScalarProducts(i,1);
    pCtr.calcCtrIterLeft(i+1);  
    if(expFlag){
      getNewAlpha(i+1,lambda,lambdaCont);
    }  
  }

  networkState.normalizeFinal(0);
  std::cout<<"STARTING LEFTSWEEP\n\n";
  for(int i=L-1;i>0;--i){
    //Step of rightsweep  
    std::cout<<"Optimizing site matrix"<<std::endl;    
    lambdaCont=lambda;
    std::chrono::steady_clock::time_point t1=std::chrono::steady_clock::now();
    optimize(i,maxIter,tol,lambda);
    std::chrono::duration<double> deltaT=std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now()-t1);
    std::cout<<"Optimization took "<<deltaT.count()<<" seconds\n\n";
    //Execute right-sided enrichment step and update the coefficient of the expansion term
    normalize(i,0,expFlag);
    //same as above for the scalar products with lower lying states
    excitedStateP.updateScalarProducts(i,-1);
    pCtr.calcCtrIterRight(i-1);
    if(expFlag){
      getNewAlpha(i-1,lambda,lambdaCont);
    }
  }
  networkState.normalizeFinal(1);
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int i, int maxIter, double tol, double &iolambda){
  //Invokes ARPACK++ to solve the eigenvalue problem
  std::complex<double> lambda;
  std::complex<double> *plambda;
  std::complex<double> *currentM;
  std::complex<double> *RTerm, *LTerm;
  void (optHMatrix::*multMv)(std::complex<double> *v, std::complex<double> *w);
  int nconv;
  //Get the projector onto the space orthogonal to any lower lying states
  //Get the current partial contractions and site matrix of the Hamiltonian
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  nconv=excitedStateP.getProjector(i);
  if(nconv){
    return -1;
  }
  plambda=&lambda;
  //Using the current site matrix as a starting point allows for much faster convergence as it has already been optimized in previous sweeps (except for the first sweep, this is where a good starting point has to be guessed)
  networkState.subMatrixStart(currentM,i);

  //Check step useful whenever something in the normalization or optimization is adjusted
#ifdef QNCHECK
  double spinCheck=0;
  double parCheck=0;
      if(check){
	measure(check,spinCheck);
	std::cout<<"Current particle number (opt): "<<spinCheck<<std::endl;
      }
      if(checkParity){
	measure(checkParity,parCheck);
	std::cout<<"Current subchain parity (opt): "<<parCheck<<std::endl;
      }
#endif
    
  double lambdaCont;
  if(pars.nQNs){// && i!=0 && i!=(L-1)){
    //For some obscure reason, ARPACK++ could previusly not handle the boundary problems with reduced dimension. 
    //They have to be solved without using the block structure. Since they have a really tiny dimension, this does not matter at all. 
    //If there is some problem with arpack at first/last site, disable QN optimized solution on those sites
    blockHMatrix BMat(RTerm, LTerm,&networkH,networkDimInfo,Dw,i,&(networkState.indexTable().getLocalIndexTable(i)),&excitedStateP,shift,&conservedQNs);
    BMat.prepareInput(currentM);
    std::complex<double> *compressedVector=BMat.getCompressedVector();
    if(BMat.dim()>1){
      ARCompStdEig<double, blockHMatrix> eigProblemBlocked(BMat.dim(),1,&BMat,&blockHMatrix::MultMvBlocked,"SR",0,tol,maxIter,compressedVector);
      nconv=eigProblemBlocked.EigenValVectors(compressedVector,plambda);
    std::cout<<"Number of iterations taken: "<<eigProblemBlocked.GetIter()<<std::endl;
    }
    if(nconv){
      //if arpack does not converge, the result is likely even worse than the input (empirical statement)
      BMat.readOutput(currentM);
    }
  }
  else{
    if(pars.nQNs){
      multMv=&optHMatrix::MultMvQNConserving;
    }
    else{
      multMv=&optHMatrix::MultMv;
    }
    //Generate matrix which is to be passed to ARPACK++
    optHMatrix HMat(RTerm,LTerm,&networkH,networkDimInfo,Dw,i,&excitedStateP,shift,&conservedQNs);
    //Note that the types given do and have to match the ones in the projector class if more than one eigenvalue is computed
    ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,multMv,"SR",0,tol,maxIter,currentM);
    //One should avoid to hit the maximum number of iterations since this can lead into a suboptimal site matrix, increasing the current energy (although usually not by a lot)
    nconv=eigProblem.EigenValVectors(currentM,plambda);
  }

  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver, number of Iterations taken: "<<maxIter<<" With tolerance "<<tol<<std::endl;
    return 1;
  }
  else{
    iolambda=real(lambda);
    std::cout<<std::setprecision(21)<<"Current energy: "<<real(lambda)-shift<<std::endl;
    return 0;
  }
}

//---------------------------------------------------------------------------------------------------//
// Functions for normalization of site matrices, creation of canonical MPS and adressing the subspace
// expansion
//---------------------------------------------------------------------------------------------------//

void network::normalize(int i, int direction, int enrichment){
  if(direction){
    if(enrichment){
      if(pars.nQNs){
	leftEnrichmentBlockwise(i);
      }
      else{
	leftEnrichment(i);
      }
    }
    else{
      networkState.leftNormalizeState(i);
    }
  }
  else{
    if(enrichment){
      if(pars.nQNs){
	rightEnrichmentBlockwise(i);
      }
      else{
	rightEnrichment(i);
      }
    }
    else{
      networkState.rightNormalizeState(i);
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void network::resetSweep(){
  double const alphaBuf=alpha;
  alpha=0.0;
  for(int j=L;j>0;--j){
    try{
      normalize(j,0,pars.nQNs);
    }
    catch(svd_failure &err){
      throw critical_error();
    }
  }
  alpha=alphaBuf;
  networkState.normalizeFinal(1);
  pCtr.Lctr.global_access(0,0,0,0)=1;
  pCtr.calcCtrFull(1);
  reSweep=1;
}

//---------------------------------------------------------------------------------------------------//
// These functions compute the standard deviation of Energy. This is used to determine the quality
// of convergence after each sweep and replaces the old-school DMRG error estimates.
// We can access the measurement class interface therefore.
//---------------------------------------------------------------------------------------------------//

double network::convergenceCheck(){
  double stdDeviation;
  double meanEnergy;
  double meanSqrEnergy;
  meanEnergy=real(expectationValue);
  meanEnergy=pow(meanEnergy,2);
  calcHSqrExpectationValue(meanSqrEnergy);
  stdDeviation=meanSqrEnergy-meanEnergy;
  std::cout<<"Current quality of convergence: "<<stdDeviation<<std::endl;
  //stdDeviation is always positive in theory.
  return std::abs(stdDeviation);
}

//---------------------------------------------------------------------------------------------------//

void network::calcHSqrExpectationValue(double &ioHsqr){
  mpo<lapack_complex_double> Hsqr(networkState.siteDim(),Dw*Dw,L);
  lapack_complex_double simpleContainer;
  for(int i=0;i<L;++i){
    getLocalDimensions(i);
    for(int si=0;si<ld;++si){
      for(int sip=0;sip<ld;++sip){
	for(int bi=0;bi<lDwR;++bi){
	  for(int bip=0;bip<lDwR;++bip){
	    for(int bim=0;bim<lDwL;++bim){
	      for(int bimp=0;bimp<lDwL;++bimp){
		simpleContainer=0;
		for(int sipp=0;sipp<ld;++sipp){
		  simpleContainer+=networkH.global_access(i,si,sipp,bi,bim)*networkH.global_access(i,sipp,sip,bip,bimp);
		}
		Hsqr.global_access(i,si,sip,bi+lDwR*bip,bim+lDwL*bimp)=simpleContainer;
	      }
	    }
	  }
	}
      }
    }
  }
  Hsqr.setUpSparse();
  measure(&Hsqr,ioHsqr,excitedStateP.nEigen());
}

//---------------------------------------------------------------------------------------------------//
// These functions are required for computing excited states.
//---------------------------------------------------------------------------------------------------//


int network::gotoNextEigen(){
  //rather self-explanatory
  //Each state is calculated using independent initial states, i.e. the converged ground state is not used as initial guess for the excited state since the projective method used for finding the excited states would map this initial state to zero
  excitedStateP.storeCurrentState(networkState);
  return excitedStateP.loadNextState(networkState);
}

//---------------------------------------------------------------------------------------------------//
// Interface function to compute the expectation value of some operator in MPO representation. 
//---------------------------------------------------------------------------------------------------//

void network::measure(mpo<lapack_complex_double> *const MPOperator, double &lambda, int iEigen){
  mps *measureState;
  if(excitedStateP.nEigen()>=iEigen){
    if(excitedStateP.nEigen()==iEigen){
      measureState=&networkState;
    }
    else{
      excitedStateP.getStoredState(measureState,iEigen);
    }
    globalMeasurement currentMeasurement(MPOperator,measureState);
    currentMeasurement.measureFull(lambda);
  }
}

//---------------------------------------------------------------------------------------------------//

void network::measureLocalOperators(localMpo<lapack_complex_double> *const MPOperator, std::vector<lapack_complex_double> &lambda, int iEigen){
  mps *measureState;
  //This is for measuring the network state during calculation, which is useful for consistency checks
  if(excitedStateP.nEigen()>=iEigen){
    if(excitedStateP.nEigen()==iEigen){
      measureState=&networkState;
    }
    else{
      excitedStateP.getStoredState(measureState,iEigen);
    }
    localMeasurementSeries currentMeasurement(MPOperator,measureState);
    currentMeasurement.measureFull(lambda);
  }
}

//---------------------------------------------------------------------------------------------------//

void network::getEntanglement(std::vector<double> &S, std::vector<std::vector<double> > &spectrum, int iEigen){
  mps *measureState;
  if(excitedStateP.nEigen()==iEigen){
    measureState=&networkState;
  }
  else{
    excitedStateP.getStoredState(measureState,iEigen);
  }
  measureState->getEntanglementEntropy(S,spectrum);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int network::locd(int i){
  return networkDimInfo.locd(i);
}

//---------------------------------------------------------------------------------------------------//

int network::locDMax(int i){
  return networkDimInfo.locDMax(i);
}

//---------------------------------------------------------------------------------------------------//
// Here come further testing functions usable in case of unexplained failure.
//---------------------------------------------------------------------------------------------------//

int network::checkQN(){
  int valid=1;
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    for(int i=0;i<L;++i){
      getLocalDimensions(i);
      for(int si=0;si<ld;++si){
	for(int ai=0;ai<lDR;++ai){
	  for(int aim=0;aim<lDL;++aim){
	    if(real((conservedQNs[iQN].QNLabel(i,ai)-conservedQNs[iQN].QNLabel(i-1,aim)-conservedQNs[iQN].QNLabel(si))) && abs(networkState.global_access(i,si,ai,aim))>0.000001){
	      //if(abs(networkState.global_access(i,si,ai,aim))>0.0001 && (si!=0 || aim!=0)){
	      std::cout<<"Violation of quantum number constraint at "<<"("<<i<<", "<<si<<", "<<ai<<", "<<aim<<"): "<<networkState.global_access(i,si,ai,aim)<<std::endl;
	      std::cout<<"QN Labels: "<<conservedQNs[iQN].QNLabel(i,ai)<<", "<<conservedQNs[iQN].QNLabel(i-1,aim)<<", "<<conservedQNs[iQN].QNLabel(si)<<std::endl;
	      return 1;
	    }
	  }
	}
      }
    }
  }
  if(valid)
    std::cout<<"Validated quantum number constraint\n";
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void network::checkContractions(int i){
  /*
  lapack_complex_double *RTerm;
  int change=0;
  if(i==8){
    getLocalDimensions(i);
    pCtr.Rctr.subContractionStart(RTerm,i);
    if(!store){
      for(int iBlock=0;iBlock<networkState.indexTable().numBlocksLP(i);++iBlock){
	for(int bi=0;bi<lDwR;++bi){
	  for(int j=0;j<networkState.indexTable().rBlockSizeLP(i,iBlock);++j){
	    for(int jp=0;jp<networkState.indexTable().rBlockSizeLP(i,iBlock);++jp){
	      if(abs(RTerm[networkState.indexTable().aiBlockIndexLP(i,iBlock,jp)+bi*lDR+lDR*lDwR*networkState.indexTable().aiBlockIndexLP(i,iBlock,j)]-backupCtr[networkState.indexTable().aiBlockIndexLP(i,iBlock,jp)+bi*lDR+lDR*lDwR*networkState.indexTable().aiBlockIndexLP(i,iBlock,j)])>1e-20)
	      change=1;
	    }
	  }
	}
      }
      if(!change){
	std::cout<<"No update of Rctr was made. Exiting process.\n";
	exit(1);
      }
      else{
	std::cout<<"Successfully updated Rctr\n!";
      }
      delete[] backupCtr;
    }
    else{
      backupCtr=new lapack_complex_double [lDR*lDR*lDwR];
      arraycpy(lDR*lDR*lDwR,RTerm,backupCtr);
    }
  }
  */
  getLocalDimensions(i);
  std::complex<double> *direct;
  std::complex<double> *plambda=new std::complex<double>;
  std::complex<double> *target=new std::complex<double>[ld*lDL*lDR];
  std::complex<double> *currentM=new std::complex<double>[ld*lDL*lDR];
  std::complex<double> *RTerm, *LTerm, *HTerm;
  std::complex<double> lambda;
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  networkH.subMatrixStart(HTerm,i);
  networkState.subMatrixStart(direct,i);
  auxiliary::arraycpy(ld*lDL*lDR,direct,currentM);
  for(int m=0;m<ld*lDL*lDR;++m){
    target[m]=0;
  }
  blockHMatrix BMat(RTerm, LTerm,&networkH,networkDimInfo,Dw,i,&(networkState.indexTable().getLocalIndexTable(i)),&excitedStateP,0,&conservedQNs);
  BMat.prepareInput(currentM);
  BMat.MultMvBlocked(BMat.getCompressedVector(),BMat.getCompressedVector());
  //ARCompStdEig<double, blockHMatrix> eigProblemBlocked(BMat.dim(),1,&BMat,&blockHMatrix::MultMvBlocked,"SR",0,tol,maxIter,BMat.getCompressedVector());
  //eigProblemBlocked.EigenValVectors(BMat.getCompressedVector(),plambda);
  BMat.readOutput(target);
  lambda=*plambda;
  optHMatrix HMat(RTerm,LTerm,&networkH,networkDimInfo,Dw,i,&excitedStateP,0,&conservedQNs);
  //ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,&optHMatrix::MultMvQNConserving,"SR",0,tol,maxIter,currentM);
  //eigProblem.EigenValVectors(currentM,plambda);
  HMat.MultMvQNConserving(currentM,currentM);
  double normb; 
  normb=cblas_dznrm2(ld*lDL*lDR,target,1);
  for(int m=0;m<ld*lDL*lDR;++m){
    target[m]-=currentM[m];
  }
  double test=cblas_dznrm2(ld*lDL*lDR,target,1);
  std::cout<<"Verification: "<<test<<std::endl;
  std::cout<<"Eigenvec norm: "<<cblas_dznrm2(ld*lDL*lDR,currentM,1)<<"\t"<<normb<<std::endl;
  if(test>1e-15){
    //std::cout<<"Eigenvalues: "<<*plambda<<"\t"<<lambda<<std::endl;
    exit(1);
  }
  delete[] target;
  delete[] currentM;
  delete plambda;
}

//---------------------------------------------------------------------------------------------------//
// Function to check whether a state is an equal weight superposition (it only computes the scalar
// product with one product state, since there are too many matching states for a full computation).
//---------------------------------------------------------------------------------------------------//

int network::checkEqualWeightState(){
  std::vector<quantumNumber> dummy;
  mps productState(networkState.getDimInfo(),dummy);
  for(int i=0;i<L;++i){
    getLocalDimensions(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  if(ai==0 && aim==0){
	    if( (i==2 && si==3) || (i==4 && si==4) || (i==3 && si==1) || (si==0 && i!=3 && i!=2)){
	      productState.global_access(i,si,ai,aim)=1;
	    }
	    else{
	      productState.global_access(i,si,ai,aim)=0;
	    }
	  }
	  else{
	    productState.global_access(i,si,ai,aim)=0;
	  }
	}
      }
    }
  }
  overlap test;
  test.loadMPS(&networkState,&productState);
  std::cout<<test.getFullOverlap()<<std::endl;
  exit(1);    
}
