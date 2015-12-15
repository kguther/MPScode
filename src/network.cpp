#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <arcomp.h>
#include <arscomp.h>
#include <time.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "mpo.h"
#include "mps.h"
#include "globalMeasurement.h"
#include "iterativeMeasurement.h"
#include "projector.h"

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//NAMING CONVENTION: i is always the site index. ALWAYS.

//TODO: MAKE USE OF RETURN VALUES WHERE INTENDED.

//---------------------------------------------------------------------------------------------------//
// 'Main' file containing the MPS ansatz solver itself (this is not the actual main.cpp, just the
// most important file)
//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
// Constructors and similar stuff for the basic network class
//---------------------------------------------------------------------------------------------------//

network::network():
  nConverged(0)
  //Empty networks dont do much, use the generate function to fill them - this allows for separation of declaration and assignment
{}

//---------------------------------------------------------------------------------------------------//

network::network(problemParameters inputpars, simulationParameters inputsimPars){
  initialize(inputpars,inputsimPars);
}

//---------------------------------------------------------------------------------------------------//

network::~network(){
  delete[] nConverged;
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary methods for construction and initialization of network objects
//---------------------------------------------------------------------------------------------------//

void network::initialize(problemParameters inputpars, simulationParameters inputsimPars){
  pars=inputpars;
  d=inputpars.d;
  D=inputsimPars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  simPars=inputsimPars;
  networkDimInfo.initialize(d,D,L);
  nConverged=new int[pars.nEigs];
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    nConverged[iEigen]=1;
  }
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access
  networkH.initialize(d,Dw,L);
  //Allocation of MPS - memory of the matrices has to be allocated exactly matching the dimension to use them with lapack and cblas
  //All states have to be stored, for reusability in later solve() with other simulation parameters
  //Somewhat unelegant way to handle loading of the stored states in the first solve()
  conservedQNs.resize(pars.nQNs);
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    conservedQNs[iQN].initialize(networkDimInfo,pars.QNconserved[iQN],pars.QNLocalList+iQN*d);
  }
  networkState.generate(d,D,L,&conservedQNs);
  excitedStateP.initialize(pars.nEigs);
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    excitedStateP.storeOrthoState(networkState,iEigen);
  }
}

//---------------------------------------------------------------------------------------------------//

void network::loadNetworkState(mps &source){
  networkState.mpsCpy(source);
}

//---------------------------------------------------------------------------------------------------//
// These functions can be employed to alter the algorithm parameters N and D during lifetime of a 
// network object. This allows for iteratively increasing D. They completely take care of all required
// updated within the network, but they should be needed only a few times.
//---------------------------------------------------------------------------------------------------//

int network::setSimParameters(simulationParameters newPars){
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
  if(Dnew<D){
    return -1;
  }
  networkState.setParameterD(Dnew);
  //All stored states have to be brought into the correct form for compatibility with current D
  excitedStateP.setParameterD(Dnew);
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    conservedQNs[iQN].setParameterD(Dnew);
  }
  //Adapt D
  D=Dnew;
  simPars.D=Dnew;
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void network::getLocalDimensions(int const i){
  lDL=networkState.locDimL(i);
  lDR=networkState.locDimR(i);
  ld=locd(i);
  lDwR=networkH.locDimR(i);
  lDwL=networkH.locDimL(i);
}

//---------------------------------------------------------------------------------------------------//
// This is the main solver, using the networkState as initial and final state and networkH as MPO
// representation of the Hamiltonian. Computed are eigenstate (initial state is overwritten) and 
// corresponding eigenvalue (output).
//---------------------------------------------------------------------------------------------------//

int network::solve(double *lambda){  //IMPORTANT TODO: ENHANCE STARTING POINT -> HUGE SPEEDUP
  int maxIter=5000;
  int stepRet;
  int cshift=-100;
  double convergenceQuality;
  double alpha;
  double tol;
  double spinCheck;
  overlap test;
  test.loadMPS(&networkState,&networkState);
  alpha=simPars.alpha;
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    pCtr.initialize(&networkH,&networkState);
    for(int i=L-1;i>0;--i){
      normalize(i,0,0);
    }
    networkState.normalizeFinal(1);
    std::cout<<"Finished right normalization\n";    
    std::cout<<"Norm of state: "<<test.getFullOverlap()<<std::endl;
    pCtr.Lctr.global_access(0,0,0,0)=1;
    //In preparation of the first sweep, generate full contraction to the right (first sweeps starts at site 0)
    pCtr.calcCtrFull(1);
    if(iEigen){
      shift=-2*abs(lambda[0]);
    }
    else{
      shift=0.0;
    }
    if(pars.nQNs){
      shift+=cshift;
    }
    std::cout<<"Computing state "<<iEigen<<std::endl;
    alpha=simPars.alpha;
    tol=simPars.tolInitial;
    //load all pairings with the current state and previous ones into the scalar products
    excitedStateP.loadScalarProducts(&networkState,iEigen);
    for(int iSweep=0;iSweep<simPars.nSweeps;++iSweep){
      if(!nConverged[iEigen]){
	break;
      }
      //actual sweep is executed here
      sweep(maxIter,tol,alpha,lambda[iEigen]);
      for(int i=L-1;i>0;--i){
      	normalize(i,0,0);
      }
      networkState.normalizeFinal(1);
      //In calcCtrIterRightBase, the second argument has to be a pointer, because it usually is an array. No call-by-reference here.
      pCtr.calcCtrIterRightBase(-1,&expectationValue);
      //Compute full scalar product - although not required, it is nice to know and requires very little computational effort
      convergenceQuality=convergenceCheck();
      if(convergenceQuality<simPars.devAccuracy){
	nConverged[iEigen]=0;
      }
      alpha*=.1;
      if(tol>simPars.tolMin){
	tol*=pow(simPars.tolMin/simPars.tolInitial,1.0/simPars.nSweeps);
      }
      std::cout<<"Quality of convergence: "<<convergenceQuality<<"\tRequired accuracy: "<<simPars.devAccuracy<<std::endl;
      excitedStateP.updateScalarProducts(0,-1);
      for(int prev=0;prev<iEigen;++prev){
	std::cout<<"Overlap with state "<<prev<<" is: "<<excitedStateP.fullOverlap(prev)<<std::endl;
      }
    }
    stepRet=gotoNextEigen();
    if(!stepRet){
      std::cout<<"LOADED STATES. PREPARED COMPUTATION OF NEXT EIGENSTATE"<<std::endl;
    }
  }
  if(pars.nQNs){
    lambda[0]-=cshift;
  }
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

void network::sweep(double const maxIter, double const tol, double const alpha, double &lambda){
  clock_t curtime;
  int errRet;
  lapack_complex_double stateNorm;
  std::cout<<"Starting rightsweep\n";
  for(int i=0;i<(L-1);++i){
    //Step of leftsweep
    std::cout<<"Optimizing site matrix"<<std::endl;
    curtime=clock();
    errRet=optimize(i,maxIter,tol,lambda);
    curtime=clock()-curtime;
    std::cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n\n";
    normalize(i,1,alpha);
    //Here, the scalar products with lower lying states are updated
    excitedStateP.updateScalarProducts(i,1);
    pCtr.calcCtrIterLeft(i+1);
  }
  networkState.normalizeFinal(0);
  std::cout<<"Starting leftsweep\n";
  for(int i=L-1;i>0;--i){
    //Step of rightsweep
    std::cout<<"Optimizing site matrix"<<std::endl;    
    curtime=clock();
    errRet=optimize(i,maxIter,tol,lambda);
    curtime=clock()-curtime;
    std::cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n\n";
    normalize(i,0,alpha);
    //same as above for the scalar products with lower lying states
    excitedStateP.updateScalarProducts(i,-1);
    pCtr.calcCtrIterRight(i-1);
  }
  networkState.normalizeFinal(1);
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int const i, int const maxIter, double const tol, double &iolambda){
  //Invokes ARPACK++ to solve the eigenvalue problem
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  arcomplex<double> *currentM;
  arcomplex<double> *RTerm, *LTerm, *HTerm;
  double spinCheck;
  void (optHMatrix::*multMV)(arcomplex<double> *v, arcomplex<double> *w);
  //Get the projector onto the space orthogonal to any lower lying states
  //Get the current partial contractions and site matrix of the Hamiltonian
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  networkH.subMatrixStart(HTerm,i);
  excitedStateP.getProjector(i);
  if(pars.nQNs){
    multMV=&optHMatrix::MultMvQNConserving;
  }
  else{
    multMV=&optHMatrix::MultMv;
  }
  //Generate matrix which is to be passed to ARPACK++
  optHMatrix HMat(RTerm,LTerm,HTerm,networkDimInfo,Dw,i,&excitedStateP,shift,&conservedQNs);
  plambda=&lambda;
  //Using the current site matrix as a starting point allows for much faster convergence as it has already been optimized in previous sweeps (except for the first sweep, this is where a good starting point has to be guessed
  networkState.subMatrixStart(currentM,i);
  measure(check,spinCheck);
  std::cout<<"Spin before optimizing: "<<spinCheck<<std::endl;
  //Note that the types given do and have to match the ones in the projector class if more than one eigenvalue is computed
  ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,multMV,"SR",0,tol,maxIter,currentM);
  //One should avoid to hit the maximum number of iterations since this can lead into a suboptimal site matrix, increasing the current energy (although usually not by a lot)
  //So far it seems that the eigensolver either converges quite fast or not at all (i.e. very slow, such that the maximum number of iterations is hit) depending strongly on the tolerance
  int nconv;
  nconv=eigProblem.EigenValVectors(currentM,plambda);
  measure(check,spinCheck);
  std::cout<<"Spin after optimizing: "<<spinCheck<<std::endl;
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver, number of Iterations taken: "<<maxIter<<" With tolerance "<<tol<<std::endl;
    return 1;
  }
  else{
    iolambda=real(lambda);
    std::cout<<setprecision(21)<<"Current energy: "<<real(lambda)-shift<<std::endl;
    return 0;
  }
}

//---------------------------------------------------------------------------------------------------//
// Functions for normalization of site matrices, creation of canonical MPS and adressing the subspace
// expansion
//---------------------------------------------------------------------------------------------------//

void network::normalize(int const i, int const direction, double const alpha){
  int enrichment=0;
  if(!pars.nQNs && alpha>1e-20){
    enrichment=1;
  }
  if(direction){
    if(enrichment){
      leftEnrichment(alpha,i);
    }
    else{
      networkState.leftNormalizeState(i);
    }
  }
  else{
    if(enrichment){
      rightEnrichment(alpha,i);
    }
    else{
      networkState.rightNormalizeState(i);
    }
  }
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
  return abs(stdDeviation);
}

//---------------------------------------------------------------------------------------------------//

void network::calcHSqrExpectationValue(double &ioHsqr){
  mpo<lapack_complex_double> Hsqr(d,Dw*Dw,L);
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
  measure(&Hsqr,ioHsqr);
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

int network::measure(mpo<lapack_complex_double> *MPOperator, double &lambda){
  globalMeasurement currentMeasurement(MPOperator,&networkState);
  lambda=currentMeasurement.measureFull();
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int network::locd(int const i){
  return d;
}

//---------------------------------------------------------------------------------------------------//

int network::locDMax(int const i){
  if(i<=L/2){
    return pow(d,i+1);
  }
  return pow(d,L-i);
}

//---------------------------------------------------------------------------------------------------//
// These functions compute the normalization of the network to the left of some site i. This is for
// testing purpose only since the correct algorithm ensures the results to be trivial.
//---------------------------------------------------------------------------------------------------//

void network::leftNormalizationMatrixFull(){
  lapack_complex_double ***psi;
  create3D(L,D,D,&psi);
  for(int i=0;i<L;++i){
    for(int ai=0;ai<D;++ai){
      for(int aip=0;aip<D;++aip){
	psi[i][ai][aip]=0.0/0.0;
      }
    }
  }
  psi[0][0][0]=1;
  std::cout<<"Printing Psi expressions\n";
  for(int i=1;i<L;++i){
    leftNormalizationMatrixIter(i,psi[0][0]);
    matrixprint(D,D,psi[i][0]);
  }
  std::cout<<"That's it\n";
  delete3D(&psi);
}

//---------------------------------------------------------------------------------------------------//

void network::leftNormalizationMatrixIter(int i, lapack_complex_double *psi){
  lapack_complex_double psiContainer;
  for(int aim=0;aim<networkState.locDimL(i);++aim){
    for(int aimp=0;aimp<networkState.locDimL(i);++aimp){
      psiContainer=0;
      for(int si=0;si<d;++si){
	for(int aimm=0;aimm<networkState.locDimL(i-1);aimm++){
	  for(int aimmp=0;aimmp<networkState.locDimL(i-1);aimmp++){
	    psiContainer+=networkState.global_access(i-1,si,aim,aimm)*conj(networkState.global_access(i-1,si,aimp,aimm))*psi[(i-1)*D*D+aimm*D+aimmp];
	  }
	}
      }
      psi[i*D*D+aim*D+aimp]=psiContainer;
    }
  }	
}

//---------------------------------------------------------------------------------------------------//

int network::checkQN(){
  int valid=1;
  for(int iQN=0;iQN<pars.nQNs;++iQN){
    for(int i=0;i<L;++i){
      getLocalDimensions(i);
      for(int si=0;si<ld;++si){
	for(int ai=0;ai<lDR;++ai){
	  for(int aim=0;aim<lDL;++aim){
	    if((conservedQNs[iQN].QNLabel(i,ai)-conservedQNs[iQN].QNLabel(i-1,aim)-conservedQNs[iQN].QNLabel(si)) && abs(networkState.global_access(i,si,ai,aim))>0.000001){
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
