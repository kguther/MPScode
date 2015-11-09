#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <arcomp.h>
#include <arscomp.h>
#include <time.h>
#include <lapacke.h>
#include <cblas.h>
#include "network.h"
#include "arraycreation.h"
#include "arrayprocessing.h"
#include "optHMatrix.h"
#include "siteoptimizer.h"
#include "mpo.h"
#include "mps.h"
#include "globalMeasurement.h"
#include "iterativeMeasurement.h"

//BEWARE: ALL MPS AND MPO NETWORK MATRICES ARE STORED WITH A CONTIGOUS COLUMN INDEX (i.e. transposed with respect to C standard, for better compatibility with LAPACK)

//NAMING CONVENTION: i is always the site index. ALWAYS.

//---------------------------------------------------------------------------------------------------//
// 'Main' file containing the MPS ansatz solver itself (this is not the actual main.cpp, just the
// most important file)
//---------------------------------------------------------------------------------------------------//

//---------------------------------------------------------------------------------------------------//
// Constructors and similar stuff for the basic network class
//---------------------------------------------------------------------------------------------------//

network::network(){
  //Empty networks dont do much, use the generate function to fill them - this allows for separation of declaration and assignment
}

//---------------------------------------------------------------------------------------------------//

network::network(parameters inputpars){
  initialize(inputpars);
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary methods for construction and initialization of network objects
//---------------------------------------------------------------------------------------------------//

void network::initialize(parameters inputpars){
  pars=inputpars;
  d=inputpars.d;
  D=inputpars.D;
  L=inputpars.L;
  Dw=inputpars.Dw;
  devAccuracy=inputpars.acc;
  nSweeps=inputpars.nSweeps;
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access
  networkH.initialize(d,Dw,L);
  //Allocation of MPS - memory of the matrices has to be allocated exactly matching the dimension to use them with lapack and cblas
  networkState.generate(d,D,L);
  alpha=1e-2;
}

//---------------------------------------------------------------------------------------------------//
// These functions can be employed to alter the algorithm parameters N and D during lifetime of a 
// network object. This allows for iteratively increasing of D. They completely take care of all required
// updated within the network, but they should be needed only a few times.
//---------------------------------------------------------------------------------------------------//

void network::setParameterNSweeps(int Nnew){
  nSweeps=Nnew;
}

//---------------------------------------------------------------------------------------------------//

void network::setParameterAlpha(double alphanew){
  alpha=alphanew;
}

//---------------------------------------------------------------------------------------------------//

int network::setParameterD(int Dnew){
  //Decreasing D is neither required nor reasonable in any context
  if(Dnew<D){
    return -1;
  }
  networkState.setParameterD(Dnew);
  //Adapt D
  D=Dnew;
  pars.D=Dnew;
  //Adapt L and R expression
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

int network::solve(double &lambda){  //IMPORTANT TODO: ENHANCE STARTING POINT -> HUGE SPEEDUP
  int errRet;
  int maxIter=400;
  double convergenceQuality;
  clock_t curtime;
  curtime=clock();
  pCtr.initialize(&networkH,&networkState);
  for(int i=L-1;i>0;--i){
    networkState.rightNormalizeState(i);
  }
  networkState.normalizeFinal(1);
  //Allocate a 4D array of dimension LxDxDwxD which is contigous in the last index for the partial contractions Lctr and Rctr
  pCtr.Lctr.global_access(0,0,0,0)=1;
  //In preparation of the first sweep, generate full contraction to the right (first sweeps starts at site 0)
  pCtr.calcCtrFull(1);
  for(int iSweep=0;iSweep<nSweeps;++iSweep){
    std::cout<<"Starting rightsweep\n";
    for(int i=0;i<(L-1);++i){
      //Step of leftsweep
      std::cout<<"Optimizing site matrix\n";
      curtime=clock();
      errRet=optimize(i,maxIter,lambda); 
      curtime=clock()-curtime;
      std::cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      //networkState.leftNormalizeState(i);
      leftEnrichment(i);
      pCtr.calcCtrIterLeft(i+1);
    }
    networkState.normalizeFinal(0);
    std::cout<<"Starting leftsweep\n";
    for(int i=L-1;i>0;--i){
      //Step of rightsweep
      std::cout<<"Optimizing site matrix\n";      
      curtime=clock();
      errRet=optimize(i,maxIter,lambda);
      curtime=clock()-curtime;
      std::cout<<"Optimization took "<<curtime<<" clicks ("<<(float)curtime/CLOCKS_PER_SEC<<" seconds)\n";
      //networkState.rightNormalizeState(i);
      rightEnrichment(i);
      pCtr.calcCtrIterRight(i-1);
    }
    networkState.normalizeFinal(1);
    pCtr.calcCtrIterRightBase(-1,&expectationValue);
    convergenceQuality=convergenceCheck();
    alpha*=.1;
  }
  std::cout<<"Quality of convergence: "<<convergenceQuality<<"\tRequired accuracy: "<<devAccuracy<<std::endl;
  if(convergenceQuality<devAccuracy){
    return 0;
  }
  else{
    return 1;
  }
}


//---------------------------------------------------------------------------------------------------//
// This function optimizes the matrices of a site, where the optimization is mapped to a large sparse
// eigenvalue problem. 
//---------------------------------------------------------------------------------------------------//

int network::optimize(int const i, int const maxIter, double &iolambda){
  //Invokes ARPACK++ to solve the eigenvalue problem
  arcomplex<double> lambda;
  arcomplex<double> *plambda;
  arcomplex<double> *currentM;
  arcomplex<double> *RTerm, *LTerm, *HTerm;
  //Get the current partial contractions and site matrix of the Hamiltonian
  pCtr.Lctr.subContractionStart(LTerm,i);
  pCtr.Rctr.subContractionStart(RTerm,i);
  networkH.subMatrixStart(HTerm,i);
  //Generate matrix which is to be passed to ARPACK++
  optHMatrix HMat(RTerm,LTerm,HTerm,pars,i);
  plambda=&lambda;
  //Using the current site matrix as a starting point allows for much faster convergence as it has already been optimized in previous sweeps (except for the first sweep, this is where a good starting point has to be guessed
  networkState.subMatrixStart(currentM,i);
  ARCompStdEig<double, optHMatrix> eigProblem(HMat.dim(),1,&HMat,&optHMatrix::MultMv,"SR",0,1e-4,maxIter,currentM);
  //One should avoid to hit the maximum number of iterations since this can lead into a suboptimal site matrix, increasing the current energy (although usually not by a lot)
  //So far it seems that the eigensolver either converges quite fast or not at all (i.e. very slow, such that the maximum number of iterations is hit) depending strongly on the tolerance
  int nconv;
  nconv=eigProblem.EigenValVectors(currentM,plambda);
  if(nconv!=1){
    std::cout<<"Failed to converge in iterative eigensolver\n";
    return 1;
  }
  else{
    iolambda=real(lambda);
    std::cout<<setprecision(21)<<"Current energy: "<<real(lambda)<<std::endl;
    return 0;
  }
}

//---------------------------------------------------------------------------------------------------//
// The implementation of the enrichment is very ugly and might still contain severe bugs which
// are irrelevant to the current testing system (Heisenberg chain). It does however, improve convergence
// speed and avoids local minima (at least thats what I expect, there are no problems with local 
// minima till now).
//---------------------------------------------------------------------------------------------------//

void network::leftEnrichment(int const i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Bnew;
  lapack_complex_double *pExpression;
  int lDRR, ldp;
  getLocalDimensions(i);
  lDRR=networkState.locDimR(i+1);
  ldp=locd(i+1);
  int MNumCols=lDR*(1+lDwL);
  int MNumRows=ld*lDL;
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  //Allocate the memory needed for the output of ZUNGBR which is more than the original matrix - also prevents possible segfault when copying (size of target is given for the same reason)
  Mnew=new lapack_complex_double[maxDim*maxDim];
  Bnew=new lapack_complex_double[ldp*lDRR*MNumCols];
  pExpression=new lapack_complex_double[MNumRows*lDR*lDwR];
  getPExpressionLeft(i,pExpression);
  for(int si=0;si<ld;++si){
    //Copy A and B while rearranging to allow for expansion
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+si*lDL+ai*MNumRows]=networkState.global_access(i,si,ai,aim);
      }
    }
    for(int air=0;air<lDRR;++air){
      for(int ai=0;ai<lDR;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=networkState.global_access(i+1,si,air,ai);
      }
    }
    //Add zeros and P-Expression to Mnew and Bnew
    for(int ai=lDR;ai<MNumCols;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+si*lDL+ai*MNumRows]=alpha*pExpression[aim+si*lDL+(ai-lDR)*MNumRows];
      }
    }
    for(int air=0;air<lDRR;++air){
      for(int ai=lDR;ai<MNumCols;++ai){
	Bnew[ai+air*MNumCols+si*lDRR*MNumCols]=0;
      }
    }
  }
  delete[] pExpression;
  //Singular Value Decomposition of Mnew=U*S*V
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  double *diags=new double[containerDim];
  lapack_complex_double *Mnewcpy=new lapack_complex_double[maxDim*maxDim];
  lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int mi=0;mi<MNumCols;++mi){
    for(int ai=0;ai<lDR;++ai){
      Mnewcpy[ai+lDR*mi]*=diags[ai];
    }
  }
  delete[] diags;
  //From here, Mnewcpy is to be treated as a lDR x MNumCols matrix
  //Postprocessing: A=U, B=S*V*Bnew
  lapack_complex_double *BStart;
  lapack_complex_double *networkB;
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  //Multiply S*V into the expanded B (Bnew) to create the normal-sized B (Bstart, direct access to the networkState mps)
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(networkB,i+1,si);
    BStart=Bnew+si*lDRR*MNumCols;
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDR,lDRR,MNumCols,&zone,Mnewcpy,lDR,BStart,MNumCols,&zzero,networkB,lDR);
  }
  delete[] Mnewcpy;
  delete[] Bnew;
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=Mnew[aim+si*lDL+ai*MNumRows];
      }
    }
  }
  delete[] Mnew;
}

//---------------------------------------------------------------------------------------------------//

void network::rightEnrichment(int const i){
  lapack_complex_double *Mnew;
  lapack_complex_double *Anew;
  lapack_complex_double *pExpression;
  int lDLL, ldm;
  getLocalDimensions(i);
  lDLL=networkState.locDimL(i-1);
  ldm=locd(i-1);
  int MNumCols=ld*lDR;
  int MNumRows=lDL*(1+lDwL);
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  Mnew=new lapack_complex_double[maxDim*maxDim];
  Anew=new lapack_complex_double[ldm*lDLL*MNumRows];
  pExpression=new lapack_complex_double[ld*lDR*lDL*lDwL];
  getPExpressionRight(i,pExpression);
  //Copy A and B while rearranging to allow for expansion
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	Mnew[aim+ai*MNumRows+si*lDR*MNumRows]=networkState.global_access(i,si,ai,aim);
      }
    }
    for(int aim=0;aim<lDL;++aim){
      for(int aimm=0;aimm<lDLL;++aimm){
	Anew[aimm+aim*ld*lDLL+si*lDLL]=networkState.global_access(i-1,si,aim,aimm);
      }
    }
    //Add zeros and P-Expression to Mnew and Bnew
    for(int ai=0;ai<lDR;++ai){
      for(int aim=lDL;aim<MNumRows;++aim){
	Mnew[aim+ai*MNumRows+si*lDR*MNumRows]=alpha*pExpression[aim-lDL+si*lDR*lDwL*lDL+ai*lDwL*lDL];
      }
    }
    for(int aim=lDL;aim<MNumRows;++aim){
      for(int aimm=0;aimm<lDLL;++aimm){
        Anew[aimm+aim*ld*lDLL+si*lDLL]=0;
      }
    }
  }
  delete[] pExpression;
  //Singular Value Decomposition of Mnew=U*S*V
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  double *diags=new double[containerDim];
  lapack_complex_double *Mnewcpy=new lapack_complex_double[maxDim*maxDim];
  lapackSVD(MNumCols,MNumRows,Mnew,Mnewcpy,diags);
  //Mnewcpy -> A, S*Mnew->Multiply to B
  //Postprocessing: Truncate S to lDR eigenvalues, U to dimension ld*lDL x lDR (from ld*lDL x ld*lDL) if neccessary
  //Truncation is carried out implicitly by only adressing part of the matrix
  for(int aim=0;aim<lDL;++aim){
    for(int mi=0;mi<MNumRows;++mi){
      Mnew[mi+aim*MNumRows]*=diags[aim];
    }
  }
  delete[] diags;
  //From here, Mnew is to be treated as a MNumRows x lDL matrix
  //Postprocessing: A=U, B=S*V*Bnew
  lapack_complex_double *AStart=new lapack_complex_double[lDLL*MNumRows];
  lapack_complex_double *networkA;
  lapack_complex_double zone=1.0;
  lapack_complex_double zzero=0.0;
  //Multiply S*V into the expanded A (AStart) to create the normal-sized A (networkA, direct access to the networkState mps)
  for(int si=0;si<ld;++si){
    networkState.subMatrixStart(networkA,i-1,si);
    for(int mi=0;mi<MNumRows;++mi){
      for(int aimm=0;aimm<lDLL;++aimm){
	AStart[aimm+mi*lDLL]=Anew[aimm+si*lDLL+mi*ld*lDLL];
      }
    }
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lDLL,lDL,MNumRows,&zone,AStart,lDLL,Mnew,MNumRows,&zzero,networkA,lDLL);
  }
  delete[] AStart;
  delete[] Mnew;
  delete[] Anew;
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int aim=0;aim<lDL;++aim){
	networkState.global_access(i,si,ai,aim)=Mnewcpy[aim+si*lDR*MNumCols+ai*MNumCols];
      }
    }
  }
  delete[] Mnewcpy;
}

//---------------------------------------------------------------------------------------------------//
// Auxiliary functions to determine the heuristic expansion term for the enrichment step.
//---------------------------------------------------------------------------------------------------//

void network::getPExpressionLeft(int const i, lapack_complex_double *pExpr){
  lapack_complex_double simpleContainer;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> innerContainer(ld,lDL,lDR,lDwL);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDL;++aim){
      for(int ai=0;ai<lDR;++ai){
	for(int bim=0;bim<lDwL;++bim){
	  simpleContainer=0;
	  for(int aimp=0;aimp<lDL;++aimp){
	    simpleContainer+=pCtr.Lctr.global_access(i,aim,bim,aimp)*networkState.global_access(i,si,ai,aimp);
	  }
	  innerContainer.global_access(si,aim,ai,bim)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int bi=0;bi<lDwR;++bi){
	for(int aim=0;aim<lDL;++aim){
	  //BEWARE: IM UNSURE WHETHER bi OR ai IS THE OUTERMOST INDEX - FORGET IT, ITS bi
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bim=0;bim<lDwL;++bim){
	      simpleContainer+=networkH.global_access(i,si,sip,bi,bim)*innerContainer.global_access(sip,aim,ai,bim);
	    }
	  }
	  pExpr[aim+lDL*si+bi*lDL*d*lDR+ai*lDL*d]=simpleContainer;    
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void network::getPExpressionRight(int const i, lapack_complex_double *pExpr){
  lapack_complex_double simpleContainer;
  getLocalDimensions(i);
  tmpContainer<lapack_complex_double> innerContainer(ld,lDL,lDR,lDwR);
  for(int si=0;si<ld;++si){
    for(int aim=0;aim<lDL;++aim){
      for(int ai=0;ai<lDR;++ai){
	for(int bi=0;bi<lDwR;++bi){
	  simpleContainer=0;
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer+=pCtr.Rctr.global_access(i,ai,bi,aip)*networkState.global_access(i,si,aip,aim);
	  }
	  innerContainer.global_access(si,aim,ai,bi)=simpleContainer;
	}
      }
    }
  }
  for(int si=0;si<ld;++si){
    for(int ai=0;ai<lDR;++ai){
      for(int bim=0;bim<lDwL;++bim){
	for(int aim=0;aim<lDL;++aim){
	  simpleContainer=0;
	  for(int sip=0;sip<ld;++sip){
	    for(int bi=0;bi<lDwR;++bi){
	      simpleContainer+=networkH.global_access(i,si,sip,bi,bim)*innerContainer.global_access(sip,aim,ai,bi);
	    }
	  }
	  pExpr[aim*lDwL+bim+ai*lDL*lDwL+si*lDL*lDwL*lDR]=simpleContainer;
	}
      }
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


/*void network::calcHSqrExpectationValue(double &ioHsqr){
  dynamicContainer<lapack_complex_double> currentP;
  dynamic5DContainer<lapack_complex_double> innerContainer;
  dynamic5DContainer<lapack_complex_double> middleContainer;
  dynamic5DContainer<lapack_complex_double> outerContainer;
  lapack_complex_double simpleContainer;
  currentP.generate(1,1,1,1);
  currentP.global_access(0,0,0,0)=1;
  std::cout<<"Starting evaluation of convergence quality\n";
  for(int i=0;i<L;++i){
    getLocalDimensions(i);
    innerContainer.generate(ld,lDwL,lDwL,lDR,lDL);
    for(int sipp=0;sipp<ld;++sipp){
      for(int bimp=0;bimp<lDwL;++bimp){
	for(int bim=0;bim<lDwL;++bim){
	  for(int aip=0;aip<lDR;++aip){
	    for(int aim=0;aim<lDL;++aim){
	      simpleContainer=0;
	      for(int aimp=0;aimp<lDL;++aimp){
		simpleContainer+=networkState.global_access(i,sipp,aip,aimp)*currentP.global_access(bim,bimp,aim,aimp);
	      }
	      innerContainer.global_access(sipp,bim,bimp,aip,aim)=simpleContainer;
	    }
	  }
	}
      }
    }
    middleContainer.generate(ld,lDwR,lDwL,lDR,lDL);
    for(int sip=0;sip<ld;++sip){
      for(int bim=0;bim<lDwL;++bim){
	for(int bip=0;bip<lDwR;++bip){
	  for(int aip=0;aip<lDR;++aip){
	    for(int aim=0;aim<lDL;++aim){
	      simpleContainer=0;
	      for(int sipp=0;sipp<ld;++sipp){
		for(int bimp=0;bimp<lDwL;++bimp){
		  simpleContainer+=innerContainer.global_access(sipp,bim,bimp,aip,aim)*networkH.global_access(i,sip,sipp,bip,bimp);
		}
	      }
	      middleContainer.global_access(sip,bip,bim,aip,aim)=simpleContainer;
	    }
	  }
	}
      }
    }
    outerContainer.generate(ld,lDwR,lDwR,lDR,lDL);
    for(int si=0;si<ld;++si){
      for(int bi=0;bi<lDwR;++bi){
	for(int bip=0;bip<lDwR;++bip){
	  for(int aip=0;aip<lDR;++aip){
	    for(int aim=0;aim<lDL;++aim){
	      simpleContainer=0;
	      for(int sip=0;sip<ld;++sip){
		for(int bim=0;bim<lDwL;++bim){
		  simpleContainer+=middleContainer.global_access(sip,bip,bim,aip,aim)*networkH.global_access(i,si,sip,bi,bim);
		}
	      }
	      outerContainer.global_access(si,bi,bip,aip,aim)=simpleContainer;
	    }
	  }
	}
      }
    }
    currentP.generate(lDwR,lDwR,lDR,lDR);
    for(int bi=0;bi<lDwR;++bi){
      for(int bip=0;bip<lDwR;++bip){
	for(int ai=0;ai<lDR;++ai){
	  for(int aip=0;aip<lDR;++aip){
	    simpleContainer=0;
	    for(int si=0;si<ld;++si){
	      for(int aim=0;aim<lDL;++aim){
		simpleContainer+=outerContainer.global_access(si,bi,bip,aip,aim)*conj(networkState.global_access(i,si,ai,aim));
	      }
	    }
	    currentP.global_access(bi,bip,ai,aip)=simpleContainer;
	  }
	}
      }
    }
  }
  ioHsqr=real(currentP.global_access(0,0,0,0));
  }*/
