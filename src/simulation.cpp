#include <fstream>
#include <sstream>
#include "simulation.h"
#include "projector.h"
#include "globalMeasurement.h"
#include "localMeasurementSeries.h"

simulation::simulation(problemParameters &parsIn, simulationParameters &simParsIn, double const J, double const g, int pathPoints, std::string &targetFile):
  simPars(simParsIn),
  pars(parsIn),
  pathLength(pathPoints),
  filePrefix(targetFile),
  parDirection(std::complex<double>(J,g))
{
  // simulation class can only use nStages=1 (by choice of algorithm)
  TensorNetwork.initialize(pars,simPars);
  E0.resize(pars.nEigs);
  convergedEigens.resize(pars.nEigs);
  double matEls;
  int L=pars.L;
  measureTask.clear();
  localMeasureTask.clear();
  parDirection*=1.0/(abs(parDirection)*100);
  particleNumber.initialize(pars.d.maxd(),2,pars.L);
  subChainParity.initialize(pars.d.maxd(),1,pars.L);
  for(int i=0;i<pars.L;++i){
    for(int bi=0;bi<1;++bi){
      for(int bim=0;bim<1;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip)*(delta(si,1)+delta(si,0)-delta(si,2)-delta(si,3));
	    subChainParity.global_access(i,si,sip,bi,bim)=matEls;
	  }
	}
      }
    }
  }
  for(int i=0;i<pars.L;++i){
    for(int bi=0;bi<2;++bi){
      for(int bim=0;bim<2;++bim){
	for(int si=0;si<pars.d.maxd();++si){
	  for(int sip=0;sip<pars.d.maxd();++sip){
	    matEls=delta(si,sip);
	    if(i!=0 && i!=L-1 && bi==1 && bim==0){
	      matEls=0;
	    }
	    if(bi==0 && bim==particleNumber.locDimL(i)-1){
	      matEls*=(delta(si,1)+delta(si,2)+2*delta(si,3));
	    }
	    particleNumber.global_access(i,si,sip,bi,bim)=matEls;
	  }
	}
      }
    }
  }
  double spinQN;
  TensorNetwork.check=&particleNumber;
  TensorNetwork.checkParity=&subChainParity;
}

//---------------------------------------------------------------------------------------------------//

void simulation::setMeasurement(mpo<std::complex<double> > &MPOperator, std::string &operatorName){
  measureTask.push_back(MPOperator);
  operatorNames.push_back(operatorName);
}

//---------------------------------------------------------------------------------------------------//

void simulation::setLocalMeasurement(localMpo<std::complex<double> > &localMPOperator, std::string &localOperatorName){
  localMeasureTask.push_back(localMPOperator);
  localOperatorNames.push_back(localOperatorName);
}

//---------------------------------------------------------------------------------------------------//

void simulation::run(){
  singleRun();
  for(int nRun=1;nRun<pathLength;++nRun){
    parDirection*=1.0/(abs(parDirection)*100)*nRun;
    singleRun();
  }
}

//---------------------------------------------------------------------------------------------------//

void simulation::singleRun(){
  int hInfo;
  double J,g;
  projector *stateRep;
  std::vector<quantumNumber> *QNs;
  mps *measureState;
  std::vector<double> expectationValues;
  std::vector<std::vector<std::complex<double> > > localExpectationValues;
  TensorNetwork.getProjector(stateRep);
  J=1+parDirection.real();
  g=1+parDirection.imag();
  hInfo=writeHamiltonian(TensorNetwork,J,g);
  if(hInfo){
    std::cout<<"Invalid bond dimension for the construction of H. Terminating process.\n";
    exit(1);
  }
  TensorNetwork.solve(&E0[0]);
  expectationValues.resize(measureTask.size());
  localExpectationValues.resize(localMeasureTask.size());
  std::ofstream ofs;
  std::ostringstream numStrJ, numStrg;
  numStrJ<<J;
  numStrg<<g;
  std::string fileName;
  fileName+=filePrefix;
  fileName+="_J_";
  fileName+=numStrJ.str();
  fileName+="_g_";
  fileName+=numStrg.str();
  fileName+=".txt";
  ofs.open(fileName.c_str());
  for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
    if(pars.nEigs==1){
      measureState=0;
    }
    else{
      stateRep->getStoredState(measureState,iEigen);
    }
    for(int iM=0;iM<measureTask.size();++iM){
      measure(&measureTask[iM],expectationValues[iM],measureState);
      ofs<<operatorNames[iM]<<"\t";
    }
    ofs<<std::endl;
    for(int iM=0;iM<measureTask.size();++iM){
      ofs<<expectationValues[iM]<<"\t";
    }
    ofs<<std::endl;
    for(int iM=0;iM<localMeasureTask.size();++iM){
      measureLocal(&localMeasureTask[iM],localExpectationValues[iM],measureState);
      ofs<<localOperatorNames[iM]<<"\t";
    }
    ofs<<std::endl;
    if(localExpectationValues.size()>0){
      for(int i=0;i<localExpectationValues[0].size();++i){
	for(int iM=0;iM<localMeasureTask.size();++iM){
	  ofs<<localExpectationValues[iM][i]<<"\t";
	}
	ofs<<std::endl;
      }
    }
  }
  TensorNetwork.resetConvergence();
}

//---------------------------------------------------------------------------------------------------//

int simulation::measure(mpo<lapack_complex_double> *MPOperator, double &expectationValue, mps *MPState){
  if(MPState){
    globalMeasurement currentMeasurement(MPOperator,MPState);
    currentMeasurement.measureFull(expectationValue);
    return 0;
  }
  return TensorNetwork.measure(MPOperator,expectationValue);
}

//---------------------------------------------------------------------------------------------------//

int simulation::measureLocal(localMpo<lapack_complex_double> *localMPOperator, std::vector<std::complex<double> > &result, mps *MPState){
  if(MPState){
    localMeasurementSeries currentMeasurement(localMPOperator,MPState);
    currentMeasurement.measureFull(result);
    return 0;
  }
  return TensorNetwork.measureLocalOperators(localMPOperator,result);
}
