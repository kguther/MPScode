#include <fstream>
#include <iostream>
#include <sstream>
#include "simulation.h"
#include "projector.h"
#include "globalMeasurement.h"
#include "localMeasurementSeries.h"

simulation::simulation(){
}

//---------------------------------------------------------------------------------------------------//

simulation::simulation(problemParameters &parsIn, simulationParameters &simParsIn, double J, double g, double WIn, int pathPoints, double stepSize, std::string const &targetFile):
  simPars(simParsIn),
  pars(parsIn),
  W(WIn),
  pathLength(pathPoints),
  filePrefix(targetFile),
  measureEE(0),
  measureES(0),
  scaling(stepSize),
  parDirection(std::complex<double>(J,g)),
  TensorNetwork(network(pars,simPars)),
  particleNumber(mpo<lapack_complex_double>(pars.d.maxd(),2,pars.L)),
  subChainParity(mpo<lapack_complex_double>(pars.d.maxd(),1,pars.L))
{
 // simulation class can only use nStages=1 (by choice of algorithm)
  E0.resize(pars.nEigs);
  dE.resize(pars.nEigs);
  convergedEigens.resize(pars.nEigs);
  double matEls;
  int L=pars.L;
  if(abs(parDirection)>1e-20){
    parDirection*=1.0/(scaling*abs(parDirection));
  }
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
  TensorNetwork.check=&particleNumber;
  TensorNetwork.checkParity=&subChainParity;
  localMeasureTask.clear();
  measureTask.clear();
}

//---------------------------------------------------------------------------------------------------//

void simulation::setMeasurement(mpo<std::complex<double> > &MPOperator, std::string &operatorName){
  //defines a measurement of some global operator MPOperator in MPO representation
  measureTask.push_back(MPOperator);
  operatorNames.push_back(operatorName);
}

//---------------------------------------------------------------------------------------------------//

void simulation::setLocalMeasurement(localMpo<std::complex<double> > &localMPOperator, std::string &localOperatorName){
  //defines a measurement for some operator depending on the site. The expectation value is computed for all sites right to the initial site. 
  localMeasureTask.push_back(localMPOperator);
  localOperatorNames.push_back(localOperatorName);
}

//---------------------------------------------------------------------------------------------------//

void simulation::setEntanglementMeasurement(){
  measureEE=1;
}

//---------------------------------------------------------------------------------------------------//

void simulation::setEntanglementSpectrumMeasurement(){
  measureEE=1;
  measureES=1;
}

//---------------------------------------------------------------------------------------------------//

void simulation::run(){
  //Solve the system for different parameters (J,g) along a straight line in radial direction 
  for(int nRun=1;nRun<pathLength+1;++nRun){
    if(abs(parDirection)>1e-20){
      parDirection*=1.0/(scaling*abs(parDirection))*nRun;
    }
    singleRun();
  }
}

//---------------------------------------------------------------------------------------------------//

void simulation::singleRun(){
  //Get nEigs eigenstates for a fixed set of parameters J,g of the hamiltonian
  int hInfo;
  double J,g;
  //Containers for measurements
  std::vector<double> expectationValues;
  std::vector<std::vector<std::complex<double> > > localExpectationValues;
  J=1+parDirection.real();
  g=1+parDirection.imag();
  std::cout<<J<<" "<<g<<std::endl;
  hInfo=writeHamiltonian(TensorNetwork,J,g,W);
  if(hInfo){
    std::cout<<"Invalid bond dimension for the construction of H. Terminating process.\n";
    exit(1);
  }
  //Actual DMRG
  TensorNetwork.solve(E0,dE);
  expectationValues.resize(measureTask.size());
  localExpectationValues.resize(localMeasureTask.size());
  //Measure previously set operators and write results into the result file
  if(localMeasureTask.size()>0 || measureTask.size()>0 || measureEE){
    std::cout<<"Measuring correlation functions\n";
    std::ofstream ofs;
    std::string finalName, fileName;
    std::ostringstream compositeName;
    compositeName<<filePrefix<<"_J_"<<J<<"_g_"<<g;
    finalName=compositeName.str();
    for(int m=0;m<finalName.length();++m){
      if(finalName[m]=='.'){
	finalName.erase(m,1);
      }
    }
    for(int iEigen=0;iEigen<pars.nEigs;++iEigen){
      if(iEigen>0){
	compositeName.str("");
	compositeName<<(iEigen+1);
	fileName=finalName+"_state_"+compositeName.str()+".txt";
      }
      else{
	fileName=finalName+".txt";
      }
      ofs.open(fileName.c_str());
      ofs<<"Values for state number "<<iEigen<<" with energy "<<E0[iEigen]<<" and energy variance "<<dE[iEigen]<<std::endl;
      //The problem parameters are written into the first lines
      ofs<<"L\tN\tsubchain parity\tJ\tg\tW\tE\tvariance of energy\n";
      ofs<<pars.L<<"\t"<<real(pars.QNconserved[0])<<"\t"<<imag(pars.QNconserved[0])<<"\t"<<J<<"\t"<<g<<"\t"<<W<<"\t"<<E0[iEigen]<<"\t"<<dE[iEigen]<<std::endl;
      //First, global measurements are performed (this is used rarely)
      for(int iM=0;iM<measureTask.size();++iM){
	TensorNetwork.measure(&measureTask[iM],expectationValues[iM],iEigen);
	ofs<<operatorNames[iM]<<"\t";
      }
      //Then, all results are written into the result file
      if(measureTask.size())
	ofs<<std::endl;
      for(int iM=0;iM<measureTask.size();++iM){
	ofs<<expectationValues[iM]<<"\t";
      }
      if(measureTask.size())
	ofs<<std::endl;
      //Now, the same is done for the local measurements (this is what is usually interesting)
      for(int iM=0;iM<localMeasureTask.size();++iM){
	TensorNetwork.measureLocalOperators(&localMeasureTask[iM],localExpectationValues[iM],iEigen);
	ofs<<localOperatorNames[iM]<<"\t";
      }
      if(measureEE){
	ofs<<"Entanglement Entropy"<<"\t";
      }
      //Knowing the initial site is useful to distinguish bulk and edge functions
      ofs<<std::endl;
      for(int iM=0;iM<localMeasureTask.size();++iM){
	ofs<<localMeasureTask[iM].currentSite()<<"\t";
      }
      if(measureEE){
	ofs<<0<<std::endl;
      }
      ofs<<std::endl;
      std::vector<double> S;
      std::vector<std::vector<double> > spec;
      if(measureEE){
	TensorNetwork.getEntanglement(S,spec,iEigen);
      }
      if(localExpectationValues.size()>0 || measureEE){
	for(int i=0;i<pars.L;++i){
	  for(int iM=0;iM<localMeasureTask.size();++iM){
	    if(i<localExpectationValues[iM].size()){
	      ofs<<abs(localExpectationValues[iM][i])<<"\t";
	    }
	    else{
	      ofs<<"\t";
	    }
	  }
	  if(measureEE && i<S.size()){
	    ofs<<S[i]<<"\t";
	  }
	  ofs<<std::endl;
	}
      }
      if(measureES){
	if(iEigen>0){
	  compositeName.str("");
	  compositeName<<(iEigen+1);
	  fileName=finalName+"_state_"+compositeName.str()+"_ES.txt";
	}
	else{
	  fileName=finalName+"_ES.txt";
	}
	ofs.close();
	ofs.open(fileName.c_str());
	for(int i=0;i<spec.size();++i){
	  for(int m=0;m<spec[i].size();++m){
	    ofs<<spec[i][m]<<"\t";
	  }
	  ofs<<std::endl;
	}
      }
      ofs.close();
    }
  }
  TensorNetwork.resetConvergence();
}
