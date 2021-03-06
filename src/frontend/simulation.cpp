#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "simulation.h"
#include "projector.h"
#include "globalMeasurement.h"
#include "localMeasurementSeries.h"
#include "exceptionClasses.h"

simulation::simulation(problemParameters &parsIn, simulationParameters &simParsIn, double J, double g, double WIn, int pathPoints, double stepSize, double deltaPIn, std::string const &targetFile, int tSiteIn, int jgScaleIn):
  simPars(simParsIn),
  pars(parsIn),
  W(WIn),
  pathLength(pathPoints),
  filePrefix(targetFile),
  measureEE(0),
  measureES(0),
  deltaP(deltaPIn),
  targetDelta(deltaPIn),
  scaling(stepSize),
  tSite(tSiteIn),
  jgScale(jgScaleIn),
  parDirection(std::complex<double>(J,g)),
  csystem(Qsystem(pars,simPars)),
  particleNumber(mpo<mpsEntryType>(pars.d.maxd(),2,pars.L)),
  subChainParity(mpo<mpsEntryType>(pars.d.maxd(),1,pars.L))
{
 // simulation class can only use nStages=1 (by choice of algorithm)
  E0.resize(pars.nEigs);
  dE.resize(pars.nEigs);
  convergedEigens.resize(pars.nEigs);
  double matEls;
  int L=pars.L;
  if(abs(parDirection)>1e-20 && jgScale==1){
    //jgScale==1 scales (J,g) normalized from (1,1), jgScale==0 does keep (J,g) constant but requires at least pathLength=2, jgScale==2 keeps J constant and scales g in steps of 0.3333 and jgScale==3 does keep (J,g) and W constant
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
  csystem.TensorNetwork.check=&particleNumber;
  csystem.TensorNetwork.checkParity=&subChainParity;
  localMeasureTask.clear();
  measureTask.clear();

  std::cout<<"Number of points per thread: "<<pathLength<<std::endl;
}

//---------------------------------------------------------------------------------------------------//

void simulation::setMeasurement(mpo<mpsEntryType > &MPOperator, std::string &operatorName){
  //defines a measurement of some global operator MPOperator in MPO representation
  measureTask.push_back(MPOperator);
  operatorNames.push_back(operatorName);
}

//---------------------------------------------------------------------------------------------------//

void simulation::setLocalMeasurement(localMpo<mpsEntryType > &localMPOperator, std::string &localOperatorName){
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
  double backupW;
  for(int nRun=1;nRun<pathLength+1;++nRun){

    if(jgScale==2){
      parDirection.imag(3.0-(nRun-1)/3.0);
    }
    if(abs(parDirection)>1e-20 && jgScale==1){
      parDirection*=1.0/(scaling*abs(parDirection))*nRun;
    }
    if(!jgScale){
      backupW=W;
      W=W/static_cast<double>(pathLength-1)*(nRun-1);
    }
    
    //Now: disorder scaling
    
    if(pathLength>1 && 0){
      deltaP=(nRun-1)*targetDelta/(pathLength-1);
    }
    else{
      deltaP=targetDelta;
    }
    std::cout<<"Starting run "<<nRun<<" out of "<<pathLength<<" runs. Using parameters\n";
    try{
      singleRun();
    }
    catch(critical_error &err){
      if(!jgScale){
	W=backupW;
      }
    }

    if(!jgScale){
      W=backupW;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void simulation::singleRun(){
  //Get nEigs eigenstates for a fixed set of parameters J,g of the hamiltonian
  int hInfo;
  double J,g;
  //Containers for measurements
  std::vector<double> expectationValues;
  std::vector<std::vector<mpsEntryType > > localExpectationValues;

  csystem.TensorNetwork.resetState();
  if(jgScale<=1){
    J=jgScale+parDirection.real();
    g=jgScale+parDirection.imag();
  }
  else{
    J=parDirection.real();
    g=parDirection.imag();
  }
  std::cout<<"J="<<J<<" g="<<g<<" W="<<W<<" t="<<pars.t<<std::endl;
  if(pars.Dw==12){
    hInfo=writeHamiltonian(csystem.TensorNetwork,J,g,W,pars.t,deltaP,tSite);
  }
  if(pathLength!=1 && simPars.nStages!=1){
    std::cout<<"Invalid simulation parameters: Staging is disabled for type-0 runs. Aborting run.\n";
    throw critical_error();
  }
  if(pars.Dw!=12){
    std::cout<<"Invalid bond dimension for the construction of H. Aborting run.\n";
    throw critical_error();
  }

  int numParities=1;
  int parities[2]={pars.QNconserved[0].imag(),-pars.QNconserved[0].imag()};
  //If disorder is enabled, we want to check both parity sectors with the same set of disorder
  if(std::abs(targetDelta)>1e-10){
    numParities=2;
  }
  
  for(int pc=0;pc<numParities;++pc){
    //If both parities are to be considered, set it now
    pars.QNconserved[0].imag(parities[pc]);
    csystem.TensorNetwork.setQuantumNumber(pars.QNconserved,pars.QNLocalList);
    //Actual DMRG
    csystem.getGroundState();
    E0=csystem.E0;
    dE=csystem.dE;
    expectationValues.resize(measureTask.size());
    localExpectationValues.resize(localMeasureTask.size());
    //Measure previously set operators and write results into the result file
    std::ofstream ofs;
    std::string finalName, fileName;
    std::ostringstream compositeName;
    std::cout<<pars.QNconserved[0]<<std::endl;
    compositeName<<filePrefix<<"_W_"<<W<<"_J_"<<J<<"_g_"<<g;
    if(std::abs(targetDelta)>1e-10){
      compositeName<<"_delta_"<<deltaP;
    }
    
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
      std::cout<<finalName<<std::endl;
      
      ofs.open(fileName.c_str());
      ofs<<"Values for state number "<<iEigen<<" with energy "<<E0[iEigen]<<" and energy variance "<<dE[iEigen]<<std::endl;
      //The problem parameters are written into the first lines
      ofs<<"L\tN\tsubchain parity\tJ\tg\tW\tE\tvariance of energy\tt\n";
      ofs<<std::setprecision(15);
      ofs<<pars.L<<"\t"<<real(pars.QNconserved[0])<<"\t"<<imag(pars.QNconserved[0])<<"\t"<<J<<"\t"<<g<<"\t"<<W<<"\t"<<E0[iEigen]<<"\t"<<dE[iEigen]<<"\t"<<pars.t<<std::endl;
      ofs<<std::setprecision(5);
      //First, global measurements are performed (this is used rarely)
      if(localMeasureTask.size()>0 || measureTask.size()>0 || measureEE){
	std::cout<<"Measuring correlation functions\n";
	for(int iM=0;iM<measureTask.size();++iM){
	  csystem.measure(&measureTask[iM],expectationValues[iM],iEigen);
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
	  csystem.measureLocalOperators(&localMeasureTask[iM],localExpectationValues[iM],iEigen);
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
	std::cout<<"Computing entanglement"<<std::endl;
	if(measureEE){
	  csystem.TensorNetwork.getEntanglement(S,spec,iEigen);
	}
	std::cout<<"Got entanglement entropy"<<std::endl;
	if(localExpectationValues.size()>0 || measureEE){
	  for(int i=0;i<pars.L;++i){
	    for(int iM=0;iM<localMeasureTask.size();++iM){
	      if(i<localExpectationValues[iM].size()){
		ofs<<std::abs(localExpectationValues[iM][i])<<"\t";
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
    csystem.TensorNetwork.resetConvergence();
  }
}
