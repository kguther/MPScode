#include "interface.h"
#include <stdlib.h>
#include "Qsystem.h"
#include "localHSpaces.h"
#include "delta.h"
#include "simulation.h"
#include "localMpo.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <math.h>
#include <vector>

interface::interface(){
  //default parameters which are used if a parameter is not specified in the file
  fileName="testRun";
  parPack.nSweeps=12;
  parPack.alphaInit=1e-3;
  parPack.arpackTol=1e-4;
  parPack.arpackTolMin=1e-8;
  parPack.L=100;
  parPack.N=parPack.L;
  parPack.scaling=100;
  //note that D=1 (or any other too small value for D) uses a fixed minimal value instead
  parPack.D=1;
  parPack.par=1;
  parPack.nGs=0;
  parPack.gsc=0;
  parPack.Jsc=-1;
  parPack.odd=0;
  parPack.rho=0.5;
  parPack.simType=0;
  parPack.nEigens=1;
  parPack.numPts=5;
  parPack.alphaMin=0;
  parPack.alphaMax=2*M_PI;
  parPack.Wsc=1;
}

//-------------------------------------------------------------------------------------------//

void interface::provideInterface(char *argv){
  /*
  char inArg;
  //choose mode
  std::cout<<"Determine correlations for a broad parameter range (c), for a single point (p) or calculate gap scaling (s)? ";
  std::cin>>inArg;
  if(inArg=='c'){
    parPack.simType=0;
  }
  if(inArg=='s'){
    parPack.simType=1;
  }
  if(inArg=='p'){
    parPack.simType=2;
  }
  std::string fN;
  std::cout<<"Enter filename ((d) for default file): ";
  std::cin>>fN;
  */
  std::string fN;
  if(argv){
    fN=argv;
    readParFile(fN);
  }
}

//-------------------------------------------------------------------------------------------//

void interface::readParFile(std::string const &fN){
  std::ifstream ifs;
  char inArg;
  char *target=new char[1000];
  int intPar;
  double fPar;
  ifs.open(fN.c_str());
  ifs.getline(target,1000);
  fileName=target;
  while(ifs.get(inArg)){
    if(inArg!=' '){
      if(inArg=='L' || inArg=='D' || inArg=='S' || inArg=='p' || inArg=='N' || inArg=='n' || inArg=='o' || inArg=='T' || inArg=='s' || inArg=='E'){
	ifs>>intPar;
	if(inArg=='L'){
	  parPack.L=intPar;
	}
	if(inArg=='D'){
	  parPack.D=intPar;
	}
	if(inArg=='S'){
	  parPack.nSweeps=intPar;
	}
	if(inArg=='N'){
	  parPack.N=intPar;
	}
	if(inArg=='p'){
	  parPack.par=intPar;
	}
	if(inArg=='o'){
	  parPack.odd=intPar;
	}
	if(inArg=='n'){
	  parPack.numPts=intPar;
	}
	if(inArg=='T'){
	  parPack.simType=intPar;
	}
	if(inArg=='s'){
	  parPack.nGs=intPar;
	}
	if(inArg=='E'){
	  parPack.nEigens=intPar;
	}
      }
      else{
	ifs>>fPar;
	if(inArg=='a'){
	  parPack.alphaInit=fPar;
	}
	if(inArg=='m'){
	  parPack.arpackTolMin=fPar;
	}
	if(inArg=='t'){
	  parPack.arpackTol=fPar;
	}
	if(inArg=='J'){
	  parPack.Jsc=fPar;
	}
	if(inArg=='g'){
	  parPack.gsc=fPar;
	}
	if(inArg=='r'){
	  parPack.rho=fPar;
	}
	if(inArg=='i'){
	  parPack.alphaMin=fPar;
	}
	if(inArg=='f'){
	  parPack.alphaMax=fPar;
	}
	if(inArg=='z'){
	  parPack.scaling=fPar;
	}
	if(inArg=='W'){
	  parPack.scaling=fPar;
	}
      }
      ifs.get(inArg);
    }
  }
  ifs.close();
  if(parPack.par!=1 && parPack.par!=-1){
    std::cout<<"Invalid parity supplied. Terminating process\n";
    exit(3);
  }
  if(parPack.N>2*parPack.L || parPack.N<0 || ((parPack.N==0 || parPack.N==2*parPack.L) && parPack.par==-1)){
    std::cout<<"Invalid particle number supplied. Terminating process\n";
    exit(3);
  }
}

//-------------------------------------------------------------------------------------------//

void interface::getScalingSerial(double J, double g){
  std::string dir="results/scaling/";
  std::vector<double> energiesGs, energiesEx, accsEx, accsGs, sizes;
  int const nEigens=2;
  int const L0=30;
  int const LMax=35;
  int N=0;
  parPack.numPts=1;
  int const nQuantumNumbers=1;
  int minimalD=(2*parPack.N>4)?2*parPack.N:4;
  int usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(parPack.N,parPack.par)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,1,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,1,1e-3,1e-4,1e-7,1e-4);
  simulation sim;
  for(int L=L0;L<=LMax;L+=5){
    N=L*parPack.rho+parPack.odd*(static_cast<int>((L*parPack.rho+1))%2)+(1-parPack.odd)*static_cast<int>(L*parPack.rho)%2;
    minimalD=(3*N>4)?3*N:4;
    usedD=(parPack.D>minimalD)?parPack.D:minimalD;
    simPars.D=usedD;
    QNValue[0]=std::complex<int>(N,parPack.par);
    pars.QNconserved=QNValue;
    pars.L=L;
    sim.generate(pars,simPars,parPack.Jsc,parPack.gsc,parPack.Wsc,parPack.numPts,parPack.scaling,fileName);
    sim.run();
    sizes.push_back(L);
    energiesGs.push_back(sim.E0[0]);
    energiesEx.push_back(sim.E0[1]);
    accsGs.push_back(sim.dE[0]);
    accsEx.push_back(sim.dE[1]);
  }
  std::ofstream ofs;
  std::ostringstream compositeName;
  compositeName<<dir<<fileName<<"_rho_"<<parPack.rho<<"_par_"<<parPack.par<<"_odd_"<<parPack.odd<<".txt";
  std::string finalName=compositeName.str();
  for(int m=0;m<finalName.length()-4;++m){
    if(finalName[m]=='.'){
      finalName.erase(m,1);
    }
  }
  ofs.open(finalName.c_str());
  ofs<<"System size\tGS energy\t excited state energy\t GS accuracy\t excited state accuracy\n";
  for(int m=0;m<energiesGs.size();++m){
    ofs<<sizes[m]<<"\t"<<energiesGs[m]<<"\t"<<energiesEx[m]<<"\t"<<accsGs[m]<<"\t"<<energiesGs[m]<<"\n";
  }
}
