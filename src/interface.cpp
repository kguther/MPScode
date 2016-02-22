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
	  parPack.Wsc=fPar;
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
