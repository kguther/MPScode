#include "interface.h"
#include "simulation.h"
#include "localMpo.h"
#include "localHSpaces.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <mpi.h>
#include <vector>

void sysSolve(double J, double g, info parPack, std::string const &fileName);

void getScaling(int L, info parPack, std::vector<double> &results, std::string const &fileName);

int main(int argc, char *argv[]){
  MPI_Init(&argc,&argv);
  std::string dir="test_results/";
  std::string type;
  int myrank, commsize;
  info necPars;
  char *fNBuf;
  int fNBufSize;
  int const dn=2;
  int const L0=30;
  int const dL=5;
  double J,g,alpha;
  int L;
  MPI_Datatype apparentTypes[dn],mpiInfo;
  int blockLengths[dn];
  MPI_Aint displacements[dn], firstAdress, secondAdress;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&commsize);
  blockLengths[0]=8;
  blockLengths[1]=6;
  MPI_Get_address(&necPars.L,&firstAdress);
  MPI_Get_address(&necPars.rho,&secondAdress);
  displacements[0]=(MPI_Aint)0;
  displacements[1]=secondAdress-firstAdress;
  apparentTypes[0]=MPI_INT;
  apparentTypes[1]=MPI_DOUBLE;
  MPI_Type_create_struct(dn,blockLengths,displacements,apparentTypes,&mpiInfo);
  MPI_Type_commit(&mpiInfo);
  if(myrank==0){
    interface settings;
    settings.provideInterface();
    necPars=settings.parPack;
    fNBufSize=settings.fileName.size();
    fNBuf=new char[fNBufSize];
    for(int m=0;m<fNBufSize;++m){
      fNBuf[m]=settings.fileName[m];
    }
  }
  MPI_Bcast(&necPars,1,mpiInfo,0,MPI_COMM_WORLD);
  MPI_Bcast(&fNBufSize,1,MPI_INT,0,MPI_COMM_WORLD);
  if(myrank!=0){
    fNBuf=new char[fNBufSize];
  }
  MPI_Bcast(fNBuf,fNBufSize,MPI_CHAR,0,MPI_COMM_WORLD);
  if(necPars.simType){
    type="_scaling";
  }
  else{
    type="_correlations";
  }
  std::ostringstream compositeName;
  compositeName<<dir<<fNBuf<<type<<"_rho_"<<necPars.rho<<"_par_"<<necPars.par<<"_odd_"<<necPars.odd;
  if(necPars.simType){
    compositeName<<"_J_"<<necPars.Jsc<<"_g_"<<necPars.gsc<<".txt";
  }
  std::string finalName=compositeName.str();
  if(necPars.simType){
    for(int m=0;m<finalName.length()-4;++m){
      if(finalName[m]=='.'){
	finalName.erase(m,1);
      }
    }
  }
  alpha=2*M_PI*(myrank/static_cast<double>(commsize));
  J=cos(alpha);
  g=sin(alpha);
  L=L0+dL*myrank;
  if(necPars.simType){
    std::vector<double> results, energies;
    getScaling(L,necPars,results,finalName);
    energies.resize(4*commsize);
    MPI_Gather(&results[0],4,MPI_DOUBLE,&energies[0],4*commsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(myrank==0){
      std::ofstream ofs;
      ofs.open(finalName.c_str());
      ofs<<"Parameters: J="<<necPars.Jsc<<" g="<<necPars.gsc<<std::endl;
      ofs<<"System size\tGS energy\t excited state energy\t GS accuracy\t excited state accuracy\n";
      for(int rk=0;rk<commsize;++rk){
	ofs<<L0+rk*dL<<"\t"<<energies[4*rk]<<"\t"<<energies[4*rk+1]<<"\t"<<energies[4*rk+2]<<"\t"<<energies[4*rk+3]<<std::endl;
      }
      ofs.close();
    }
  }
  else{
    sysSolve(J,g,necPars,finalName);
  }
  MPI_Finalize();
  return 0;
}

//-------------------------------------------------------------------------------------------//

void getScaling(int L, info parPack, std::vector<double> &results, std::string const &fileName){
  int const nEigens=2;
  int const N=2*L*parPack.rho+parPack.odd*(static_cast<int>((2*L*parPack.rho+1))%2)+(1-parPack.odd)*static_cast<int>(2*L*parPack.rho)%2;
  int const nQuantumNumbers=1;
  int minimalD=(2*N>4)?2*N:4;
  int usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::ofstream ofs;
  std::complex<int> QNValue[1]={std::complex<int>(N,parPack.par)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,1,parPack.alphaInit,1e-4,parPack.arpackTolMin,parPack.arpackTol);
  simulation sim(pars,simPars,parPack.Jsc,parPack.gsc,parPack.numPts,fileName);

  sim.run();
  results.clear();
  results.push_back(sim.E0[0]);
  results.push_back(sim.E0[1]);
  results.push_back(sim.dE[0]);
  results.push_back(sim.dE[1]);
}

//-------------------------------------------------------------------------------------------//

void sysSolve(double J, double g, info parPack, std::string const &fileName){
  int const L=parPack.L;
  int const nEigens=1;
  int const numPoints=1;
  int measureEdge=1;
  int measureBulk=1;
  //The required bond dimension for the perturbed system seems to be greater than that of the unperturbed system
  int const nQuantumNumbers=1;
  int const minimalD=(2*parPack.N>4)?2*parPack.N:4;
  int const usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(parPack.N,parPack.par)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,parPack.L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //simulationParameters simPars(100,5,2,1e-4,1e-8,1e-9,1e-2);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,1,parPack.alphaInit,1e-8,parPack.arpackTolMin,parPack.arpackTol);

  simulation sim(pars,simPars,J,g,numPoints,fileName);
  int parityQNs[4]={1,-1,-1,1};
  int const bulkStart=parPack.L/3;

  localMpo<lapack_complex_double> greensFunction(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> densityCorrelation(pars.d.maxd(),1,L,1,0);
  localMpo<lapack_complex_double> localDensity(pars.d.maxd(),1,L,1,0);
  localMpo<lapack_complex_double> interChainCorrelation(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> superconductingOrder(pars.d.maxd(),1,L,1,parityQNs);
  localMpo<lapack_complex_double> interChainDensityCorrelation(pars.d.maxd(),1,L,1,0);
  localMpo<lapack_complex_double> bulkGreensFunction(pars.d.maxd(),1,L,bulkStart,parityQNs);
  localMpo<lapack_complex_double> bulkDensityCorrelation(pars.d.maxd(),1,L,bulkStart,0);
  localMpo<lapack_complex_double> bulkInterChainCorrelation(pars.d.maxd(),1,L,bulkStart,parityQNs);
  localMpo<lapack_complex_double> bulkSuperconductingOrder(pars.d.maxd(),1,L,bulkStart,parityQNs);
  localMpo<lapack_complex_double> bulkInterChainDensityCorrelation(pars.d.maxd(),1,L,bulkStart,0);
  std::string gFName="Intrachain correlation";
  std::string dCName="Intrachain density correlation";
  std::string lDName="Local density";
  std::string iCDCName="Interchain density correlation";
  std::string iCCName="Interchain hopping correlation";
  std::string scName="Interchain pairwise correlation";
  std::string bgFName="Bulk correlation function";
  std::string bdCName="Bulk intrachain density correlation";
  std::string biCDCName="Bulk interchain density correlation";
  std::string biCCName="Bulk interchain hopping correlation";
  std::string bscName="Bulk interchain pairwise correlation";  
  //Define some interesting operators in MPO representation. These are mostly correlation functions which are product operators and therefore have Dw=1
  for(int i=0;i<L;++i){
    for(int si=0;si<pars.d.maxd();++si){
      for(int sip=0;sip<pars.d.maxd();++sip){
	greensFunction.global_access(i,si,sip,0,0)=delta(si,sip);
        densityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	localDensity.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	superconductingOrder.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkGreensFunction.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkInterChainCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkSuperconductingOrder.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkInterChainDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
      }
    }
  }
  for(int si=0;si<pars.d.maxd();++si){
    for(int sip=0;sip<pars.d.maxd();++sip){
      greensFunction.global_access(1,si,sip,0,0)=bMatrix(sip,si);
      greensFunction.global_access(0,si,sip,0,0)=bMatrix(si,sip)*(delta(sip,2)-delta(sip,3));
      densityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      densityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,2)+delta(si,3));
      localDensity.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainCorrelation.global_access(1,si,sip,0,0)=delta(si,1)*delta(sip,2);
      interChainCorrelation.global_access(0,si,sip,0,0)=delta(si,1)*delta(sip,2);
      superconductingOrder.global_access(1,si,sip,0,0)=delta(si,3)*delta(sip,0);
      superconductingOrder.global_access(0,si,sip,0,0)=delta(si,0)*delta(sip,3);
      bulkGreensFunction.global_access(bulkStart,si,sip,0,0)=bMatrix(sip,si);
      bulkGreensFunction.global_access(bulkStart-1,si,sip,0,0)=bMatrix(si,sip)*(delta(sip,2)-delta(sip,3));
      bulkDensityCorrelation.global_access(bulkStart,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      bulkDensityCorrelation.global_access(bulkStart-1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      bulkInterChainDensityCorrelation.global_access(bulkStart,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      bulkInterChainDensityCorrelation.global_access(bulkStart-1,si,sip,0,0)=delta(si,sip)*(delta(si,2)+delta(si,3));
      bulkInterChainCorrelation.global_access(bulkStart,si,sip,0,0)=delta(si,1)*delta(sip,2);
      bulkInterChainCorrelation.global_access(bulkStart-1,si,sip,0,0)=delta(si,1)*delta(sip,2);
      bulkSuperconductingOrder.global_access(bulkStart,si,sip,0,0)=delta(si,3)*delta(sip,0);
      bulkSuperconductingOrder.global_access(bulkStart-1,si,sip,0,0)=delta(si,0)*delta(sip,3);
    }
  }
  //The hamiltonian is subchain parity conserving, thus, the expectation vaue of interChainCorrelation is zero
  //The hamiltonian is particle number conserving, thus, the expectation value of superconductingOrder is zero
  sim.setLocalMeasurement(localDensity,lDName);
  if(measureEdge){
  sim.setLocalMeasurement(greensFunction,gFName);
  sim.setLocalMeasurement(interChainCorrelation,iCCName);
  sim.setLocalMeasurement(densityCorrelation,dCName);
  sim.setLocalMeasurement(interChainDensityCorrelation,iCDCName);
  sim.setLocalMeasurement(superconductingOrder,scName);
  }
  if(measureBulk){
  sim.setLocalMeasurement(bulkGreensFunction,bgFName);
  sim.setLocalMeasurement(bulkInterChainCorrelation,biCCName);
  sim.setLocalMeasurement(bulkDensityCorrelation,bdCName);
  sim.setLocalMeasurement(bulkInterChainDensityCorrelation,biCDCName);
  sim.setLocalMeasurement(bulkSuperconductingOrder,bscName);
  }
  clock_t curtime;
  curtime=clock();
  sim.run();
  curtime=clock()-curtime;
  std::cout<<"\nTotal simulation took "<<(float)curtime/CLOCKS_PER_SEC<<" seconds\n\n";
  cout<<setprecision(21);
  for(int mi=0;mi<nEigens;++mi){
    std::cout<<"Obtained energy of state "<<mi<<" as: "<<sim.E0[mi]<<std::endl;
  }
}

