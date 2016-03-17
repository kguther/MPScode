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

void sysSolve(info const &parPack, std::string const &fileName);
void sysScan(double J, double g, info const &parPack, std::string const &fileName);
void getScaling(int L, info const &parPack, double *results, std::string const &fileName);
void sysSetMeasurements(simulation &sim, int d, int L);
//results has to be at least of size 4 (in the sense of a C array)

int main(int argc, char *argv[]){
  //Here, the parameters are distributed via MPI to the processes. Each process then individuall solves the system for a specific set of parameters - great paralellization.
  //There are currently two settings: scaling and correlation. The former computes the behaivour of the gap with increasing system size and the latter computes correlations etc across the parameter space for fixed system size
  MPI_Init(&argc,&argv);
  std::string dir="results/";
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

  //Defines the MPI_Datatype used for transferring the input parameters
  MPI_Aint displacements[dn], firstAdress, secondAdress;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&commsize);
  blockLengths[0]=10;
  blockLengths[1]=12;
  MPI_Get_address(&necPars.L,&firstAdress);
  MPI_Get_address(&necPars.rho,&secondAdress);
  displacements[0]=(MPI_Aint)0;
  displacements[1]=secondAdress-firstAdress;
  apparentTypes[0]=MPI_INT;
  apparentTypes[1]=MPI_DOUBLE;
  MPI_Type_create_struct(dn,blockLengths,displacements,apparentTypes,&mpiInfo);
  MPI_Type_commit(&mpiInfo);

  //The main process now gets the input and prepares the broadcast of parameters
  //The interface class mainly reads parameter files
  if(myrank==0){
    interface settings;
    if(argc==2){
      settings.provideInterface(argv[1]);
    }
    else{
      settings.provideInterface(0);
    }
    necPars=settings.parPack;
    fNBufSize=settings.fileName.size();
    fNBuf=new char[fNBufSize+1];
    for(int m=0;m<fNBufSize;++m){
      fNBuf[m]=settings.fileName[m];
    }
    //C-strings are null-terminated. This is really necessary, strange names can be given if omitted.
    fNBuf[fNBufSize]='\0';
  }
  //Here, the parameters are broadcasted
  MPI_Bcast(&necPars,1,mpiInfo,0,MPI_COMM_WORLD);
  MPI_Bcast(&fNBufSize,1,MPI_INT,0,MPI_COMM_WORLD);
  if(myrank!=0){
    fNBuf=new char[fNBufSize+1];
  }
  MPI_Bcast(fNBuf,fNBufSize+1,MPI_CHAR,0,MPI_COMM_WORLD);
  
  //The output filename is generated
  if(necPars.simType==1){
    type="_scaling";
  }
  if(necPars.simType==0){
      type="_correlations";
      necPars.par=2*(myrank%2)-1;
  }
  if(necPars.simType==2){
    type="_point";
  }
  std::ostringstream compositeName;
  compositeName<<dir<<fNBuf<<type;
  if(necPars.simType==1){
    compositeName<<"_rho_"<<necPars.rho<<"_par_"<<necPars.par<<"_odd_"<<necPars.odd<<"_J_"<<necPars.Jsc<<"_g_"<<necPars.gsc<<".txt";
  }
  else{
    compositeName<<"_L_"<<necPars.L<<"_N_"<<necPars.N<<"_p_"<<necPars.par;
  }
  std::string finalName=compositeName.str();
  //Only type-1 runs do not use the simulation output, where the filename is generated. There, the name is generated here
  if(necPars.simType==1){
    for(int m=0;m<finalName.length()-4;++m){
      if(finalName[m]=='.'){
	finalName.erase(m,1);
      }
    }
  }
  std::cout<<finalName<<std::endl;
  //Each process calculates its own couplings/system size
  int const range=4;
  int const redRank=myrank/2;
  alpha=(necPars.alphaMax-necPars.alphaMin)*((redRank%range)/static_cast<double>(range))+necPars.alphaMin;
  J=cos(alpha);
  g=sin(alpha);
  int WStage=redRank/range;
  if(necPars.simType==0){
    necPars.Wsc+=WStage*0.2*pow(-1,WStage);
  }
  L=L0+dL*myrank;
  //And evaluates its computation
  if(necPars.simType==1){
    double *energies=new double[4*commsize];
    double results[4];
    getScaling(L,necPars,results,finalName);
    //rcvcount is the number of objects recieved PER PROCESS, not in total
    MPI_Gather(results,4,MPI_DOUBLE,energies,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //For scaling, all results are written in the same file. This is done by the poor man's solution: only the main process writes.
    if(myrank==0){
      std::ofstream ofs;
      ofs.open(finalName.c_str());
      ofs<<"Parameters: J="<<necPars.Jsc<<" g="<<necPars.gsc<<std::endl;
      ofs<<"Filling: "<<necPars.rho<<" subchain parity="<<necPars.par<<std::endl;
      ofs<<"System size\tGS energy\t excited state energy\t GS accuracy\t excited state accuracy\n";
      for(int rk=0;rk<commsize;++rk){
	ofs<<L0+rk*dL<<"\t"<<energies[4*rk]<<"\t"<<energies[4*rk+1]<<"\t"<<energies[4*rk+2]<<"\t"<<energies[4*rk+3]<<std::endl;
      }
      ofs.close();
      delete[] energies;
    }
  }
  if(necPars.simType==0){
    sysScan(J,g,necPars,finalName);
  }
  if(necPars.simType==2){
    sysSolve(necPars,finalName);
  }
  MPI_Finalize();
  return 0;
}

//-------------------------------------------------------------------------------------------//

void getScaling(int L, info const &parPack, double *results, std::string const &fileName){
  int const nEigens=2;
  int const N=2*L*parPack.rho+parPack.odd*(static_cast<int>((2*L*parPack.rho+1))%2)+(1-parPack.odd)*static_cast<int>(2*L*parPack.rho)%2;
  int const nQuantumNumbers=1;
  int minimalD=(2*N>4)?2*N:4;
  int usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::ofstream ofs;
  std::complex<int> QNValue[1]={std::complex<int>(N,parPack.par)};
  std::complex<int> QNList[4]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,parPack.nStages,parPack.alphaInit,parPack.acc,parPack.arpackTolMin,parPack.arpackTol);
  simulation sim(pars,simPars,parPack.Jsc,parPack.gsc,parPack.Wsc,parPack.numPts,parPack.scaling,parPack.delta,fileName);
  sim.run();
  results[0]=sim.E0[0];
  results[1]=sim.E0[1];
  results[2]=sim.dE[0];
  results[3]=sim.dE[1];
}

//-------------------------------------------------------------------------------------------//

void sysScan(double J, double g, info const &parPack, std::string const &fileName){
  int const L=parPack.L;
  int const nEigens=parPack.nEigens;
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
  simulationParameters simPars(usedD,parPack.nSweeps,1,parPack.alphaInit,parPack.acc,parPack.arpackTolMin,parPack.arpackTol);
  simulation sim(pars,simPars,J,g,parPack.Wsc,parPack.numPts,parPack.scaling,parPack.delta,fileName);
  sysSetMeasurements(sim,pars.d.maxd(),L);
}

//-------------------------------------------------------------------------------------------//

void sysSolve(info const &parPack, std::string const &fileName){
  int const nEigens=parPack.nEigens;
  int const nQuantumNumbers=1;
  int const numPoints=1;
  int const nGs=parPack.nGs;
  int minimalD=(2*parPack.N>4)?2*parPack.N:4;
  int usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(parPack.N,parPack.par)};
  std::complex<int> QNList[8]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,parPack.L,12,nEigens,nQuantumNumbers,QNValue,QNList);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,parPack.nStages,parPack.alphaInit,parPack.acc,parPack.arpackTolMin,parPack.arpackTol);
  simulation sim(pars,simPars,parPack.Jsc,parPack.gsc,parPack.Wsc,numPoints,parPack.scaling,parPack.delta,fileName);

  int const d=pars.d.maxd();
  int const L=parPack.L;
  localMpo<lapack_complex_double> gamma(d,1,L,L/4,0);
  std::string gammaName="Second Order ";
  std::string fName;
  std::ostringstream cGName;
  double theta;
  for(int m=0;m<nGs;++m){
    theta=m*2*M_PI/static_cast<double>(nGs);
    std::cout<<theta<<std::endl;
    writePhasedSecondOrder(gamma,theta);
    cGName<<gammaName<<theta;
    fName=cGName.str();
    sim.setLocalMeasurement(gamma,fName);
    cGName.str("");
  }
  sysSetMeasurements(sim,d,L);
}

//-------------------------------------------------------------------------------------------//

void sysSetMeasurements(simulation &sim, int d, int L){
  int const bulkStart=(L/4>2)?L/4:3;
  int parityQNs[4]={1,-1,-1,1};
  localMpo<lapack_complex_double> greensFunction(d,1,L,1,parityQNs);
  localMpo<lapack_complex_double> densityCorrelation(d,1,L,1,0);
  localMpo<lapack_complex_double> localDensity(d,1,L,0,0);
  localMpo<lapack_complex_double> localDensityB(d,1,L,0,0);
  localMpo<lapack_complex_double> localDensityProd(d,1,L,1,0);
  localMpo<lapack_complex_double> interChainCorrelation(d,1,L,1,parityQNs);
  localMpo<lapack_complex_double> superconductingOrder(d,1,L,1,0);
  localMpo<lapack_complex_double> interChainDensityCorrelation(d,1,L,1,0);
  localMpo<lapack_complex_double> bulkGreensFunction(d,1,L,bulkStart,parityQNs);
  localMpo<lapack_complex_double> bulkDensityCorrelation(d,1,L,bulkStart,0);
  localMpo<lapack_complex_double> bulkInterChainCorrelation(d,1,L,bulkStart,parityQNs);
  localMpo<lapack_complex_double> bulkSuperconductingOrder(d,1,L,bulkStart,0);
  localMpo<lapack_complex_double> bulkInterChainDensityCorrelation(d,1,L,bulkStart,0);
  localMpo<lapack_complex_double> bulkSuperConductingCorrelation(d,1,L,bulkStart,0,2);
  std::string gFName="Intrachain correlation";
  std::string dCName="Intrachain density correlation";
  std::string lDName="Local density";
  std::string lDOName="Local density B";
  std::string lDPName="Local density product";
  std::string iCDCName="Interchain density correlation";
  std::string iCCName="Interchain hopping correlation";
  std::string scName="Interchain pairwise correlation";
  std::string bgFName="Bulk correlation function";
  std::string bdCName="Bulk intrachain density correlation";
  std::string biCDCName="Bulk interchain density correlation";
  std::string biCCName="Bulk interchain hopping correlation";
  std::string bscName="Bulk interchain pairwise correlation";  
  std::string pscName="Bulk superconducting corrleation";
  //Define some interesting operators in MPO representation. These are mostly correlation functions which are product operators and therefore have Dw=1
  for(int i=0;i<L;++i){
    for(int si=0;si<d;++si){
      for(int sip=0;sip<d;++sip){
	greensFunction.global_access(i,si,sip,0,0)=delta(si,sip);
        densityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	localDensity.global_access(i,si,sip,0,0)=delta(si,sip);
	localDensityB.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	superconductingOrder.global_access(i,si,sip,0,0)=delta(si,sip);
	interChainDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkGreensFunction.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkInterChainCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkSuperconductingOrder.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkInterChainDensityCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
	localDensityProd.global_access(i,si,sip,0,0)=delta(si,sip);
	bulkSuperConductingCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
      }
    }
  }
  for(int si=0;si<d;++si){
    for(int sip=0;sip<d;++sip){
      greensFunction.global_access(1,si,sip,0,0)=bMatrix(sip,si);
      greensFunction.global_access(0,si,sip,0,0)=bMatrix(si,sip)*(delta(sip,2)-delta(sip,3));
      densityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      densityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(1,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      interChainDensityCorrelation.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,2)+delta(si,3));
      localDensity.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,1)+delta(si,3));
      localDensityB.global_access(0,si,sip,0,0)=delta(si,sip)*(delta(si,2)+delta(si,3));
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
      localDensityProd.global_access(1,si,sip,0,0)=delta(si,sip)*delta(si,3);
      bulkSuperConductingCorrelation.global_access(bulkStart-3,si,sip,0,0)=aMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
      bulkSuperConductingCorrelation.global_access(bulkStart-2,si,sip,0,0)=aMatrix(sip,si);
      bulkSuperConductingCorrelation.global_access(bulkStart-1,si,sip,0,0)=aMatrix(si,sip)*(delta(sip,1)-delta(sip,3));
      bulkSuperConductingCorrelation.global_access(bulkStart,si,sip,0,0)=aMatrix(si,sip);
    }
  }
  sim.setLocalMeasurement(localDensity,lDName);
  sim.setLocalMeasurement(localDensityB,lDOName);
  sim.setLocalMeasurement(localDensityProd,lDPName);
  sim.setLocalMeasurement(greensFunction,gFName);
  sim.setLocalMeasurement(interChainCorrelation,iCCName);
  sim.setLocalMeasurement(densityCorrelation,dCName);
  sim.setLocalMeasurement(interChainDensityCorrelation,iCDCName);
  sim.setLocalMeasurement(superconductingOrder,scName);
  sim.setLocalMeasurement(bulkGreensFunction,bgFName);
  sim.setLocalMeasurement(bulkInterChainCorrelation,biCCName);
  sim.setLocalMeasurement(bulkDensityCorrelation,bdCName);
  sim.setLocalMeasurement(bulkInterChainDensityCorrelation,biCDCName);
  sim.setLocalMeasurement(bulkSuperconductingOrder,bscName);
  sim.setLocalMeasurement(bulkSuperConductingCorrelation,pscName);
  sim.setEntanglementSpectrumMeasurement();
  sim.run();
}
