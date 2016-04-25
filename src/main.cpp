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

void sysSolve(info const &parPack, std::string const &fileName, std::vector<double> &energies);
void getScaling(int L, info const &parPack, double *results, std::string const &fileName);
void sysSetMeasurements(simulation &sim, int d, int L, int meas);
void getFileName(info const &necPars, char *fNBuf, int commsize, int myrank, std::string &finalName);
//results has to be at least of size 4 (in the sense of a C array)

int main(int argc, char *argv[]){
  //Here, the parameters are distributed via MPI to the processes. Each process then individuall solves the system for a specific set of parameters - great paralellization.
  //There are currently two settings: scaling and correlation. The former computes the behaivour of the gap with increasing system size and the latter computes correlations etc across the parameter space for fixed system size
  MPI_Init(&argc,&argv);
  int myrank, commsize;
  info necPars;
  char *fNBuf;
  int fNBufSize;
  int const dn=2;
  int const L0=30;
  int const dL=5;
  double alpha;
  int L;
  MPI_Datatype apparentTypes[dn],mpiInfo;
  int blockLengths[dn];

  //Defines the MPI_Datatype used for transferring the input parameters
  MPI_Aint displacements[dn], firstAdress, secondAdress;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&commsize);
  blockLengths[0]=13;
  blockLengths[1]=14;
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
  std::string finalName;
  //Each process calculates its own couplings/system size
  

  int const range=4;
  int const redRank=myrank/2;
  int WStage=redRank/range;
  if(necPars.simType==0){
    necPars.Wsc+=(1+WStage)/2*0.2*pow(-1,WStage);
    necPars.par=2*(myrank%2)-1;
  }
  alpha=(necPars.alphaMax-necPars.alphaMin)*((redRank%range)/static_cast<double>(range))+necPars.alphaMin;
  if(necPars.simType==0){
    necPars.Jsc=cos(alpha);
    necPars.gsc=sin(alpha);
  }

  std::vector<double> energies;
  L=L0+dL*myrank;
  //And evaluates its computation
  if(necPars.simType==1){
    for(int pty=-1;pty<2;pty+=2){
      for(int varL=35;varL<110;varL+=3){
	necPars.L=varL;
	necPars.N=necPars.rho*necPars.L*2;
	necPars.par=pty;
	getFileName(necPars,fNBuf,commsize,myrank,finalName);
	sysSolve(necPars,finalName,energies);
      }
    }
    /*
    double *energies=new double[4*commsize];
    double results[4];
    //getScaling(L,necPars,results,finalName);
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
    }
    delete[] energies;
    */
  }
 
  if(necPars.simType==3){
    necPars.nEigens=2;
    /*
    switch(myrank){
    case 0:
      necPars.tImag=0;
      necPars.tPos=0;
      break;
    case 1:
      necPars.tReal=0;
      necPars.tPos=0;
      break;
    case 2:
      necPars.tImag=0;
      necPars.tPos=necPars.L/2;
      break;
    case 3:
      necPars.tReal=0;
      necPars.tPos=necPars.L/2;
      break;
    case 4:
      necPars.tReal=0;
      necPars.tPos=-1;
      break;
    default:
      necPars.tImag=0;
      necPars.tPos=-1;
      }
    */
    necPars.tReal=0;
    necPars.tPos=myrank+12;
  }
  //The output filename is generated
  if(necPars.simType!=1){
    getFileName(necPars,fNBuf,commsize,myrank,finalName);
  }
  if(necPars.simType==2 || necPars.simType==0){
    sysSolve(necPars,finalName,energies);
  }
  if(necPars.simType==3){
    necPars.numPts=1;
    sysSolve(necPars,finalName,energies);
    energies.resize(4);
    double *energyBuf=new double[4*commsize];
    double results[4];
    results[0]=energies[0];
    results[1]=energies[1];
    results[2]=energies[2];
    results[3]=energies[3];
    //rcvcount is the number of objects recieved PER PROCESS, not in total
    MPI_Gather(results,4,MPI_DOUBLE,energyBuf,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //For scaling, all results are written in the same file. This is done by the poor man's solution: only the main process writes.
    if(myrank==0){
      std::ofstream ofs;
      ofs.open("SB_local_sweep_B.txt");
      ofs<<"Parameters: J="<<necPars.Jsc<<" g="<<necPars.gsc<<std::endl;
      ofs<<"System size: "<<necPars.L<<" filling: "<<necPars.rho<<" subchain parity="<<necPars.par<<std::endl;
      ofs<<"SB Position\tGS energy\t excited state energy\t GS accuracy\t excited state accuracy\n";
      for(int rk=0;rk<commsize;++rk){
	ofs<<rk<<"\t"<<energyBuf[4*rk]<<"\t"<<energyBuf[4*rk+1]<<"\t"<<energyBuf[4*rk+2]<<"\t"<<energyBuf[4*rk+3]<<std::endl;
      }
      ofs.close();
    }
    delete[] energyBuf;

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
  std::complex<int> QNValue[1]={std::complex<int>(N,parPack.par)};
  std::complex<int> QNList[4]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};
  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,L,parPack.Dw,nEigens,nQuantumNumbers,QNValue,QNList);
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

void sysSolve(info const &parPack, std::string const &fileName, std::vector<double> &energies){
  int const nQuantumNumbers=1;
  int minimalD=(2*parPack.N>4)?2*parPack.N:4;
  int usedD=(parPack.D>minimalD)?parPack.D:minimalD;
  std::complex<int> QNValue[1]={std::complex<int>(parPack.N,parPack.par)};
  std::complex<int> QNList[4]={std::complex<int>(0,1),std::complex<int>(1,1),std::complex<int>(1,-1),std::complex<int>(2,-1)};

  localHSpaces localHilbertSpaceDims(4);
  problemParameters pars(localHilbertSpaceDims,parPack.L,parPack.Dw,parPack.nEigens,nQuantumNumbers,QNValue,QNList,parPack.tReal,parPack.tImag);
  //Arguments of simPars: D, NSweeps, NStages, alpha (initial value), accuracy threshold, minimal tolerance for arpack, initial tolerance for arpack
  simulationParameters simPars(usedD,parPack.nSweeps,parPack.nStages,parPack.alphaInit,parPack.acc,parPack.arpackTolMin,parPack.arpackTol);
  simulation sim(pars,simPars,parPack.Jsc,parPack.gsc,parPack.Wsc,parPack.numPts,parPack.scaling,parPack.delta,fileName,parPack.tPos);
  if(parPack.simType==2){
    int const d=pars.d.maxd();
    int const L=parPack.L;
    localMpo<lapack_complex_double> gamma(d,1,L,L/4,0);
    std::string gammaName="Second Order ";
    std::string fName;
    std::ostringstream cGName;
    double theta;
    for(int m=0;m<parPack.nGs;++m){
      theta=m*2*M_PI/static_cast<double>(parPack.nGs);
      writePhasedSecondOrder(gamma,theta);
      cGName<<gammaName<<theta;
      fName=cGName.str();
      sim.setLocalMeasurement(gamma,fName);
      cGName.str("");
    }
  }
  int const doMeas=(parPack.simType==1 || parPack.simType==3 || parPack.simType==0)?0:1;
  sysSetMeasurements(sim,pars.d.maxd(),parPack.L,doMeas);
  energies.resize(2*sim.E0.size());
  for(int iE=0;iE<sim.E0.size();++iE){
    energies[iE]=sim.E0[iE];
    energies[iE+sim.E0.size()]=sim.dE[iE];
  }
}

//-------------------------------------------------------------------------------------------//

void sysSetMeasurements(simulation &sim, int d, int L, int meas){
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
  localMpo<lapack_complex_double> bulkICSuperConductingCorrelation(d,1,L,bulkStart,0,2);
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
  std::string pscName="Bulk interchain superconducting corrleation";
  std::string picscName="Bulk intrachain superconducting corrleation";
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
	bulkICSuperConductingCorrelation.global_access(i,si,sip,0,0)=delta(si,sip);
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
      bulkSuperConductingCorrelation.global_access(bulkStart-3,si,sip,0,0)=bMatrix(sip,si)*(delta(sip,0)-delta(sip,1));
      bulkSuperConductingCorrelation.global_access(bulkStart-2,si,sip,0,0)=bMatrix(sip,si);
      bulkSuperConductingCorrelation.global_access(bulkStart-1,si,sip,0,0)=aMatrix(si,sip)*(delta(sip,1)-delta(sip,3));
      bulkSuperConductingCorrelation.global_access(bulkStart,si,sip,0,0)=aMatrix(si,sip);
      bulkICSuperConductingCorrelation.global_access(bulkStart-3,si,sip,0,0)=aMatrix(sip,si)*(delta(sip,0)-delta(sip,2));
      bulkICSuperConductingCorrelation.global_access(bulkStart-2,si,sip,0,0)=aMatrix(sip,si);
      bulkICSuperConductingCorrelation.global_access(bulkStart-1,si,sip,0,0)=aMatrix(si,sip)*(delta(sip,1)-delta(sip,3));
      bulkICSuperConductingCorrelation.global_access(bulkStart,si,sip,0,0)=aMatrix(si,sip);
    }
  }

  
  sim.setLocalMeasurement(localDensity,lDName);
  sim.setLocalMeasurement(greensFunction,gFName);
  sim.setLocalMeasurement(bulkGreensFunction,bgFName);
  if(meas){
    sim.setLocalMeasurement(localDensityB,lDOName);
    sim.setLocalMeasurement(localDensityProd,lDPName);
    sim.setLocalMeasurement(interChainCorrelation,iCCName);
    sim.setLocalMeasurement(densityCorrelation,dCName);
    sim.setLocalMeasurement(interChainDensityCorrelation,iCDCName);
    sim.setLocalMeasurement(superconductingOrder,scName);
    sim.setLocalMeasurement(bulkInterChainCorrelation,biCCName);
    sim.setLocalMeasurement(bulkDensityCorrelation,bdCName);
    sim.setLocalMeasurement(bulkInterChainDensityCorrelation,biCDCName);
    sim.setLocalMeasurement(bulkSuperconductingOrder,bscName);
    sim.setLocalMeasurement(bulkSuperConductingCorrelation,pscName);
    sim.setLocalMeasurement(bulkICSuperConductingCorrelation,picscName);
  }
  sim.setEntanglementSpectrumMeasurement();
  
  sim.run();
}

//-------------------------------------------------------------------------------------------//

void getFileName(info const &necPars, char *fNBuf, int commsize, int myrank, std::string &finalName){
  std::string dir="results/";
  std::string type;
  if(necPars.simType==1){
    type="_scaling";
  }
  if(necPars.simType==0){
    type="_correlations";
  }
  if(necPars.simType==2){
    type="_point";
  }
  if(necPars.simType==3){
    type="_position_scan_SB_2";
    /*
    switch(myrank){
    case 0:
    case 1:
      type="_0";
      break;
    case 2:
    case 3:
      type="_mid";
      break;
    default:
      type="_global";
    }
      */
  }
  std::ostringstream compositeName;
  compositeName<<dir<<fNBuf<<type;
  if(symmetryBroken(necPars)){
    compositeName<<"_single_hop_tR_"<<necPars.tReal<<"_tI_"<<necPars.tImag;
  }
  if(necPars.simType==1 || necPars.simType==3){
    //compositeName<<"_run_"<<myrank<<"_rho_"<<necPars.rho<<"_par_"<<necPars.par<<"_odd_"<<necPars.odd<<"_J_"<<necPars.Jsc<<"_g_"<<necPars.gsc<<".txt";
    compositeName<<"_run_"<<myrank<<"_L_"<<necPars.L<<"_N_"<<necPars.N<<"_p_"<<necPars.par;
  }
  else{
    if(necPars.simType==2 && commsize>1){
      compositeName<<"_run_"<<myrank;
    }
    compositeName<<"_L_"<<necPars.L<<"_N_"<<necPars.N<<"_p_"<<necPars.par;
  }
  finalName=compositeName.str();
  //Only type-1 runs do not use the simulation output, where the filename is generated. There, the name is generated here
  if(necPars.simType==1){
    for(int m=0;m<finalName.length()-4;++m){
      if(finalName[m]=='.'){
	finalName.erase(m,1);
      }
    }
  }
  std::cout<<finalName<<std::endl;
}
