#ifndef PARAMETER_CLASSES
#define PARAMETER_CLASSES

class problemParameters{
 public:
  problemParameters(){}
 problemParameters(int din, int Lin, int Dwin, int nEigsin=1): d(din),L(Lin),Dw(Dwin),nEigs(nEigsin){}
  int L, d, Dw, nEigs;
};

//---------------------------------------------------------------------------------------------------//

//THE CONVERGENCE CHECK MAY OR MAY NOT FAIL FOR EXCITED STATES, SINCE THEY CAN CONTAIN NOTABLE FROM THE GROUND STATE AT SOME POINT IN THE ALGORITHM. 
//THEREFORE, SET accin=0 (OR SOME VERY LOW VALUE) FOR EXCITED STATE SEARCH UNTIL FIXED
//POSSIBLE FIX: THRESHOLD VALUE FOR OVERLAP WITH GROUNDSTATE -> EXCLUDES STATES TOO SIMILAR TO THE GROUND STATE

class simulationParameters{
 public:
 simulationParameters(int Din=100, int Nin=4, int nStagesin=2, double admixture=1e-2 ,double accin=1e-8, double tolMinin=1e-8, double tolInitialin=1e-4): D(Din), nSweeps(Nin), nStages(nStagesin),devAccuracy(accin),tolMin(tolMinin),tolInitial(tolInitialin), alpha(admixture){}
  int D, nStages, nSweeps;
  double devAccuracy, tolMin, tolInitial, alpha;
};

#endif
