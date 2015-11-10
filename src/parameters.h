#ifndef PARAMETER_CLASSES
#define PARAMETER_CLASSES

class problemParameters{
 public:
  problemParameters(){}
 problemParameters(int din, int Lin, int Dwin, int nEigsin=1): d(din),L(Lin),Dw(Dwin),nEigs(nEigsin){}
  int L, d, Dw, nEigs;
};

class simulationParameters{
 public:
 simulationParameters(int Din=100, int Nin=4, int nStagesin=2, double admixture=1e-2 ,double accin=1e-8, double tolMinin=1e-8, double tolInitialin=1e-4): D(Din), nSweeps(Nin), nStages(nStagesin),devAccuracy(accin),tolMin(tolMinin),tolInitial(tolInitialin), alpha(admixture){}
  int D, nStages, nSweeps;
  double devAccuracy, tolMin, tolInitial, alpha;
};

#endif
