#ifndef PARAMETERS
#define PARAMETERS

class parameters;

class parameters{
 public:
  parameters();
  parameters(int din, int Din, int Lin, int Dwin, int Nin, int nEigsin=1, int nStagesin=2);
  int D, L, nSweeps, d, Dw, nEigs, nStages;
};

#endif
