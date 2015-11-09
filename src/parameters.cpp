#include "parameters.h"

parameters::parameters(){
}

parameters::parameters(int din, int Din, int Lin, int Dwin, int Nin, int nEigsin, int nStagesin, double accin){
  nSweeps=Nin;
  L=Lin;
  D=Din;
  d=din;
  Dw=Dwin;
  nEigs=nEigsin;
  nStages=nStagesin;
  acc=accin;
}
