#include "parameters.h"

parameters::parameters(){
}

parameters::parameters(int din, int Din, int Lin, int Dwin, int Nin){
  N=Nin;
  L=Lin;
  D=Din;
  d=din;
  Dw=Dwin;
}
