#include "quantumNumber.h"

quantumNumber::quantumNumber(){
}

void quantumNumber::initialize(int const din, int const Lin, int const Nin, int *QNlocin){
  d=din;
  L=Lin;
  N=Nin;
  QNloc=QNlocin;
  QNlocMax=QNloc[0];
  for(int si=1;si<d;++si){
    if(QNloc[si]>QNlocMax){
      QNlocMax=QNloc[si];
    }
  }
}

int quantumNumber::QNLabel(int const i, int const ai, int const lDR){
  if(i>=L/2){
    return ai;
  }
  return ai;
}

int quantumNumber::QNLabel(int const si){
  return QNloc[si];
}
