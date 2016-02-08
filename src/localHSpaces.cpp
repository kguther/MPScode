#include "localHSpaces.h"

localHSpaces::localHSpaces(){
}

localHSpaces::localHSpaces(int const d):
  constantDimension(1)
{
  localHDims.push_back(d);
}

localHSpaces::localHSpaces(std::vector<int> d):
  constantDimension(0),
  localHDims(d)
{}

int localHSpaces::locd(int const i) const{
  if(constantDimension){
    return localHDims[0];
  }
  return localHDims[i];
}

int localHSpaces::maxd() const{
  int dmax=1;
  for(int i=0;i<localHDims.size();++i){
    if(localHDims[i]>dmax){
      dmax=localHDims[i];
    }
  }
  return dmax;
}
