#include <math.h>
#include "dimensionTable.h"

dimensionTable::dimensionTable(){
}

dimensionTable::dimensionTable(int const din, int const Din, int const Lin):
  dpar(din),
  Dpar(Din),
  Lpar(Lin)
{
  getIcrit();
}

void dimensionTable::initialize(int const din, int const Din, int const Lin){
  dpar=din;
  Lpar=Lin;
  setParameterD(Din);
}

void dimensionTable::setParameterD(int const Dnew){
  Dpar=Dnew;
  getIcrit();
}

//---------------------------------------------------------------------------------------------------//

void dimensionTable::getIcrit(){
  icrit=Lpar/2;//In case chain is too short, this is the correct value (trust me)
  for(int j=0;j<Lpar/2;++j){
    if(locDMax(j)>Dpar){
      icrit=j;
      break;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDimL(int const i){
  if(i<=icrit){
    return locDMax(i-1);
  }
  if(i<=Lpar-icrit-1){
    return Dpar;
  }
  return locDMax(i);
}

//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDimR(int const i){
  if(i<icrit){
    return locDMax(i);
  }
  if(i<=Lpar-icrit-2){
    return Dpar;
  }
  return locDMax(i+1);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int dimensionTable::locd(int const i){
  return dpar;
}

//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDMax(int const i){
  if(i<=Lpar/2){
    //return i+2;
    return pow(dpar,i+1);
  }
  //return L-i;
  return pow(dpar,Lpar-i);
}
