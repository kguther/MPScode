#include <math.h>
#include <iostream>
#include "dimensionTable.h"

dimensionTable::dimensionTable(){
}

//---------------------------------------------------------------------------------------------------//

dimensionTable::dimensionTable(int const Din, int const Lin, localHSpaces din):
  dpars(din),
  Dpar(Din),
  Lpar(Lin)
{
  dpar=dpars.maxd();
  getDMaxTable();
  getIcrit();
}

//---------------------------------------------------------------------------------------------------//

void dimensionTable::initialize(int const Din, int const Lin, localHSpaces din){
  dpars=din;
  Lpar=Lin;
  dpar=dpars.maxd();
  getDMaxTable();
  setParameterD(Din);
}

void dimensionTable::getDMaxTable(){
  int lDM;
  for(int i=0;i<=Lpar;++i){
    lDM=1;
    if(i<Lpar/2+1){
      for(int j=1;j<=i;++j){
	lDM*=dpars.locd(j);
      }
      DMaxTable.push_back(lDM);
    }
    else{
      for(int j=Lpar;j>i;--j){
	lDM*=dpars.locd(j);
      }
      DMaxTable.push_back(lDM);
    }
  }
}

//---------------------------------------------------------------------------------------------------//

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
// Functions returning the local bond dimensions of a MPS with the given parameters at site i.
//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDimL(int const i) const{
  if(i<=icrit){
    return locDMax(i-1);
  }
  if(i<=Lpar-icrit-1){
    return Dpar;
  }
  return locDMax(i-1);
}

//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDimR(int const i) const{
  if(i<icrit){
    return locDMax(i);
  }
  if(i<=Lpar-icrit-2){
    return Dpar;
  }
  return locDMax(i);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int dimensionTable::locd(int const i) const{
  return dpars.locd(i);
}

//---------------------------------------------------------------------------------------------------//

int dimensionTable::locDMax(int const i) const{
  return DMaxTable[i+1];
}
