#include "stateArray.h"
#include "arrayprocessing.h"
#include "arraycreation.h"

stateArray::stateArray(){
  createStateArray(1,1,1,&state_array_access_structure);
}

//---------------------------------------------------------------------------------------------------//

stateArray::stateArray(int const din, int const Din, int const Lin){
  initialize(din,Din,Lin);
}

//---------------------------------------------------------------------------------------------------//

stateArray::~stateArray(){
  deleteStateArray(&state_array_access_structure);
}

//---------------------------------------------------------------------------------------------------//

void stateArray::mpsCpy(stateArray &source){
  deleteStateArray(&state_array_access_structure);
  initialize(source.siteDim(),source.maxDim(),source.length());
  int lDL, lDR;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    for(int si=0;si<d;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  state_array_access_structure[i][si][ai][aim]=source.global_access(i,si,ai,aim);
	}
      }
    }
  }
}


//---------------------------------------------------------------------------------------------------//

void stateArray::generate(int din, int Din, int Lin){
  deleteStateArray(&state_array_access_structure);
  initialize(din,Din,Lin);
}

//---------------------------------------------------------------------------------------------------//

void stateArray::initialize(int din, int Din, int Lin){
  d=din;
  D=Din;
  L=Lin;
  createStateArray(d,D,L,&state_array_access_structure);
  getIcrit();
}

//---------------------------------------------------------------------------------------------------//

int stateArray::setParameterD(int Dnew){
  if(Dnew<D){
    return -1;
  }
  lapack_complex_double ****newNetworkState;
  //Copy the content of the current state into the larger array (which is initialized with zero)
  createStateArray(d,Dnew,L,&newNetworkState);
  int lDL, lDR;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    for(int si=0;si<d;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  newNetworkState[i][si][ai][aim]=state_array_access_structure[i][si][ai][aim];
	}
      }
    }
  }
  //Replace the current state
  deleteStateArray(&state_array_access_structure);
  state_array_access_structure=newNetworkState;
  D=Dnew;
  getIcrit();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

void stateArray::getIcrit(){
  icrit=L/2;//In case chain is too short, this is the correct value (trust me)
  for(int j=0;j<L/2;++j){
    if(locDMax(j)>D){
      icrit=j;
      break;
    }
  }
}

//---------------------------------------------------------------------------------------------------//
// These functions only return the column/row dimension of the matrices of the i-th site
// Simple, but essential. Using local matrix dimensions make lots of stuff easier.
// NAMING CONVENTION: Any variable-name starting with l indicates a local dimension
// NAMING CONVENTION: Any variable-name of a local dimension ending with R is a local row dimension,
// any ending with L is a local column dimension
//---------------------------------------------------------------------------------------------------//

int stateArray::locDimL(int const i){
  if(i<=icrit){
    return locDMax(i-1);
  }
  if(i<=L-icrit-1){
    return D;
  }
  return locDMax(i);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDimR(int const i){
  if(i<icrit){
    return locDMax(i);
  }
  if(i<=L-icrit-2){
    return D;
  }
  return locDMax(i+1);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int stateArray::locd(int const i){
  return d;
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDMax(int const i){
  if(i<=L/2){
    return pow(d,i+1);
  }
  return pow(d,L-i);
}
