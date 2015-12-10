#include "stateArray.h"
#include "arrayprocessing.h"
#include "arraycreation.h"

stateArray::stateArray(){
  dimInfo.initialize(1,1,1);
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
  dimInfo.initialize(din,Din,Lin);
  createStateArray(d,D,L,&state_array_access_structure);
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
  dimInfo.setParameterD(Dnew);
  return 0;
}

//---------------------------------------------------------------------------------------------------//
// These functions only return the column/row dimension of the matrices of the i-th site
// Simple, but essential. Using local matrix dimensions make lots of stuff easier.
// NAMING CONVENTION: Any variable-name starting with l indicates a local dimension
// NAMING CONVENTION: Any variable-name of a local dimension ending with R is a local row dimension,
// any ending with L is a local column dimension
//---------------------------------------------------------------------------------------------------//

int stateArray::locDimL(int const i){
  return dimInfo.locDimL(i);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDimR(int const i){
  return dimInfo.locDimR(i);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int stateArray::locd(int const i){
  return dimInfo.locd(i);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDMax(int const i){
  return dimInfo.locDMax(i);
}

void stateArray::createStateArray(int const d, int const D, int const L, lapack_complex_double *****array){
  int dimR, dimL, lD, rD;
  dimR=0;
  dimL=0;
  for(int i=0;i<L;i++){
    dimR+=d*locDimR(i);
    dimL+=d*locDimR(i)*locDimL(i);
  }
  create2D(L,d,array);
  (*array)[0][0]=new lapack_complex_double*[dimR];
  for(int i=0;i<L;i++){
    rD=locDimR(i);
    if(i>0){
      (*array)[i][0]=(*array)[i-1][0]+d*locDimR(i-1);
    }
    for(int si=1;si<d;si++){
      (*array)[i][si]=(*array)[i][si-1]+rD;
    }
  }
  (*array)[0][0][0]=new lapack_complex_double[dimL];
  for(int i=0;i<L;i++){
    lD=locDimL(i);
    rD=locDimR(i);
    if(i>0){
      (*array)[i][0][0]=(*array)[i-1][0][0]+d*locDimL(i-1)*locDimR(i-1);
    }
    for(int si=0;si<d;si++){
      if(si>0){
	(*array)[i][si][0]=(*array)[i][si-1][0]+lD*rD;
      }
      for(int ai=1;ai<rD;ai++){
	(*array)[i][si][ai]=(*array)[i][si][ai-1]+lD;
      }
    }
  }
  for(int i=0;i<L;i++){
    lD=locDimL(i);
    rD=locDimR(i);
    for(int si=0;si<d;si++){
      for(int ai=0;ai<rD;ai++){
	for(int aim=0;aim<lD;aim++){
	  (*array)[i][si][ai][aim]=0;
	}
      }
    }
  }
}
