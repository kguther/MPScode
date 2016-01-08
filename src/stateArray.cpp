#include "stateArray.h"
#include "arrayprocessing.h"
#include "arraycreation.h"

stateArray::stateArray(){
  dimInfo.initialize(1,1,1);
  createStateArray(1,1,&state_array_access_structure);
}

//---------------------------------------------------------------------------------------------------//

stateArray::stateArray(dimensionTable &dimInfoIn){
  initialize(dimInfoIn);
}

//---------------------------------------------------------------------------------------------------//

stateArray::~stateArray(){
  deleteStateArray(&state_array_access_structure);
}

//---------------------------------------------------------------------------------------------------//

void stateArray::mpsCpy(stateArray &source){
  deleteStateArray(&state_array_access_structure);
  initialize(source.dimInfo);
  int lDL, lDR, ld;
  for(int i=0;i<L;++i){
    lDL=locDimL(i);
    lDR=locDimR(i);
    ld=locd(i);
    for(int si=0;si<ld;++si){
      for(int ai=0;ai<lDR;++ai){
	for(int aim=0;aim<lDL;++aim){
	  state_array_access_structure[i][si][ai][aim]=source.global_access(i,si,ai,aim);
	}
      }
    }
  }
}


//---------------------------------------------------------------------------------------------------//

void stateArray::generate(dimensionTable &dimInfoIn){
  deleteStateArray(&state_array_access_structure);
  initialize(dimInfoIn);
}

//---------------------------------------------------------------------------------------------------//

void stateArray::initialize(dimensionTable &dimInfoIn){
  D=dimInfoIn.D();
  L=dimInfoIn.L();
  dimInfo=dimInfoIn;
  createStateArray(D,L,&state_array_access_structure);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::setParameterD(int Dnew){
  if(Dnew<D){
    return -1;
  }
  lapack_complex_double ****newNetworkState;
  dimensionTable backupDimInfo=dimInfo;
  //Copy the content of the current state into the larger array (which is initialized with zero)
  dimInfo.setParameterD(Dnew);
  createStateArray(Dnew,L,&newNetworkState);
  int lDL, lDR, ld;
  for(int i=0;i<L;++i){
    lDL=backupDimInfo.locDimL(i);
    lDR=backupDimInfo.locDimR(i);
    ld=backupDimInfo.locd(i);
    for(int si=0;si<ld;++si){
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

void stateArray::createStateArray(int const D, int const L, lapack_complex_double *****array){
  int dimR, dimL, lD, rD, dimd;
  dimR=0;
  dimL=0;
  dimd=0;
  for(int i=0;i<L;i++){
    dimR+=locd(i)*locDimR(i);
    dimL+=locd(i)*locDimR(i)*locDimL(i);
    dimd+=locd(i);
  }
  (*array)=new lapack_complex_double ***[L];
  (*array)[0]=new lapack_complex_double **[dimd];
  for(int i=1;i<L;++i){
    (*array)[i]=(*array)[i-1]+locd(i);
  }
  (*array)[0][0]=new lapack_complex_double*[dimR];
  for(int i=0;i<L;i++){
    rD=locDimR(i);
    if(i>0){
      (*array)[i][0]=(*array)[i-1][0]+locd(i)*locDimR(i-1);
    }
    for(int si=1;si<locd(i);++si){
      (*array)[i][si]=(*array)[i][si-1]+rD;
    }
  }
  (*array)[0][0][0]=new lapack_complex_double[dimL];
  for(int i=0;i<L;i++){
    lD=locDimL(i);
    rD=locDimR(i);
    if(i>0){
      (*array)[i][0][0]=(*array)[i-1][0][0]+locd(i)*locDimL(i-1)*locDimR(i-1);
    }
    for(int si=0;si<locd(i);++si){
      if(si>0){
	(*array)[i][si][0]=(*array)[i][si-1][0]+lD*rD;
      }
      for(int ai=1;ai<rD;++ai){
	(*array)[i][si][ai]=(*array)[i][si][ai-1]+lD;
      }
    }
  }
  for(int i=0;i<L;i++){
    lD=locDimL(i);
    rD=locDimR(i);
    for(int si=0;si<locd(i);si++){
      for(int ai=0;ai<rD;ai++){
	for(int aim=0;aim<lD;aim++){
	  (*array)[i][si][ai][aim]=0;
	}
      }
    }
  }
}
