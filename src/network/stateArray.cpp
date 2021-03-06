#include "stateArray.h"

stateArray::stateArray(){
}

//---------------------------------------------------------------------------------------------------//

stateArray::stateArray(dimensionTable const &dimInfoIn){
  initialize(dimInfoIn);
}

//---------------------------------------------------------------------------------------------------//

void stateArray::initialize(dimensionTable const &dimInfoIn){
  dimInfo=dimInfoIn;
  D=dimInfoIn.D();
  L=dimInfoIn.L();
  createStateArray(L);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::setParameterD(int Dnew){
  //Copy the content of the current state into the larger array (which is initialized with zero)
  dimInfo.setParameterD(Dnew);
  int info=0;
  std::vector<int> newDims(3,0);
  for(int i=0;i<L;++i){
    newDims[0]=locd(i);
    newDims[1]=locDimR(i);
    newDims[2]=locDimL(i);
    info+=stateArrayAccessStructure[i].setParameterDims(newDims);
  }
  //Replace the current state
  D=Dnew;
  return info;
}

//---------------------------------------------------------------------------------------------------//

int stateArray::setParameterL(int Lnew){
  dimInfo.setParameterL(Lnew);
  std::vector<int> locDims(3,0);
  int startPos;
  int const deltaL=Lnew-L;

  if(deltaL==0){
    return 0;
  }
  if(deltaL<0){
    startPos=L/2+deltaL/2;
    int endPos=L/2-deltaL/2+deltaL%2;
    stateArrayAccessStructure.erase(stateArrayAccessStructure.begin()+startPos,stateArrayAccessStructure.begin()+endPos);
  }
  else{
    startPos=L/2-deltaL/2;
    baseTensor<mpsEntryType > dummyTensor;
    std::vector<int> localDims(3,0);
    int cPos;
    for(int dI=0;dI<deltaL;++dI){
      cPos=startPos+dI+1;
      localDims[0]=locd(cPos);
      localDims[1]=locDimR(cPos);
      localDims[2]=locDimL(cPos);
      dummyTensor=baseTensor<mpsEntryType >(localDims);
      stateArrayAccessStructure.insert(stateArrayAccessStructure.begin()+cPos,dummyTensor);
    }
  }

  for(int i=0;i<dimInfo.L();++i){
    //std::cout<<stateArrayAccessStructure[i].lDL()<<"x"<<stateArrayAccessStructure[i].lDR()<<std::endl;
  }
  
  L=Lnew;
  return 0;
}

//---------------------------------------------------------------------------------------------------//

std::vector<int> stateArray::getIndexVec(int si, int ai, int aim) const{
  std::vector<int> result(3,0);
  result[0]=si;
  result[1]=ai;
  result[2]=aim;
  return result;
}

//---------------------------------------------------------------------------------------------------//
// These functions only return the column/row dimension of the matrices of the i-th site
// Simple, but essential. Using local matrix dimensions make lots of stuff easier.
// NAMING CONVENTION: Any variable-name starting with l indicates a local dimension
// NAMING CONVENTION: Any variable-name of a local dimension ending with R is a local row dimension,
// any ending with L is a local column dimension
//---------------------------------------------------------------------------------------------------//

int stateArray::locDimL(int i) const{
  return dimInfo.locDimL(i);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDimR(int i) const{
  return dimInfo.locDimR(i);
}

//---------------------------------------------------------------------------------------------------//
// These are placeholder functions to allow for the dimension of the on-site Hilbert space to be
// site dependent. This allows for the implementation of wire networks. Currently, they are just
// returning a fixed dimension.
//---------------------------------------------------------------------------------------------------//

int stateArray::locd(int i) const{
  return dimInfo.locd(i);
}

//---------------------------------------------------------------------------------------------------//

int stateArray::locDMax(int i) const{
  return dimInfo.locDMax(i);
}

void stateArray::createStateArray(int L){
  std::vector<int> locDims(3,0);
  stateArrayAccessStructure.resize(L);
  for(int i=0;i<L;++i){
    locDims[0]=locd(i);
    locDims[1]=locDimR(i);
    locDims[2]=locDimL(i);
    stateArrayAccessStructure[i].generate(locDims);
  }
}
