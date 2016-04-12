#ifndef BASE_TENSOR_CLASS
#define BASE_TENSOR_CLASS

#include <vector>

//---------------------------------------------------------------------------------------------------//
// Basic class for a tensor of an MPS (or MPO, although MPOs are stored in another scheme for sake of efficiency). It can have an arbitrary dimension, but this comes at the cost of an ugly direct access. Therefore, when performance is required, getPtr should be used instead.
//---------------------------------------------------------------------------------------------------//

template<typename T>
class baseTensor{
 public:
  baseTensor();
  baseTensor(std::vector<int> const &dims);
  T& operator()(std::vector<int> const &indices);
  const T& operator()(std::vector<int> const &indices) const;
  void generate(std::vector<int> const &dims);
  void getPtr(T *&target, int si=0){target=&(entries[0])+si*factors[0];}
  void getPtr(T const *&target, int si=0)const {target=&(entries[0])+si*factors[0];}
  int setParameterDims(std::vector<int> const &dimsNew);
 private:
  //TODO: replace T* with std::vector<T>
  std::vector<T> entries;
  std::vector<int> dimensions;
  std::vector<int> factors;
  int containerSize;
  void initialize();
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>::baseTensor():containerSize(0){
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>::baseTensor(std::vector<int> const &dims):dimensions(dims){
  initialize();
}

//---------------------------------------------------------------------------------------------------//
// Direct access function for the baseTensor class. Not very fast, so it should be avoided when fast access is required, instead, use getPtr to get C-style access.
//---------------------------------------------------------------------------------------------------//

template<typename T>
const T& baseTensor<T>::operator()(std::vector<int> const &indices) const{
  int position=0;
  int numArgs=(indices.size()>dimensions.size())?dimensions.size():indices.size();
  for(int m=0;m<numArgs;++m){
    position+=factors[m]*indices[m];
  }
  return entries[position];
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
T& baseTensor<T>::operator()(std::vector<int> const &indices){
  return const_cast<T&>(static_cast<baseTensor<T> const&>(*this)(indices));
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void baseTensor<T>::generate(std::vector<int> const &dims){
  dimensions=dims;
  initialize();
}

//---------------------------------------------------------------------------------------------------//
// This feature is used to ramp up the (virtual) bond dimension over the course of calculation.
//---------------------------------------------------------------------------------------------------//

template<typename T>
int baseTensor<T>::setParameterDims(std::vector<int> const &dimsNew){
  int backupCSize=containerSize;
  std::vector<int> backupFactors=factors;
  std::vector<int> backupDims=dimensions;
  std::vector<T> backupEntries=entries;
  dimensions=dimsNew;
  initialize();
  if(backupCSize>containerSize || backupDims.size()!=dimensions.size()){
    containerSize=backupCSize;
    factors=backupFactors;
    entries=backupEntries;
    dimensions=backupDims;
    return 1;
  }
  for(int m=0;m<containerSize;++m){
    entries[m]=0;
  }
  std::vector<int> indices;
  int mReduced, position, backupPosition;
  //This is more sophisticated than I thought. Increasing the size of a Tensor requires to translate the absolute position in the source entries-array to a position in the target entries-array, where both have different factors-arrays.
  for(int m=0;m<backupCSize;++m){
    indices.clear();
    mReduced=m;
    for(int iDim=0;iDim<dimensions.size();++iDim){
      //Get virtual bond indices to translate absolute position in the entries-arrays for different dimensions to comparable variables
      indices.push_back(mReduced/backupFactors[iDim]);
      mReduced=mReduced%(backupFactors[iDim]);
    }
    backupPosition=0;
    position=0;
    for(int jDim=0;jDim<dimensions.size();++jDim){
      position+=factors[jDim]*indices[jDim];
      backupPosition+=backupFactors[jDim]*indices[jDim];
    }
    entries[position]=backupEntries[backupPosition];
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void baseTensor<T>::initialize(){
  containerSize=1;
  factors.resize(dimensions.size());
  for(int m=0;m<dimensions.size();++m){
    containerSize*=dimensions[m];
    factors[m]=1;
    for(int k=(m+1);k<dimensions.size();++k){
      factors[m]*=dimensions[k];
    }
  }
  //deleting the memory has to be called manually because there is a case where initialization and deallocating are reversed (using a backup pointer)
  entries.resize(containerSize);
  for(int m=0;m<containerSize;++m){
    entries[m]=0;
  }
}

#endif
