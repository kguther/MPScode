#ifndef BASE_TENSOR_CLASS
#define BASE_TENSOR_CLASS

#include <vector>

template<typename T>
class baseTensor{
 public:
  baseTensor();
  baseTensor(std::vector<int> const &dims);
  baseTensor(baseTensor<T> const &source);
  ~baseTensor();
  baseTensor& operator=(baseTensor<T> const &source);
  T& operator()(std::vector<int> const &indices);
  const T& operator()(std::vector<int> const &indices) const;
  void generate(std::vector<int> const &dims);
  void getPtr(T *&target, int si=0){target=entries+si*factors[0];}
  void getPtr(T const *&target, int si=0)const {target=entries+si*factors[0];}
  int setParameterDims(std::vector<int> const &dimsNew);
 private:
  T *entries;
  std::vector<int> dimensions;
  std::vector<int> factors;
  int containerSize;
  void tensorCpy(baseTensor<T> const &source);
  void initialize();
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>::baseTensor():entries(0),containerSize(0){
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>::baseTensor(std::vector<int> const &dims):dimensions(dims){
  initialize();
}

//---------------------------------------------------------------------------------------------------//  

template<typename T>
baseTensor<T>::baseTensor(baseTensor<T> const &source):entries(0),containerSize(0){
  tensorCpy(source);
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>::~baseTensor(){
  delete[] entries;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
baseTensor<T>& baseTensor<T>::operator=(baseTensor<T> const &source){
  tensorCpy(source);
  return *this;
}

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

template<typename T>
T& baseTensor<T>::operator()(std::vector<int> const &indices){
  return const_cast<T&>(static_cast<baseTensor<T> const&>(*this)(indices));
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void baseTensor<T>::tensorCpy(baseTensor<T> const &source){
  if(&source!=this && source.entries){
    containerSize=source.containerSize;
    dimensions=source.dimensions;
    factors=source.factors;
    delete[] entries;
    entries=new T[containerSize];
    for(int m=0;m<containerSize;++m){
      entries[m]=source.entries[m];
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void baseTensor<T>::generate(std::vector<int> const &dims){
  dimensions=dims;
  delete[] entries;
  initialize();
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
int baseTensor<T>::setParameterDims(std::vector<int> const &dimsNew){
  int backupCSize=containerSize;
  std::vector<int> backupFactors=factors;
  std::vector<int> backupDims=dimensions;
  T *backupEntries=entries;
  dimensions=dimsNew;
  initialize();
  if(backupCSize>containerSize || backupDims.size()!=dimensions.size()){
    containerSize=backupCSize;
    factors=backupFactors;
    delete[] entries;
    entries=backupEntries;
    dimensions=backupDims;
    return 1;
  }
  for(int m=0;m<containerSize;++m){
    entries[m]=0;
  }
  std::vector<int> indices;
  int mReduced, position, backupPosition;
  for(int m=0;m<backupCSize;++m){
    indices.clear();
    mReduced=m;
    for(int iDim=0;iDim<dimensions.size();++iDim){
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
  delete[] backupEntries;
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
  entries=new T[containerSize];
  for(int m=0;m<containerSize;++m){
    entries[m]=0;
  }
}

#endif
