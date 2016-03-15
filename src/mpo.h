#ifndef MPO_OPERATOR
#define MPO_OPERATOR

#include <vector>
#include <iostream>

//REMARK: The MPO container consists only of same-sized arrays, even if d is site dependent. In this case, the maximal d has to be used as an input parameter. This creates some unused memory between matrices, but speeds up access to single matrices (which is what we want to do). This works because no matrix operations using external libraries have to be applied to H

template<typename T>
class mpo{
 public:
  mpo();
  mpo(int const din, int const Dwin, int const Lin);
  const T& operator()(int i, int si, int sip, int bi, int bip)const{return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  T& operator()(int i, int si, int sip, int bi, int bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  T* operator()(int i, int si=0, int sip=0){return Qoperator+i*Dw*Dw*d*d+sip*Dw*Dw+si*Dw*Dw*d;}
  const T& global_access(int const i, int const si, int const sip, int const bi, int const bip) const{return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  T& global_access(int const i, int const si, int const sip, int const bi, int const bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  int locDimL(int const i) const;
  int locDimR(int const i) const;
  int maxDim() const{return Dw;}
  int maxlocd() const{return d;}
  int length() const{return L;}
  void shift(T const delta);
  void subMatrixStart(T *&pStart, int const i, int const si=0, int const sip=0);
  void setUpSparse();
  void sparseSubMatrixStart(T *&pStart, int const i){pStart=&(sparseOperator[i*d*d*Dw*Dw]);}
  void biSubIndexArrayStart(int *&target, int const i){target=&(biIndices[i*d*d*Dw*Dw]);}
  void bimSubIndexArrayStart(int *&target, int const i){target=&(bimIndices[i*d*d*Dw*Dw]);}
  void siSubIndexArrayStart(int *&target, int const i){target=&(siIndices[i*d*d*Dw*Dw]);}
  void sipSubIndexArrayStart(int *&target, int const i){target=&(sipIndices[i*d*d*Dw*Dw]);}
  int numEls(int const i) const {return nNzero[i];}
 protected:
  void setUpSiteSparse(int const i);
  int d, Dw, L;
  std::vector<int> nNzero;
  std::vector<T> sparseOperator, Qoperator;
  std::vector<int> siIndices, sipIndices, biIndices, bimIndices;
};

template<typename T>
mpo<T>::mpo(){
}

template<typename T>
mpo<T>::mpo(int const din, int const Dwin, int const Lin):
d(din),
  Dw(Dwin),
  L(Lin)
{
  Qoperator.resize(d*d*Dw*Dw*L);
}

template<typename T>
void mpo<T>::shift(T const delta){
  int lDwL;
  for(int i=0;i<L;++i){
    lDwL=locDimL(i);
    for(int si=0;si<d;++si){
      global_access(i,si,si,0,lDwL)+=delta/2.0;
    }
  }
}

template<typename T>
void mpo<T>::setUpSparse(){
  nNzero.resize(L);
  sparseOperator.resize(L*Dw*Dw*d*d);
  biIndices.resize(L*Dw*Dw*d*d);
  bimIndices.resize(L*Dw*Dw*d*d);
  siIndices.resize(L*Dw*Dw*d*d);
  sipIndices.resize(L*Dw*Dw*d*d);
  for(int i=0;i<L;++i){
    setUpSiteSparse(i);
  }
}

template<typename T>
void mpo<T>::setUpSiteSparse(int const i){
  nNzero[i]=0;
  int const lDwR=locDimR(i);
  int const lDwL=locDimL(i);
  int threshold=1e-10;
  for(int si=0;si<d;++si){
    for(int sip=0;sip<d;++sip){
      for(int bi=0;bi<lDwR;++bi){
	for(int bim=0;bim<lDwL;++bim){
	  if(abs(global_access(i,si,sip,bi,bim))>threshold){
	    sparseOperator[nNzero[i]+i*Dw*Dw*d*d]=global_access(i,si,sip,bi,bim);
	    biIndices[nNzero[i]+i*Dw*Dw*d*d]=bi;
	    bimIndices[nNzero[i]+i*Dw*Dw*d*d]=bim;
	    siIndices[nNzero[i]+i*Dw*Dw*d*d]=si;
	    sipIndices[nNzero[i]+i*Dw*Dw*d*d]=sip;
	    ++(nNzero[i]);
	  }
	}
      }
    }
  }
}

template<typename T>
void mpo<T>::subMatrixStart(T *&pStart, int const i, int const si, int const sip){
  pStart=&(Qoperator[i*Dw*Dw*d*d+sip*Dw*Dw+si*Dw*Dw*d]);
}

template<typename T>
int mpo<T>::locDimL(int const i) const{
  if(i==0){
    return 1;
  }
  return Dw;
}

template<typename T>
int mpo<T>::locDimR(int const i) const{
  if(i==(L-1)){
    return 1;
  }
  return Dw;
}

#endif
