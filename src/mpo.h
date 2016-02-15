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
  mpo(mpo<T> const &source);
  ~mpo();
  mpo& operator=(mpo<T> const &source);
  T& global_read(int const i, int const si, int const sip, int const bi, int const bip) const{return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  T& global_access(int const i, int const si, int const sip, int const bi, int const bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  int locDimL(int const i) const;
  int locDimR(int const i) const;
  int maxDim() const{return Dw;}
  int maxlocd() const{return d;}
  int length() const{return L;}
  void shift(T const delta);
  void subMatrixStart(T *&pStart, int const i, int const si=0, int const sip=0);
  void initialize(int const din, int const Dwin, int const Lin);
  void setUpSparse();
  void sparseSubMatrixStart(T *&pStart, int const i){pStart=sparseOperator+i*d*d*Dw*Dw;}
  void biSubIndexArrayStart(int *&target, int const i){target=biIndices+i*d*d*Dw*Dw;}
  void bimSubIndexArrayStart(int *&target, int const i){target=bimIndices+i*d*d*Dw*Dw;}
  void siSubIndexArrayStart(int *&target, int const i){target=siIndices+i*d*d*Dw*Dw;}
  void sipSubIndexArrayStart(int *&target, int const i){target=sipIndices+i*d*d*Dw*Dw;}
  int numEls(int const i) const {return nNzero[i];}
  int *nNzero;
 protected:
  void setUpSiteSparse(int const i);
  void mpoCpy(mpo<T> const &source);
  int d, Dw, L;
  T *Qoperator;
  T* sparseOperator;
  int *siIndices, *sipIndices, *biIndices, *bimIndices;
};

template<typename T>
mpo<T>::mpo(){
  Qoperator=0;
  nNzero=0;
  sparseOperator=0;
  biIndices=0;
  bimIndices=0;
  siIndices=0;
  sipIndices=0;
}

template<typename T>
mpo<T>::mpo(int const din, int const Dwin, int const Lin){
  Qoperator=0;
  initialize(din,Dwin,Lin);
}

template<typename T>
mpo<T>::mpo(mpo<T> const  &source){
  mpoCpy(source);
}

template<typename T>
mpo<T>& mpo<T>::operator=(mpo<T> const  &source){
  mpoCpy(source);
  return *this;
}

template<typename T>
void mpo<T>::mpoCpy(mpo<T> const &source){
  Qoperator=0;
  initialize(source.maxlocd(),source.maxDim(),source.length());
  for(int i=0;i<L;++i){
    for(int si=0;si<d;++si){
      for(int sip=0;sip<d;++sip){
	for(int bi=0;bi<locDimR(i);++bi){
	  for(int bim=0;bim<locDimL(i);++bim){
	    global_access(i,si,sip,bi,bim)=source.global_read(i,si,sip,bi,bim);
	  }
	}
      }
    }
  }
  setUpSparse();
}

template<typename T>
mpo<T>::~mpo(){
  delete[] Qoperator;
  delete[] nNzero;
  delete[] sparseOperator;
  delete[] biIndices;
  delete[] bimIndices;
  delete[] siIndices;
  delete[] sipIndices;
}

template<typename T>
void mpo<T>::initialize(int const din, int const Dwin, int const Lin){
  d=din; 
  Dw=Dwin; 
  L=Lin;
  delete[] Qoperator;
  Qoperator=new T[Lin*d*Dw*d*Dw];
  nNzero=0;
  sparseOperator=0;
  biIndices=0;
  bimIndices=0;
  siIndices=0;
  sipIndices=0;
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
  delete[] sparseOperator;
  sparseOperator=0;
  delete[] biIndices;
  biIndices=0;
  delete[] bimIndices;
  bimIndices=0;
  delete[] siIndices;
  siIndices=0;
  delete[] sipIndices;
  sipIndices=0;
  delete[] nNzero;
  nNzero=0;
  nNzero=new int[L];
  sparseOperator=new T[L*Dw*Dw*d*d];
  biIndices=new int[L*Dw*Dw*d*d];
  bimIndices=new int[L*Dw*Dw*d*d];
  siIndices=new int[L*Dw*Dw*d*d];
  sipIndices=new int[L*Dw*Dw*d*d];
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
  pStart=Qoperator+i*Dw*Dw*d*d+sip*Dw*Dw+si*Dw*Dw*d;
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
