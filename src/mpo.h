#ifndef MPO_OPERATOR
#define MPO_OPERATOR

#include <vector>
#include <iostream>

//---------------------------------------------------------------------------------------------------//
// MPO class for handling of sparse mpos.
//---------------------------------------------------------------------------------------------------//
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
  int loadTI(mpo<T> const &source);
  void subMatrixStart(T *&pStart, int const i, int const si=0, int const sip=0);
  void setUpSparse();
  void sparseSubMatrixStart(T *&pStart, int const i){pStart=&(sparseOperator[i*d*d*Dw*Dw]);}
  void biSubIndexArrayStart(int *&target, int const i){target=&(biIndices[i*d*d*Dw*Dw]);}
  void bimSubIndexArrayStart(int *&target, int const i){target=&(bimIndices[i*d*d*Dw*Dw]);}
  void siSubIndexArrayStart(int *&target, int const i){target=&(siIndices[i*d*d*Dw*Dw]);}
  void sipSubIndexArrayStart(int *&target, int const i){target=&(sipIndices[i*d*d*Dw*Dw]);}
  int numEls(int const i) const {return nNzero[i];}
  int setParameterL(int Lnew);
 protected:
  void setUpSiteSparse(int const i);
  int d, Dw, L;
  std::vector<int> nNzero;
  std::vector<T> sparseOperator, Qoperator;
  std::vector<int> siIndices, sipIndices, biIndices, bimIndices;
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
mpo<T>::mpo(){
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
mpo<T>::mpo(int const din, int const Dwin, int const Lin):
d(din),
  Dw(Dwin),
  L(Lin)
{
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access
  Qoperator.resize(d*d*Dw*Dw*L);
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
int mpo<T>::loadTI(mpo<T> const &source){
  //Load a mpo by filling it with a sitematrix from another mpo (and treaing the ends seperately).
  if(source.d!=d || source.Dw!=Dw){
    //Only makes sense when target and source have the same dimensions
    return -1;
  }
  int const siteDim=d*d*Dw*Dw;
  int const sourceEndOffset=(source.L-1)*siteDim;
  for(int m=0;m<siteDim;++m){
    Qoperator[m]=source.Qoperator[m];
    Qoperator[m+(L-1)*siteDim]=source.Qoperator[m+sourceEndOffset];
  }
  for(int i=1;i<L-1;++i){
    for(int m=0;m<siteDim;++m){
      Qoperator[m+i*siteDim]=source.Qoperator[m+siteDim];
    }
  }
  setUpSparse();
  return 0;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void mpo<T>::shift(T const delta){
  //Deprecated
  int lDwL;
  for(int i=0;i<L;++i){
    lDwL=locDimL(i);
    for(int si=0;si<d;++si){
      global_access(i,si,si,0,lDwL)+=delta/2.0;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void mpo<T>::setUpSparse(){
  //Generates the index tables for the treatment of the mpo as sparse matrix
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

//---------------------------------------------------------------------------------------------------//

template<typename T>
void mpo<T>::setUpSiteSparse(int const i){
  //Extract the nonzero entries at site i
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

//---------------------------------------------------------------------------------------------------//

template<typename T>
int mpo<T>::setParameterL(int Lnew){
  //Deprecated
  if(Lnew==L){
    return 0;
  }
  if(Lnew<L){
    return -1;
  }
  std::vector<T> buffer;
  int const deltaL=L-Lnew;
  buffer.resize(deltaL*d*d*Dw*Dw);
  for(int iInsert=0;iInsert<deltaL;++iInsert){
    for(int m=0;m<d*d*Dw*Dw;++m){
      buffer[m+iInsert*d*d*Dw*Dw]=Qoperator[m+d*d*Dw*Dw*L/2];
    }
  }
  Qoperator.insert(Qoperator.begin()+L/2+deltaL/2,buffer.begin(),buffer.end());
  L=Lnew;
}
    

//---------------------------------------------------------------------------------------------------//

template<typename T>
void mpo<T>::subMatrixStart(T *&pStart, int const i, int const si, int const sip){
  //Allows for fast access by setting a pointer directly to some subarray start
  pStart=&(Qoperator[i*Dw*Dw*d*d+sip*Dw*Dw+si*Dw*Dw*d]);
}

//---------------------------------------------------------------------------------------------------//
// These are the local bond dimensions of the mpo. While comparable to those of an mps in principle,
// there is no need for normalizability of the mpo. Therefore, the local dimension is constant except
// for the first and the last sites.
//---------------------------------------------------------------------------------------------------//


template<typename T>
int mpo<T>::locDimL(int const i) const{
  if(i==0){
    return 1;
  }
  return Dw;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
int mpo<T>::locDimR(int const i) const{
  if(i==(L-1)){
    return 1;
  }
  return Dw;
}

#endif
