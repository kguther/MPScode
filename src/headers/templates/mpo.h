#ifndef MPO_OPERATOR
#define MPO_OPERATOR

#include <vector>
#include "templates/mpoSiteTensor.h"

//---------------------------------------------------------------------------------------------------//
// MPO class for handling of sparse mpos.
//---------------------------------------------------------------------------------------------------//
//REMARK: The MPO container consists only of same-sized arrays, even if d is site dependent. In this case, the maximal d has to be used as an input parameter. This creates some unused memory between matrices, but speeds up access to single matrices (which is what we want to do). This works because no matrix operations using external libraries have to be applied to H

template<typename T>
class mpo{
 public:
  mpo();
  //Usually, the MPO is created supplying the (maximum) local hilbert space dimension din, the MPO Bond dimension DWin and the system size Lin
  mpo(int din, int Dwin, int Lin);
  //It may then be filled with entries manually, the sparse form is then generated internally
  const T& operator()(int i, int si, int sip, int bi, int bip)const{return Qoperator[i](si,sip,bi,bip);}
  T& operator()(int i, int si, int sip, int bi, int bip){return Qoperator[i](si,sip,bi,bip);}
  T* operator()(int i, int si=0, int sip=0){
    T* pStart;
    Qoperator[i].subMatrixStart(pStart,si,sip);
    return pStart;
  }
  const T& global_access(int i, int si, int sip, int bi, int bip) const{return Qoperator[i].globalAccess(si,sip,bi,bip);}
  T& global_access(int i, int si, int sip, int bi, int bip){return Qoperator[i].globalAccess(si,sip,bi,bip);}
  int locDimL(int i) const;
  int locDimR(int i) const;
  int maxDim() const{return Dw;}
  int maxlocd() const{return d;}
  int length() const{return L;}
  void shift(T const delta);
  int loadTI(mpo<T> const &source);

  //Submatrix access functions
  void subMatrixStart(T *&pStart, int i, int si=0, int sip=0){Qoperator[i].subMatrixStart(pStart,si,sip);}
  void subMatrixStart(T const*&pStart, int i, int si=0, int sip=0)const {Qoperator[i].subMatrixStart(pStart,si,sip);}
  
  //Sparse structure access functions
  void sparseSubMatrixStart(T const*&pStart, int i)const {Qoperator[i].sparseSubMatrixStart(pStart);}
  void biSubIndexArrayStart(int const*&target, int i)const {Qoperator[i].biSubIndexArrayStart(target);}
  void bimSubIndexArrayStart(int const*&target, int i)const {Qoperator[i].bimSubIndexArrayStart(target);}
  void siSubIndexArrayStart(int const*&target, int i)const {Qoperator[i].siSubIndexArrayStart(target);}
  void sipSubIndexArrayStart(int const*&target, int i)const {Qoperator[i].sipSubIndexArrayStart(target);}

  void setUpSparse();
  int numEls(int const i) const {return Qoperator[i].numEls();}
  mpoSiteTensor<T>const & getSiteTensor(int i){return Qoperator[i];}
  
 protected:
  void setUpSiteSparse(int const i);
  int d, Dw, L;
  std::vector<mpoSiteTensor<T> > Qoperator;
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
  Qoperator.resize(L);
  for(int i=0;i<L;++i){
    Qoperator[i]=mpoSiteTensor<T>(d,Dw,locDimL(i),locDimR(i));
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
int mpo<T>::loadTI(mpo<T> const &source){
  //Load a mpo by filling it with a sitematrix from another mpo (and treaing the ends seperately).
  if(source.d!=d || source.Dw!=Dw){
    //Only makes sense when target and source have the same dimensions
    return -1;
  }
  for(int i=1;i<L-1;++i){
    Qoperator[i]=source.Qoperator[1];
  }
  Qoperator.back()=source.Qoperator.back();
  Qoperator.front()=source.Qoperator.front();
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
  for(int i=0;i<L;++i){
    setUpSiteSparse(i);
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void mpo<T>::setUpSiteSparse(int const i){
  //Extract the nonzero entries at site i
  Qoperator[i].setUpSparse();
}  

//---------------------------------------------------------------------------------------------------//
// These are the local bond dimensions of the mpo. While comparable to those of an mps in principle,
// there is no need for normalizability of the mpo. Therefore, the local dimension is constant except
// for the first and the last sites.
//---------------------------------------------------------------------------------------------------//


template<typename T>
int mpo<T>::locDimL(int i) const{
  if(i==0){
    return 1;
  }
  return Dw;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
int mpo<T>::locDimR(int i) const{
  if(i==(L-1)){
    return 1;
  }
  return Dw;
}

#endif
