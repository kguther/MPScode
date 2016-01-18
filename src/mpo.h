#ifndef MPO_OPERATOR
#define MPO_OPERATOR


//REMARK: The MPO container consists only of same-sized arrays, even if d is site dependent. In this case, the maximal d has to be used as an input parameter. This creates some unused memory between matrices, but speeds up access to single matrices (which is what we want to do). This works because no matrix operations using external libraries have to be applied to H

template<typename T>
class mpo{
 public:
  mpo();
  mpo(int const din, int const Dwin, int const Lin);
  ~mpo();
  T& global_access(int const i, int const si, int const sip, int const bi, int const bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  int locDimL(int const i);
  int locDimR(int const i);
  int maxDim() const{return Dw;}
  int length() const{return L;}
  void shift(T const delta);
  T& sparse_access(int const si, int const sip, int const bi, int const n){return Qoperator[bimIndex[n+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d]+bi*Dw+si*Dw*Dw+sip*Dw*Dw*d];}
  int sparse_nNzero(int const si, int const sip, int const bi){return nNzero[bi+sip*Dw+si*Dw*d];}
  void subMatrixStart(T *&pStart, int const i, int const si=0, int const sip=0);
  void initialize(int const din, int const Dwin, int const Lin);
  void setUpSiteSparse(int const i);
 private:
  int d, Dw, L;
  int *nNzero, *bimIndex;
  T *Qoperator;
};

template<typename T>
mpo<T>::mpo(){
  Qoperator=0;
}

template<typename T>
mpo<T>::mpo(int const din, int const Dwin, int const Lin){
  initialize(din,Dwin,Lin);
}

template<typename T>
mpo<T>::~mpo(){
  delete[] Qoperator;
  //delete[] nNzero;
  //delete[] bimIndex;
}

template<typename T>
void mpo<T>::initialize(int const din, int const Dwin, int const Lin){
  d=din; 
  Dw=Dwin; 
  L=Lin;
  Qoperator=new T[Lin*d*Dw*d*Dw];
  //nNzero=new int[Dw*d*d];
  //bimIndex=new int[Dw*Dw*d*d];
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
void mpo<T>::setUpSiteSparse(int const i){
  int threshold=1e-10;
  for(int si=0;si<d;si++){
    for(int sip=0;sip<d;sip++){
      for(int bi=0;bi<Dw;bi++){
	nNzero[bi+sip*Dw+si*Dw*d]=0;
	for(int bim=0;bim<Dw;bim++){
	  if(abs(global_access(i,si,sip,bi,bim))>threshold){
	    bimIndex[nNzero[bi+sip*Dw+si+Dw*d]+bi*Dw+sip*Dw*Dw+si*Dw*d*Dw]=bim;
	    (nNzero[bi+sip*Dw+si*Dw*d])++;
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
int mpo<T>::locDimL(int const i){
  if(i==0){
    return 1;
  }
  return Dw;
}

template<typename T>
int mpo<T>::locDimR(int const i){
  if(i==(L-1)){
    return 1;
  }
  return Dw;
}

#endif
