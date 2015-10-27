#ifndef MPO_OPERATOR
#define MPO_OPERATOR

template<typename T>
class mpo{
 public:
  mpo();
  mpo(int din, int Dwin, int Lin);
  ~mpo();
  T& global_access(int i, int si, int sip, int bi, int bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  T& sparse_access(int si, int sip, int bi, int n){return Qoperator[bimIndex[n+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d]+bi*Dw+si*Dw*Dw+sip*Dw*Dw*d];}
  int sparse_nNzero(int si, int sip, int bi){return nNzero[bi+sip*Dw+si*Dw*d];}
  void subMatrixStart(int i, T **pStart);
  void initialize(int din, int Dwin, int Lin);
  void setUpSiteSparse(int i);
 private:
  int d, Dw;
  int *nNzero, *bimIndex;
  T *Qoperator;
};

template<typename T>
mpo<T>::mpo(){
  Qoperator=0;
}

template<typename T>
mpo<T>::mpo(int din, int Dwin, int Lin){
  initialize(din,Dwin,Lin);
}

template<typename T>
mpo<T>::~mpo(){
  delete[] Qoperator;
  delete[] nNzero;
  delete[] bimIndex;
}

template<typename T>
void mpo<T>::initialize(int din, int Dwin, int Lin){
  d=din; 
  Dw=Dwin; 
  Qoperator=new T[Lin*d*Dw*d*Dw];
  nNzero=new int[Dw*d*d];
  bimIndex=new int[Dw*Dw*d*d];
}

template<typename T>
void mpo<T>::setUpSiteSparse(int i){
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
void mpo<T>::subMatrixStart(int i, T **pStart){
  *pStart=Qoperator+i*Dw*Dw*d*d;
}

#endif
