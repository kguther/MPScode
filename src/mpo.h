#ifndef MPO_OPERATOR
#define MPO_OPERATOR

template<typename T>
class mpo{
 public:
  mpo();
  mpo(int din, int Dwin, int Lin);
  ~mpo();
  T& global_access(int i, int si, int sip, int bi, int bip){return Qoperator[bip+bi*Dw+sip*Dw*Dw+si*Dw*Dw*d+i*Dw*Dw*d*d];}
  void subMatrixStart(int i, T **pStart);
  void initialize(int din, int Dwin, int Lin);
 private:
  int d, Dw;
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
}

template<typename T>
void mpo<T>::initialize(int din, int Dwin, int Lin){
  d=din; 
  Dw=Dwin; 
  Qoperator=new T[Lin*d*Dw*d*Dw];
}

template<typename T>
void mpo<T>::subMatrixStart(int i, T **pStart){
  *pStart=Qoperator+i*Dw*Dw*d*d;
}

#endif
