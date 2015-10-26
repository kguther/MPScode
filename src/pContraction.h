#ifndef PARTIAL_CONTRACTION
#define PARTIAL_CONTRACTION

template<typename T>
class pContraction{
 public:
  pContraction();
  pContraction(int Lin, int Din, int Dwin);
  ~pContraction();
  T& global_access(int i, int ai, int bi, int aip){return pctr[aip+bi*D+ai*D*Dw+i*D*Dw*D];}
  void subContractionStart(int i, T **pStart);
  void initialize(int Lin, int Din, int Dwin);
 private:
  int D, Dw, L;
  T *pctr;
};

template<typename T>
pContraction<T>::pContraction(){
  pctr=0;
}

template<typename T>
pContraction<T>::pContraction(int Lin, int Din, int Dwin){
  initialize(Lin, Din, Dwin);
}

template<typename T>
pContraction<T>::~pContraction(){
  delete[] pctr;
}

template<typename T>
void pContraction<T>::subContractionStart(int i, T **pStart){
  *pStart=pctr+i*D*D*Dw;
}

template<typename T>
void pContraction<T>::initialize(int Lin, int Din, int Dwin){
  delete[] pctr;
  L=Lin;
  D=Din;
  Dw=Dwin;
  pctr=new T[L*D*Dw*D];
}

#endif
