#ifndef TMP_CONTAINER
#define TMP_CONTAINER

#include<complex>

#define complex_type std::complex<double>

class tmpContainer{
 public:
  tmpContainer(int Lin, int D1in, int D2in, int D3in){L=Lin; D1=D1in; D2=D2in; D3=D3in; container=new complex_type[Lin*D1in*D2in*D3in];}
  ~tmpContainer(){delete[] container;}
  complex_type& access(int i, int ai1, int ai2, int ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
 private:
  int L, D1, D2, D3;
  complex_type *container;
};

#endif
