#ifndef TMP_CONTAINER
#define TMP_CONTAINER

//---------------------------------------------------------------------------------------------------//
// Basic class for container arrays. It offers a contigous one-dimensional array with a 4D access
// function, i.e. access looks like a 4D array is adressed.
//---------------------------------------------------------------------------------------------------//

template<typename T>
class tmpContainer{
 public:
  tmpContainer(){container=0;}
  tmpContainer(int Lin, int D1in, int D2in, int D3in){initializeContainer(Lin,D1in,D2in,D3in);}
  ~tmpContainer(){delete[] container;}
  T& global_access(int i, int ai1, int ai2, int ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
 protected:
  int L, D1, D2, D3;
  T *container;
  void initializeContainer(int Lin, int D1in, int D2in, int D3in){L=Lin; D1=D1in; D2=D2in; D3=D3in; container=new T[Lin*D1in*D2in*D3in];}
};

#endif
