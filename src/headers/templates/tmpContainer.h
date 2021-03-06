#ifndef TMP_CONTAINER
#define TMP_CONTAINER

//---------------------------------------------------------------------------------------------------//
// Basic class for container arrays. It offers a contigous one-dimensional array with a 4D access
// function, i.e. access looks like a 4D array is adressed.
//---------------------------------------------------------------------------------------------------//

#include <vector>

template<typename T>
class tmpContainer{
 public:
 tmpContainer():L(0),D1(0),D2(0),D3(0){}
 tmpContainer(int Lin, int D1in, int D2in, int D3in):L(Lin),D1(D1in),D2(D2in),D3(D3in){container.resize(L*D1*D2*D3);}
  const T& operator()(int i, int ai1, int ai2, int ai3)const {return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
  T& operator()(int i, int ai1, int ai2, int ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
  T& global_access(int const i, int const ai1, int const ai2, int const ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
  void getPtr(T *&target){target=&(container[0]);}
  void initializeContainer(int Lin, int D1in, int D2in, int D3in){L=Lin; D1=D1in; D2=D2in; D3=D3in; container.resize(L*D1*D2*D3);}
 protected:
  int L, D1, D2, D3;
  std::vector<T> container;
 private:
  //Copying tmpContainers in a meaningful way can easily lead to a memory overflow, therefore, it is forbidden.
  tmpContainer(tmpContainer<T> const &doNot);
  tmpContainer& operator=(tmpContainer<T> const &doNot);
};

//---------------------------------------------------------------------------------------------------//
// Some more obscure containers which might be useful someday.
//---------------------------------------------------------------------------------------------------//

template<typename T>
class dynamic5DContainer{
 public:
  dynamic5DContainer(){container=0;}
  ~dynamic5DContainer(){delete[] container;}
  void generate(int d1, int d2, int d3, int d4, int d5){delete[] container; D1=d1; D2=d2; D3=d3; D4=d4; D5=d5; container=new T[d1*d2*d3*d4*d5];}
  T& global_access(int const a1, int const a2, int const a3, int const a4, int const a5){return container[a5+a4*D5+a3*D4*D5+a2*D3*D4*D5+a1*D2*D3*D4*D5];}
  void getPtr(T *&target){target=container;}
 private:
  T *container;
  int D1,D2,D3,D4,D5;
};

#endif
