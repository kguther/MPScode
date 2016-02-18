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
  const T& operator()(int i, int ai1, int ai2, int ai3)const {return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
  T& operator()(int i, int ai1, int ai2, int ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
  T& global_access(int const i, int const ai1, int const ai2, int const ai3){return container[ai3+ai2*D3+ai1*D3*D2+i*D1*D2*D3];}
 protected:
  int L, D1, D2, D3;
  T *container;
  void initializeContainer(int Lin, int D1in, int D2in, int D3in){L=Lin; D1=D1in; D2=D2in; D3=D3in; container=new T[Lin*D1in*D2in*D3in];}
 private:
  //Copying tmpContainers in a meaningful way can easily lead to a memory overflow, therefore, it is forbidden.
  tmpContainer(tmpContainer<T> const &doNot);
  tmpContainer& operator=(tmpContainer<T> const &doNot);
};

//---------------------------------------------------------------------------------------------------//
// Some more obscure containers which might be useful someday.
//---------------------------------------------------------------------------------------------------//

template<typename T>
class dynamicContainer: public tmpContainer<T>{
 public: 
  dynamicContainer(){tmpContainer<T>();}
  void generate(int d1, int d2, int d3, int d4){delete[] this->container; this->initializeContainer(d1,d2,d3,d4);}
};

template<typename T>
class dynamic5DContainer{
 public:
  dynamic5DContainer(){container=0;}
  ~dynamic5DContainer(){delete[] container;}
  void generate(int d1, int d2, int d3, int d4, int d5){delete[] container; D1=d1; D2=d2; D3=d3; D4=d4; D5=d5; container=new T[d1*d2*d3*d4*d5];}
  T& global_access(int const a1, int const a2, int const a3, int const a4, int const a5){return container[a5+a4*D5+a3*D4*D5+a2*D3*D4*D5+a1*D2*D3*D4*D5];}
 private:
  T *container;
  int D1,D2,D3,D4,D5;
};

#endif
