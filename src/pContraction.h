#ifndef PARTIAL_CONTRACTION
#define PARTIAL_CONTRACTION

#include "tmpContainer.h"

//---------------------------------------------------------------------------------------------------//
// Class for storage of partial contractions of the network. Basically, this is a 4D container with
// the additional functionality of access to the partial contractions up to a specific site. 
//---------------------------------------------------------------------------------------------------//

template<typename T>
class pContraction: public tmpContainer<T>{
 public:
  pContraction();
  pContraction(int Lin, int Din, int Dwin);
  void subContractionStart(int i, T **pStart);
  void initialize(int Lin, int Din, int Dwin);
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
pContraction<T>::pContraction(){
  this->container=0;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
pContraction<T>::pContraction(int Lin, int Din, int Dwin){
  tmpContainer(Lin, Din, Dwin, Din);
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void pContraction<T>::subContractionStart(int i, T **pStart){
  *pStart=this->container+i*this->D1*this->D3*this->D2;
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
void pContraction<T>::initialize(int Lin, int Din, int Dwin){
  //Initialization handy for updating D
  delete[] this->container;
  this->initializeContainer(Lin,Din,Dwin,Din);
}

#endif
