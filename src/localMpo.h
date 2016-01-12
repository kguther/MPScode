#ifndef MPO_CLASS_FOR_LOCAL_OPERATORS
#define MPO_CLASS_FOR_LOCAL_OPERATORS

#include "mpo.h"

template<typename T>
class localMpo: public mpo<T>{
 public:
 localMpo():mpo(){}
 localMpo(int const din, int const Dwin, int const Lin):mpo(din,Dwin,Lin){}
  void stepRight();
};

template<typename T>
void localMpo<T>::stepRight(){
}

#endif
