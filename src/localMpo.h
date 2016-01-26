#ifndef MPO_CLASS_FOR_LOCAL_OPERATORS
#define MPO_CLASS_FOR_LOCAL_OPERATORS

#include "mpo.h"

template<typename T>
class localMpo: public mpo<T>{
 public:
 localMpo():mpo(){}
 localMpo(int const din, int const Dwin, int const Lin, int const initialSite, int *fermionicParities=0):mpo<T>(din,Dwin,Lin),i(initialSite),fermionicSign(fermionicParities){}
  void stepRight();
  int currentSite()const {return i;}
 private: 
  int i;
  int *fermionicSign;
  int fermionicSignFunction(int const si);
};

template<typename T>
void localMpo<T>::stepRight(){
  int lDwR=this->locDimR(i);
  int lDwL=this->locDimL(i);
  for(int si=0;si<this->d;++si){
    for(int sip=0;sip<this->d;++sip){
      for(int bi=0;bi<lDwR;++bi){
	for(int bim=0;bim<lDwL;++bim){
	  this->global_access(i+1,si,sip,bi,bim)=this->global_access(i,si,sip,bi,bim);
	  this->global_access(i,si,sip,bi,bim)=(si==sip)?fermionicSignFunction(si):0;
	}
      }
    }
  }
  ++i;
}

template<typename T>
int localMpo<T>::fermionicSignFunction(int const si){
  if(fermionicSign){
    return fermionicSign[si];
  }
  return 1;
}

#endif
