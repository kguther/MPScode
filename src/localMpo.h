#ifndef MPO_CLASS_FOR_LOCAL_OPERATORS
#define MPO_CLASS_FOR_LOCAL_OPERATORS

#include "mpo.h"

template<typename T>
class localMpo: public mpo<T>{
 public:
 localMpo():mpo<T>(){}
 localMpo(int const din, int const Dwin, int const Lin, int const initialSite, int *fermionicParities=0, int blockSize=1):mpo<T>(din,Dwin,Lin),i(initialSite),fermionicSign(fermionicParities),size(blockSize){}
  void initializeLocal(int const din, int const Dwin, int const Lin, int const initialSite, int *fermionicParites=0){this->initialize(din,Dwin,Lin);i=initialSite;fermionicSign=fermionicParites;}
  void stepRight();
  int currentSite()const {return i;}
 private: 
  int i, size;
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
	  for(int m=0;m<size;++m){
	    this->global_access(i+1-m,si,sip,bi,bim)=this->global_access(i-m,si,sip,bi,bim);
	  }
	  this->global_access(i+1-size,si,sip,bi,bim)=(si==sip)?fermionicSignFunction(si):0;
	}
      }
    }
  }
  this->setUpSiteSparse(i);
  this->setUpSiteSparse(i+1);
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
