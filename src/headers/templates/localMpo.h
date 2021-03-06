#ifndef MPO_CLASS_FOR_LOCAL_OPERATORS
#define MPO_CLASS_FOR_LOCAL_OPERATORS

#include "templates/mpo.h"

//---------------------------------------------------------------------------------------------------//
// A 'localMpo' is a mpo wich only contains operators acting on a few sites. Usually, it can be seperated into two parts, each acting on spatially different sites. This class allows to move the part operating on sites with higher indices to be moved to the right, a useful feature when measuring correlation functions. Example: particle number per site n_i
//---------------------------------------------------------------------------------------------------//

template<typename T>
class localMpo: public mpo<T>{
 public:
 localMpo():mpo<T>(){}
  //The localMpo has three additional arguments in the constructor. Apart from that, it behaves just like a normal mpo for application purposes. 

  //The first argument is the initial site, from which the 'measurement sweep' is executed (usually 0 or 1), the second argument is the fermionic parity of the local basis states (as C-array) and the last one is the number of sites on which the sweeping part is evaluated (i.e. 1 for n_i and 2 for p-wave superconducting correlations, usually, it is 1)
 localMpo(int const din, int const Dwin, int const Lin, int const initialSite, int *fermionicParities=0, int blockSize=1):mpo<T>(din,Dwin,Lin),i(initialSite),fermionicSign(fermionicParities),size(blockSize){}

//---------------------------------------------------------------------------------------------------//

  void initializeLocal(int const din, int const Dwin, int const Lin, int const initialSite, int *fermionicParites=0){this->initialize(din,Dwin,Lin);i=initialSite;fermionicSign=fermionicParites;}
  void stepRight();
  int currentSite()const {return i;}
  int width()const {return size;}
 private: 
  int i, size;
  int *fermionicSign;
  int fermionicSignFunction(int const si);
};

//---------------------------------------------------------------------------------------------------//

template<typename T>
void localMpo<T>::stepRight(){
  //Special feature of localMpo: moves size MPO-matrices at size i one site to the right (and fills up the vacancy with identity)
  int lDwR=this->locDimR(i);
  int lDwL=this->locDimL(i);
  for(int si=0;si<this->maxlocd();++si){
    for(int sip=0;sip<this->maxlocd();++sip){
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
  for(int m=0;m<=size;++m){
    this->setUpSiteSparse(i+1-m);
  }
  ++i;
}

template<typename T>
int localMpo<T>::fermionicSignFunction(int const si){
  //When moving matrices, one has to take into account for signs occuring due to fermionic statistics if the matrices moved contain unpaired fermionic operators. This is indicated by the fermionicSign flag.
  if(fermionicSign){
    return fermionicSign[si];
  }
  return 1;
}

#endif
