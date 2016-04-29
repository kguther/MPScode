#ifndef MPO_SITE_TENSOR
#define MPO_SITE_TENSOR

#include <vector>

template<typename T>
class mpoSiteTensor{
 public:
  mpoSiteTensor(){}
  //Allocation of Hamiltonian MPO - square matrices are used since no library matrix functions have to be applied - this allows for faster access  
 mpoSiteTensor(int dIn, int DwIn, int DwL, int DwR):d(dIn),Dw(DwIn),lDwL(DwL),lDwR(DwR){siteMatrix.resize(d*d*Dw*Dw);}
  const T& operator()(int si, int sip, int bi, int bim)const {return siteMatrix[bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si];}
  T& operator()(int si, int sip, int bi, int bim){return siteMatrix[bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si];}
  const T& globalAccess(int si, int sip, int bi, int bim)const {return siteMatrix[bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si];}
  T& globalAccess(int si, int sip, int bi, int bim){return siteMatrix[bim+Dw*bi+Dw*Dw*sip+Dw*Dw*d*si];}
  void subMatrixStart(T *&pStart, int si=0, int sip=0){pStart=&(siteMatrix[sip*Dw*Dw+si*Dw*Dw*d]);}
  void subMatrixStart(T const*&pStart, int si=0, int sip=0)const{pStart=&(siteMatrix[sip*Dw*Dw+si*Dw*Dw*d]);}
  void setUpSparse();
  void sparseSubMatrixStart(T const*&pStart)const{pStart=&(sparseOperator[0]);}
  void sparseSubMatrixStart(T *&pStart){pStart=&(sparseOperator[0]);}
  void biSubIndexArrayStart(int const*&target)const {target=&(biIndices[0]);}
  void bimSubIndexArrayStart(int const*&target)const {target=&(bimIndices[0]);}
  void siSubIndexArrayStart(int const*&target)const {target=&(siIndices[0]);}
  void sipSubIndexArrayStart(int const*&target)const {target=&(sipIndices[0]);}

  int numEls()const {return nNzero;}
 private:
  std::vector<int> siIndices, sipIndices, biIndices, bimIndices;
  int nNzero, Dw, d, lDwL, lDwR;
  std::vector<T> sparseOperator, siteMatrix;
};

template<typename T>
void mpoSiteTensor<T>::setUpSparse(){
  sparseOperator.resize(Dw*Dw*d*d);
  biIndices.resize(Dw*Dw*d*d);
  bimIndices.resize(Dw*Dw*d*d);
  siIndices.resize(Dw*Dw*d*d);
  sipIndices.resize(Dw*Dw*d*d);
  nNzero=0;
  int const threshold=1e-10;
  for(int si=0;si<d;++si){
    for(int sip=0;sip<d;++sip){
      for(int bi=0;bi<lDwR;++bi){
	for(int bim=0;bim<lDwL;++bim){
	  if(abs(globalAccess(si,sip,bi,bim))>threshold){
	    sparseOperator[nNzero]=globalAccess(si,sip,bi,bim);
	    biIndices[nNzero]=bi;
	    bimIndices[nNzero]=bim;
	    siIndices[nNzero]=si;
	    sipIndices[nNzero]=sip;
	    ++(nNzero);
	  }
	}
      }
    }
  }
}

#endif
