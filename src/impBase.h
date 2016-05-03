#ifndef BASE_CLASS_FOR_MP_CONSTRUCTS
#define BASE_CLASS_FOR_MP_CONSTRUCTS

#include <complex>
#include "twositeQNOrderMatrix.h"
#include "basisQNOrderMatrix.h"
#include "dimensionTable.h"
#include "pseudoQuantumNumber.h"
#include "baseTensor.h"

class impBase{
 public:
  impBase(){}
 impBase(dimensionTable const &dimInfoBaseIn):dimInfoBase(dimInfoBaseIn){}
  virtual ~impBase(){}
  virtual void subMatrixStart(std::complex<double> *&pStart, int i, int si=0)=0;
  virtual int addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN, std::vector<std::complex<int> > const &source)=0;
  virtual int refineQN(int i, std::vector<std::complex<int> > const& newQN)=0;
  virtual int currentSite()const =0;
  virtual int internalSite()const=0;
  int maxDim()const {return dimInfoBase.D();}
  int length()const {return dimInfoBase.L();}
  dimensionTable const& getDimInfo()const {return dimInfoBase;}
  virtual twositeQNOrderMatrix const& centralIndexTable()const=0;
  virtual basisQNOrderMatrix const& indexTable()const=0;
  virtual pseudoQuantumNumber* getConservedQNs(int i)=0;
  virtual baseTensor<std::complex<double> > const& getSiteTensor(int i)=0;
 protected:
  dimensionTable dimInfoBase;
};

#endif
