#ifndef BASE_CLASS_FOR_MP_CONSTRUCTS
#define BASE_CLASS_FOR_MP_CONSTRUCTS

#include <complex>
#include "mkl_complex_defined.h"
#include "twositeQNOrderMatrix.h"

class impBase{
 public:
  virtual ~impBase(){}
  virtual void subMatrixStart(lapack_complex_double *&pStart, int i, int si=0)=0;
  virtual void addSite(int Lnew, int i, std::vector<std::complex<int> > const &targetQN)=0;
  virtual int refineQN(int i, std::vector<std::complex<int> > const& newQN)=0;
  virtual int currentSite()const =0;
  virtual twositeQNOrderMatrix const& centralIndexTable()=0;
};

#endif
