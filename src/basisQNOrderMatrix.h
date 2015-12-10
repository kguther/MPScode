#ifndef MATRIX_TO_CREATE_QN_BLOCK_ORDERING
#define MATRIX_TO_CREATE_QN_BLOCK_ORDERING

#include "mkl_complex_defined.h"
#include "quantumNumber.h"

class basisQNOrderMatrix{
 public:
  basisQNOrderMatrix(int const i, int const dimin, quantumNumber *conservedQNsin);
  lapack_complex_double *matrix;
 private:
  quantumNumber *conservedQNs;
  int dim;
};

#endif
