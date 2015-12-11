#ifndef ARRAYPROCESSING
#define ARRAYPROCESSING

#include "mkl_complex_defined.h"

using namespace std;

void upperdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout);
void lowerdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout, int ldim=0);
void arraycpy(int const dim1, int const dim2, lapack_complex_double *arraysource, lapack_complex_double *arraytarget);
void transp(int const dim1, int const dim2, lapack_complex_double *array);
void matrixprint(int const n, int const m, lapack_complex_double *array);
void lapackSVD(int dim1, int dim2, lapack_complex_double *U, lapack_complex_double *V, double *diags);

#endif
