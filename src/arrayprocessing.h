#ifndef ARRAYPROCESSING
#define ARRAYPROCESSING

#include <lapacke.h>

using namespace std;

void upperdiag(int dim1, int dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout); //for lapack postprocessing
void lowerdiag(int dim1, int dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout);
void arraycpy(int dim1, int dim2, lapack_complex_double *arraysource, lapack_complex_double *arraytarget);
void transp(int dim1, int dim2, lapack_complex_double *array);
void matrixprint(int n, int m, lapack_complex_double *array);

#endif
