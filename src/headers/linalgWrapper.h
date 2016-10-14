#ifndef WRAPPER_FUNCTIONS_FOR_LAPACK_AND_BLAS
#define WRAPPER_FUNCTIONS_FOR_LAPACK_AND_BLAS

#include "mkl_complex_defined.h"
#include "mpstype.h"

void cblas_gemm(CBLAS_ORDER const order, CBLAS_TRANSPOSE const transA, CBLAS_TRANSPOSE const transB, int const m, int const n, int const K, void const *alpha, void const *A, int const lda, void const *B, int const ldb, void const *beta, void *C, int const ldc);
double cblas_norm(int const n, void const *x, int const incx);
void cblas_trmm(CBLAS_ORDER const order, CBLAS_SIDE const side, CBLAS_UPLO const uplo, CBLAS_TRANSPOSE const transA, CBLAS_DIAG const diag, int const m, int const n, void const *alpha, void const *A, int const lda, void *B, int const ldb);
void cblas_scal(int const n, void const *alpha, void *X, int const incx);
void cblas_axpy(int const n, void const *alpha, void const *x, int const incx, void *y, int const incy);
void lapacke_svd(char job, int blockDimL, int blockDimR, void *M, double *diags, void *U, void *VT, int i);

#endif
