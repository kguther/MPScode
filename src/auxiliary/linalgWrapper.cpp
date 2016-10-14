#include "linalgWrapper.h"
#include "exceptionClasses.h"
#include <memory>
#include <iostream>

void cblas_gemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, int const m, int const n, int const K, void const *alpha, void const *A, int const lda, void const *B, int const ldb, void const *beta, void *C, int const ldc){
#ifdef REAL_MPS_ENTRIES
  cblas_dgemm(order,transA,transB,m,n,K,*(static_cast<double const *>(alpha)),static_cast<double const *>(A),lda,static_cast<double const *>(B),ldb,*(static_cast<double const *>(beta)),static_cast<double *>(C),ldc);
#else
  cblas_zgemm(order,transA,transB,m,n,K,alpha,A,lda,B,ldb,beta,C,ldc);
#endif
}

double cblas_norm(int const n, void const *x, int const incx){
#ifdef REAL_MPS_ENTRIES
  return cblas_dnrm2(n,static_cast<double const *>(x),incx);
#else
  return cblas_dznrm2(n,x,incx);
#endif
}

void cblas_trmm(CBLAS_ORDER const order, CBLAS_SIDE const side, CBLAS_UPLO const uplo, CBLAS_TRANSPOSE const transA, CBLAS_DIAG const diag, int const m, int const n, void const *alpha, void const *A, int const lda, void *B, int const ldb){
#ifdef REAL_MPS_ENTRIES
  cblas_dtrmm(order,side,uplo,transA,diag,m,n,*(static_cast<double const *>(alpha)),static_cast<double const *>(A),lda,static_cast<double *>(B),ldb);
#else
  cblas_ztrmm(order,side,uplo,transA,diag,m,n,alpha,A,lda,B,ldb);
#endif  
}

void cblas_scal(int const n, void const *alpha, void *X, int const incx){
#ifdef REAL_MPS_ENTRIES
  cblas_dscal(n,*(static_cast<double const *>(alpha)),static_cast<double *>(X),incx);
#else
  cblas_zscal(n,alpha,X,incx);
#endif
}

void cblas_axpy(int const n, void const *alpha, void const *x, int const incx, void *y, int const incy){
#ifdef REAL_MPS_ENTRIES
  cblas_daxpy(n,*(static_cast<double const *>(alpha)),static_cast<double const *>(x),incx,static_cast<double *>(y),incy);
#else
  cblas_zaxpy(n,alpha,x,incx,y,incy);
#endif
}

void lapacke_svd(char job, int blockDimL, int blockDimR, void *M, double *diags, void *U, void *VT, int i){
  int const containerDim=(blockDimL>blockDimR)?blockDimR:blockDimL;
#ifdef USE_MKL
#ifdef REAL_MPS_ENTRIES
  int info=LAPACKE_dgesdd(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<double*>(M),blockDimL,diags,static_cast<double*>(U),blockDimL,static_cast<double*>(VT),blockDimR);
#else
  int info=LAPACKE_zgesdd(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<std::complex<double>*>(M),blockDimL,diags,static_cast<std::complex<double>*>(U),blockDimL,static_cast<std::complex<double>*>(VT),blockDimR);
#endif
  if(info>0){
    std::unique_ptr<double[]> workBufP(new double[containerDim]);
    double *workBuf=workBufP.get();
#ifdef REAL_MPS_ENTRIES
    info=LAPACKE_dgesvd(LAPACK_COL_MAJOR,job,job,blockDimL,blockDimR,static_cast<double*>(M),blockDimL,diags,static_cast<double*>(U),blockDimL,static_cast<double*>(VT),blockDimR,workBuf);
#else
    info=LAPACKE_zgesvd(LAPACK_COL_MAJOR,job,job,blockDimL,blockDimR,static_cast<std::complex<double>*>(M),blockDimL,diags,static_cast<std::complex<double>*>(U),blockDimL,static_cast<std::complex<double>*>(VT),blockDimR,workBuf);
#endif
  }
#endif
  //There seems to be a bug in liblapacke providing a wrong size for the work array, which can lead to a segfault. This bug is not present in the mkl implementation
  //Maybe add a workspace query to get rid of the explicit specification, which is rather bug-prone
#ifndef USE_MKL
  int const maxDim=(blockDimR>blockDimL)?blockDimR:blockDimL;
  int lwork=-1;
  int const lrwork=(5*containerDim*(1+containerDim)>(containerDim*(1+2*(maxDim+containerDim))))?5*containerDim*(1+containerDim):containerDim*(1+2*(maxDim+containerDim));
  std::unique_ptr<int[]> iworkP(new int[8*containerDim]);
  std::unique_ptr<double[]> rworkP(new double[lrwork]);
  std::unique_ptr<mpsEntryType[]> workP(new mpsEntryType[1]);
  mpsEntryType *work=workP.get();
  double *rwork=rworkP.get();
  int *iwork=iworkP.get();
#ifdef REAL_MPS_ENTRIES
  int info=LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<double*>(M),blockDimL,diags,static_cast<double*>(U),blockDimL,static_cast<double*>(VT),blockDimR,work,lwork,iwork);
  lwork=static_cast<int>(work[0]);
#else
  int info=LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<std::complex<double>*>(M),blockDimL,diags,static_cast<std::complex<double>*>(U),blockDimL,static_cast<std::complex<double>*>(VT),blockDimR,work,lwork,rwork,iwork);
  lwork=(containerDim*(containerDim+2)+maxDim);
#endif
  workP.reset(new mpsEntryType[lwork]);
  work=workP.get();
#ifdef REAL_MPS_ENTRIES
  info=LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<double*>(M),blockDimL,diags,static_cast<double*>(U),blockDimL,static_cast<double*>(VT),blockDimR,work,lwork,iwork);
#else
  info=LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,job,blockDimL,blockDimR,static_cast<std::complex<double>*>(M),blockDimL,diags,static_cast<std::complex<double>*>(U),blockDimL,static_cast<std::complex<double>*>(VT),blockDimR,work,lwork,rwork,iwork);
#endif
  if(info)
    std::cout<<"SVD FAILURE/n";
  if(info>0){
    //In case of failure of the divide and conquer algorithm, use the zgesvd routine
    //Failure with info<0 is due to wrong input and cannot be remedied by using another algorithm
#ifdef REAL_MPS_ENTRIES
    info=LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,job,job,blockDimL,blockDimR,static_cast<double*>(M),blockDimL,diags,static_cast<double*>(U),blockDimL,static_cast<double*>(VT),blockDimR,work,lwork);
#else
    info=LAPACKE_zgesvd_work(LAPACK_COL_MAJOR,job,job,blockDimL,blockDimR,static_cast<std::complex<double>*>(M),blockDimL,diags,static_cast<std::complex<double>*>(U),blockDimL,static_cast<std::complex<double>*>(VT),blockDimR,work,lwork,rwork);
#endif
  }
#endif
  if(info){
    throw svd_failure(i,1);
  }
}
