#include "arrayprocessing.h"
#include "mkl_complex_defined.h"
#include <iostream>

using namespace std;
namespace auxiliary{

void matrixprint(int n, int m, arcomplex<double> *array){
  //The second argument is always the contigous index
  //Prints col major
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      cout<<array[j+m*i]<<"\t";
    }
    cout<<"\n";
  }
 cout<<"END OF MATRIX\n";
}

//---------------------------------------------------------------------------------------------------//

//outdated and not well written, avoid usage (especially since there are zgesdd and zgesvd)
void lapackSVD(int MNumCols, int MNumRows, arcomplex<double> *Mnew, arcomplex<double> *Mnewcpy, double *diags){
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  char uplo=(MNumRows>=MNumCols)?'U':'L';
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  lapack_int info;
  double *offdiags=new double[containerDim-1];
  arcomplex<double> *QContainer=new arcomplex<double>[MNumRows*MNumRows];
  arcomplex<double> *PContainer=new arcomplex<double>[MNumCols*MNumCols];
  info=LAPACKE_zgebrd(LAPACK_COL_MAJOR,MNumRows,MNumCols,Mnew,MNumRows,diags,offdiags,QContainer,PContainer);
  //Okay, this is really nasty: Mnewcpy has to have a leading dimension of MNumCols, although it is (at input) of size MNumRows x MNumCols. zungbr does resize it, but it does not change the leading dimension -> manually take care of the correct leading dimension of the input matrix
  for(int mi=0;mi<MNumCols;++mi){
    for(int mip=0;mip<MNumRows;++mip){
      Mnewcpy[mip+mi*MNumCols]=Mnew[mip+mi*MNumRows];
    }
  }
  //I do not really know what ZUNGBR does with the uninitialized entries which have to be allocated since the matrices in the SVD can be larger than the original one - take CARE
  info=LAPACKE_zungbr(LAPACK_COL_MAJOR,'Q',MNumRows,MNumRows,MNumCols,Mnew,MNumRows,QContainer);
  info=LAPACKE_zungbr(LAPACK_COL_MAJOR,'P',MNumCols,MNumCols,MNumRows,Mnewcpy,MNumCols,PContainer);
  delete[] QContainer;
  delete[] PContainer;
  info=LAPACKE_zbdsqr(LAPACK_COL_MAJOR,uplo,containerDim,MNumCols,MNumRows,0,diags,offdiags,Mnewcpy,MNumCols,Mnew,MNumRows,0,1);
  delete[] offdiags;
}
}

