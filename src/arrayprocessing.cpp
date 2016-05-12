#include "arrayprocessing.h"
#include <iostream>

using namespace std;
namespace auxiliary{

void upperdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout, int ldim){
  //for lapack postprocessing: extracts the upper triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  if(ldim==0){
    ldim=dim2;
  }
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(j<i){
	arrayout[j+dim2*i]=0.0;
      }
      else{
	arrayout[j+dim2*i]=arrayin[j+ldim*i];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void lowerdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout, int ldim){ 
  //for lapack postprocessing: extracts the lower triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  if(ldim==0){
    ldim=dim2;
  }
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(i<j){
	arrayout[j+dim2*i]=0.0;
      }
      else{
	arrayout[j+dim2*i]=arrayin[j+ldim*i];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void arraycpy(int const dim1, int const dim2, lapack_complex_double *arraysource, lapack_complex_double *arraytarget){
  for(int i=0;i<dim1;++i){
    for(int j=0;j<dim2;++j){
      arraytarget[j+dim2*i]=arraysource[j+dim2*i];
    }
  }
}

//---------------------------------------------------------------------------------------------------//
void arraycpy(int const dim, lapack_complex_double *arraysource, lapack_complex_double *arraytarget){
  for(int i=0;i<dim;++i){
    arraytarget[i]=arraysource[i];
  }
}

//---------------------------------------------------------------------------------------------------//

void matrixprint(int const n, int const m, lapack_complex_double *array){
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

void transp(int const dim1, int const dim2, lapack_complex_double *array){                                                
  //This switches from contigous column index to contigous row index and vice versa - quite effortive, try to avoid usage
  lapack_complex_double *tmp=new lapack_complex_double [dim1*dim2];
  arraycpy(dim1,dim2,array,tmp);
  for(int j=0;j<dim2;j++){
    for(int i=0;i<dim1;i++){
      array[i+j*dim1]=tmp[j+i*dim2];
    }
  }
  delete[] tmp;
}

//---------------------------------------------------------------------------------------------------//

void lapackSVD(int MNumCols, int MNumRows, lapack_complex_double *Mnew, lapack_complex_double *Mnewcpy, double *diags){
  int containerDim=(MNumRows>MNumCols)?MNumCols:MNumRows;
  char uplo=(MNumRows>=MNumCols)?'U':'L';
  int maxDim=(MNumRows>MNumCols)?MNumRows:MNumCols;
  lapack_int info;
  double *offdiags=new double[containerDim-1];
  lapack_complex_double *QContainer=new lapack_complex_double[MNumRows*MNumRows];
  lapack_complex_double *PContainer=new lapack_complex_double[MNumCols*MNumCols];
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

//---------------------------------------------------------------------------------------------------//
  void matrixMult(int dimL, int dimR, int dimContract, lapack_complex_double *a, lapack_complex_double *b, lapack_complex_double *result){
    for(int m=0;m<dimR*dimL;++m){
      result[m]=0;
    }
    for(int mL=0;mL<dimL;++mL){
      for(int mR=0;mR<dimR;++mR){
	for(int mContract=0;mContract<dimContract;++mContract){
	  result[mL+mR*dimL]+=a[mL+mContract*dimL]*b[mContract+mR*dimContract];
	}
      }
    }
  }

}
