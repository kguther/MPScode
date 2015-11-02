#include "arrayprocessing.h"
#include <iostream>

void upperdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout){
  //for lapack postprocessing: extracts the upper triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(j<i){
	arrayout[j+dim2*i]=lapack_make_complex_double(0.0,0.0);
      }
      else{
	arrayout[j+dim2*i]=arrayin[j+dim2*i];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void lowerdiag(int const dim1, int const dim2, lapack_complex_double *arrayin, lapack_complex_double *arrayout){ 
  //for lapack postprocessing: extracts the lower triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(i<j){
	arrayout[j+dim2*i]=lapack_make_complex_double(0.0,0.0);
      }
      else{
	arrayout[j+dim2*i]=arrayin[j+dim2*i];
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void arraycpy(int const dim1, int const dim2, lapack_complex_double *arraysource, lapack_complex_double *arraytarget){
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      arraytarget[j+dim2*i]=arraysource[j+dim2*i];
    }
  }
}

//---------------------------------------------------------------------------------------------------//

void matrixprint(int const n, int const m, lapack_complex_double *array){
  //The second argument is always the contigous index
for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
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
