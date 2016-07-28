#ifndef ARRAYPROCESSING
#define ARRAYPROCESSING

#include <arcomp.h>
#include <memory>

namespace auxiliary{

 void matrixprint(int n, int m, arcomplex<double> *array);
 void lapackSVD(int dim1, int dim2, arcomplex<double> *U, arcomplex<double> *V, double *diags);

//---------------------------------------------------------------------------------------------------//

template<typename T>
  void upperdiag(int dim1, int dim2, T *arrayin, T *arrayout, int ldim=0){
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

template<typename T>
  void lowerdiag(int dim1, int dim2, T *arrayin, T *arrayout, int ldim=0){
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

template<typename T>
 void arraycpy(int dim1, int dim2, T const *arraysource, T *arraytarget){
  for(int i=0;i<dim1;++i){
    for(int j=0;j<dim2;++j){
      arraytarget[j+dim2*i]=arraysource[j+dim2*i];
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
  void arraycpy(int dim, T const *arraysource, T *arraytarget){
  for(int i=0;i<dim;++i){
    arraytarget[i]=arraysource[i];
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
  void transp(int dim1, int dim2, T *array){
  //This switches from contigous column index to contigous row index and vice versa - quite effortive, try to avoid usage
  std::unique_ptr<T> tmpp(new T[dim1*dim2]);
  T* tmp=tmpp.get();
  arraycpy(dim1,dim2,array,tmp);
  for(int j=0;j<dim2;j++){
    for(int i=0;i<dim1;i++){
      array[i+j*dim1]=tmp[j+i*dim2];
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T>
  void matrixMult(int dimR, int dimL, int dimContract, T *a, T *b, T *result){
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
#endif
