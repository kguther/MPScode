#ifndef ARRAYPROCESSING
#define ARRAYPROCESSING 

#ifndef complex_type_plh
#define complex_type_plh lapack_complex_double
#endif

using namespace std;

void upperdiag(int dim1, int dim2, complex_type_plh *arrayin, complex_type_plh *arrayout); //for lapack postprocessing
void lowerdiag(int dim1, int dim2, complex_type_plh *arrayin, complex_type_plh *arrayout);
void arraycpy(int dim1, int dim2, complex_type_plh *arraysource, complex_type_plh *arraytarget);
void transp(int dim1, int dim2, complex_type_plh  **array);
void matrixprint(int m, complex_type_plh *array);

void upperdiag(int dim1, int dim2, complex_type_plh *arrayin, complex_type_plh *arrayout){ //for lapack postprocessing: extracts the upper triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(j<i){
	arrayout[j+dim2*i]=0;
      }
      else{
	arrayout[j+dim2*i]=arrayin[j+dim2*i];
      }
    }
  }
}

void lowerdiag(int dim1, int dim2, complex_type_plh *arrayin, complex_type_plh *arrayout){ //for lapack postprocessing: extracts the lower triangular part of some matrix arrayin, pastes it into arrayout and fills up arrayout with zeros
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      if(i<j){
	(arrayout)[j+dim2*i]=0;
      }
      else{
	(arrayout)[j+dim2*i]=arrayin[j+dim2*i];
      }
    }
  }
}

void arraycpy(int dim1, int dim2, complex_type_plh *arraysource, complex_type_plh *arraytarget){
  for(int i=0;i<dim1;i++){
    for(int j=0;j<dim2;j++){
      arraytarget[j+dim2*i]=arraysource[j+dim2*i];
    }
  }
}

void conjugate(int dim1, int dim2, complex_type_plh **array){
  complex_type_plh *tmp=new complex_type_plh[dim1*dim2];
  for(int j=0;j<dim2;j++){}
}


void matrixprint(int m, complex_type_plh *array){
for(int i=0;i<m;i++){
    for(int j=0;j<m;j++){
      cout<<i<<","<<j<<"\t"<<array[j+m*i]<<endl;
    }
  }
}

void transp(int dim1, int dim2, complex_type_plh *array){                                                //This switches from contigous column index to contigous row index and vice versa - quite effortive, try to avoid usage
  complex_type_plh *tmp;
  arraycpy(dim1,dim2,array,tmp);
  for(int j=0;j<dim2;j++){
    for(int i=0;i<dim1;i++){
      array[i+j*dim1]=tmp[j+i*dim2];
    }
  }
}


#endif
