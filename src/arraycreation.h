#ifndef ARRAYCREATIONCUSTOM
#define ARRAYCREATIONCUSTOM
using namespace std; //BEWARE: FOR LAPACK ACCESS ALWAYS USE [0] AS SECOND INDEX FROM RIGHT, i.e. to access the third matrix in a 3D array in a lapack function, use array[2][0] - thats how lapack works

//------------------------------------------------------------------------------------------------------------------//
// This file contains template functions to dynamically allocate C arrays with 2 to 5 dimensions. The generated
// arrays are contigous, can be accessed like normal, i.e. static C arrays and have all other properties of
// generic C arrays. Use the corresponding delete function to free the memory when not needed anymore. 
// Save for the stateArray functions, none of them is used directly anymore.
//-----------------------------------------------------------------------------------------------------------------//

template<typename T> void create2D(const int dim1, const int dim2, T ***array);
template<typename T> void create3D(const int dim1, const int dim2, const int dim3, T ****array);
template<typename T> void create4D(const int dim1, const int dim2, const int dim3, const int dim4, T *****array);
template<typename T> void create5D(const int dim1, const int dim2, const int dim3, const int dim4, const int dim5,  T ******array);
template<typename T> void delete2D(T ***array);
template<typename T> void delete3D(T ****array);
template<typename T> void delete4D(T *****array);
template<typename T> void delete5D(T ******array);
template<typename T> void deleteStateArray(T *****array);

//---------------------------------------------------------------------------------------------------//

template<typename T> void create2D(const int dim1, const int dim2, T ***array){
  (*array)=new T*[dim1];
  (*array)[0]=new T[dim1*dim2];
  for(int i=1;i<dim1;i++){
    (*array)[i]=(*array)[i-1]+dim2;
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T> void create3D(const int dim1, const int dim2, const int dim3, T ****array){
  //iteratively builds up the 3D array from a 2D array of pointers - create2D is called by reference, thus, the adress of ***array is handed over
  create2D(dim1,dim2,array);
  (*array)[0][0]=new T[dim1*dim2*dim3];
  for(int i=0;i<dim1;i++){
    if(i>0){
      (*array)[i][0]=(*array)[i-1][0]+dim2*dim3;
    }
    for(int j=1;j<dim2;j++){
      (*array)[i][j]=(*array)[i][j-1]+dim3;
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T> void create4D(const int dim1, const int dim2, const int dim3, const int dim4, T *****array){
  create3D(dim1,dim2,dim3,array);
  (*array)[0][0][0]=new T[dim1*dim2*dim3*dim4];
  for(int i=0;i<dim1;i++){
    if(i>0){
      (*array)[i][0][0]=(*array)[i-1][0][0]+dim2*dim3*dim4;
    }
    for(int j=0;j<dim2;j++){
      if(j>0){
	(*array)[i][j][0]=(*array)[i][j-1][0]+dim3*dim4;
      }
      for(int k=1;k<dim3;k++){
	(*array)[i][j][k]=(*array)[i][j][k-1]+dim4;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T> void create5D(const int dim1, const int dim2, const int dim3, const int dim4, const int dim5, T ******array){
  create4D(dim1,dim2,dim3,dim4,array);
  (*array)[0][0][0][0]=new T[dim1*dim2*dim3*dim4*dim5];
  for(int i=0;i<dim1;i++){
    if(i>0){
      (*array)[i][0][0][0]=(*array)[i-1][0][0][0]+dim2*dim3*dim4*dim5;
    }
    for(int j=0;j<dim2;j++){
      if(j>0){
	(*array)[i][j][0][0]=(*array)[i][j-1][0][0]+dim3*dim4*dim5;
      }
      for(int k=0;k<dim3;k++){
	if(k>0){
	  (*array)[i][j][k][0]=(*array)[i][j][k-1][0]+dim4*dim5;
	}
	for(int l=1;l<dim4;l++){
	  (*array)[i][j][k][l]=(*array)[i][j][k][l-1]+dim5;
	}
      }
    }
  }
}

//---------------------------------------------------------------------------------------------------//

template<typename T> void delete2D(T ***array){
  delete[] (*array)[0];
  delete *array;
}

template<typename T> void delete3D(T ****array){
  delete[] (*array)[0][0];
  delete[] (*array)[0];
  delete *array;
}

template <typename T> void delete4D(T *****array){
  if(*array){
    if(**array){
      if(***array){
	if(****array){
	  delete[] (*array)[0][0][0];
	}
	delete[] (*array)[0][0];
      }
      delete[] (*array)[0];
    }
    delete *array;
  }
}

template <typename T> void delete5D(T ******array){
  delete[] (*array)[0][0][0][0];
  delete[] (*array)[0][0][0];
  delete[] (*array)[0][0];
  delete[] (*array)[0];
  delete *array;
}

template <typename T> void deleteStateArray(T *****array){
  delete4D(array);
}

#endif
