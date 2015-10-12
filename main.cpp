#include <iostream>
#include <complex>

using namespace std;

#define D 100
#define L 100

double solve(complex<double> ****cstate, complex<double> ****Hamiltonian, parameters pars);


int main(int *argc, char *argv[]){
  complex<double> ***test;
  create3D(2,2,2,&test);
  test[1][1][1]=1;
  cout<<test[1][1][1]<<endl;
  delete3D(&test);
  return 0;
}

double solve(complex<double> ****cstate, complex<double> ****Hamiltonian, parameters pars){return 0;}

