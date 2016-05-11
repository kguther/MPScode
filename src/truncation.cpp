//---------------------------------------------------------------------------------------------------//
// Order sortDatas with respect to lambdas, as it is used to truncate to the largest SVs
//---------------------------------------------------------------------------------------------------//

#include "truncation.h"

namespace auxiliary{

  bool compareSortData(sortData const &a, sortData const &b){
    return a.lambda>b.lambda;
  }

  bool compareSortDataFull(sortData const &a, sortData const &b){
    double const tol=1e-10;
    double buf=(a.lambda>b.lambda)?(a.lambda-b.lambda):(b.lambda-a.lambda);
    if(buf>tol)
      return a.lambda>b.lambda;
    if(a.QN.real()!=b.QN.real())
      return a.QN.real()>b.QN.real();
    return a.QN.imag()>b.QN.imag();
  }

  bool compareSortDataQNBased(sortData const &a, sortData const &b){
    if(a.QN.real()!=b.QN.real()){
      return a.QN.real()>b.QN.real();
    }
    return a.QN.imag()>b.QN.imag();  
    //The ordering has to be strictly deterministic, such that two arrays of sortData with the same entries of lambdas and QNs are ordered in the same way, even if some lambdas are degenerate
  }

}
