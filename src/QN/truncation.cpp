//---------------------------------------------------------------------------------------------------//
// Order sortDatas with respect to lambdas, as it is used to truncate to the largest SVs
//---------------------------------------------------------------------------------------------------//

#include "truncation.h"
#include <cmath>

namespace auxiliary{

  bool compareSortData(sortData const &a, sortData const &b){
    double const tol=1e-12;
    if(std::abs(a.lambda-b.lambda)>tol)
      return a.lambda>b.lambda;
    //if exactly one of a and b is unassigned, the assigned one shall have precedence -> the first lDR entries in left/rightEnrichmentBlockwise are always assigned
    return a.QN.size()>b.QN.size();
  }

  bool compareSortDataFull(sortData const &a, sortData const &b){
    double const tol=1e-10;
    double buf=(a.lambda>b.lambda)?(a.lambda-b.lambda):(b.lambda-a.lambda);
    if(buf>tol)
      return a.lambda>b.lambda;
    if(a.QN[0].real()!=b.QN[0].real())
      return a.QN[0].real()>b.QN[0].real();
    return a.QN[0].imag()>b.QN[0].imag();
  }

  bool compareSortDataQNBased(sortData const &a, sortData const &b){
    if(a.QN[0].real()!=b.QN[0].real()){
      return a.QN[0].real()>b.QN[0].real();
    }
    return a.QN[0].imag()>b.QN[0].imag();  
    //The ordering has to be strictly deterministic, such that two arrays of sortData with the same entries of lambdas and QNs are ordered in the same way, even if some lambdas are degenerate
  }

}
