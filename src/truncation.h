#ifndef AUXILIARY_FUNCTIONS_FOR_SV_BASED_TRUNCATION
#define AUXILIARY_FUNCTIONS_FOR_SV_BASED_TRUNCATION

#include <complex>
#include <vector>

namespace auxiliary{

struct sortData{
  double lambda;
  std::vector<std::complex<int> > QN;
  int index;
  int indexExp;
};

 bool compareSortData(sortData const &a, sortData const &b);
 bool compareSortDataFull(sortData const &a, sortData const &b);
 bool compareSortDataQNBased(sortData const &a, sortData const &b);

}

#endif
