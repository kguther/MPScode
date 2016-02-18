#ifndef DIMENSION_TABLE
#define DIMENSION_TABLE

#include <vector>
#include "localHSpaces.h"

//---------------------------------------------------------------------------------------------------//
// This really simple, really useful class contains the information on the local bond dimensions of an MPS.
//---------------------------------------------------------------------------------------------------//

class dimensionTable{
 public:
  dimensionTable();
  dimensionTable(int const Din, int const Lin, localHSpaces din);
  int d() const {return dpar;}
  int D() const {return Dpar;}
  int L() const {return Lpar;}
  int iCrit() const {return icrit;}
  localHSpaces const & locdTable() const {return dpars;}
  int locDimR(int const i) const;
  int locDimL(int const i) const;
  int locDMax(int const i) const;
  int locd(int const i) const;
  void setParameterD(int Dnew);
 private:
  int dpar, Dpar, Lpar;
  localHSpaces dpars;
  int icrit;
  void getIcrit();
  void getDMaxTable();
  std::vector<int> DMaxTable;
};

#endif
