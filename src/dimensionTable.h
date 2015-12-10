#ifndef DIMENSION_TABLE
#define DIMENSION_TABLE

class dimensionTable{
 public:
  dimensionTable();
  dimensionTable(int const din, int const Din, int const Lin);
  int d() const {return dpar;}
  int D() const {return Dpar;}
  int L() const {return Lpar;}
  int locDimR(int const i);
  int locDimL(int const i);
  int locDMax(int const i);
  int locd(int const i);
  void initialize(int const din, int const Din, int const Lin);
  void setParameterD(int Dnew);
 private:
  int dpar, Dpar, Lpar;
  int icrit;
  void getIcrit();
};

#endif
