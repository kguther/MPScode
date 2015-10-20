#ifndef PARAMETERS
#define PARAMETERS

class parameters;

class parameters{
 public:
  parameters();
  parameters(int din, int Din, int Lin, int Dwin, int Nin);
  parameters(int din, int Din, int Lin, int Dwin, int Nin, int neigsin);
  int D, L, N, d, Dw, neigs;
};

#endif
