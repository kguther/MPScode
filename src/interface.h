#ifndef CONSOLE_UI
#define CONSOLE_UI

#include <string>

struct info{
  int L;
  int D;
  int par;
  int N;
  int nSweeps;
  int simType;
  int odd;
  int numPts;
  double rho;
  double alphaInit;
  double arpackTol;
  double arpackTolMin;
  double Jsc, gsc;
};

class interface{
 public: 
  interface();
  void provideInterface();
  void getScalingSerial(double J, double g);
  info parPack;
  std::string fileName;
 private:
  void readParFile(std::string const &fN);
};

#endif
