#ifndef CONSOLE_UI
#define CONSOLE_UI

#include <string>

//-------------------------------------------------------------------------------------------//
// The info struct contains all currently setable parameters.
//-------------------------------------------------------------------------------------------//

struct info{
  int L;
  int D;
  int par;
  int N;
  int nSweeps;
  int simType;
  int odd;
  int numPts;
  int nEigens;
  int nGs;
  double rho;
  double alphaInit;
  double arpackTol;
  double arpackTolMin;
  double Jsc, gsc, Wsc;
  double alphaMin, alphaMax;
  double scaling;
};

//-------------------------------------------------------------------------------------------//
// The interface class is used for input of parameters. This is done via parameter files
//-------------------------------------------------------------------------------------------//

class interface{
 public: 
  interface();
  void provideInterface(char *argv);
  void getScalingSerial(double J, double g);
  info parPack;
  std::string fileName;
 private:
  void readParFile(std::string const &fN);
};

#endif
