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
  int nStages;
  int nGs;
  int Dw;
  int tPos;
  double rho;
  double alphaInit;
  double arpackTol;
  double arpackTolMin;
  double Jsc, gsc, Wsc;
  double alphaMin, alphaMax;
  double scaling;
  double delta;
  double acc;
  double tReal, tImag;
};

int symmetryBroken(info parPack);

//-------------------------------------------------------------------------------------------//
// The interface class is used for input of parameters. This is done via parameter files
//-------------------------------------------------------------------------------------------//

class interface{
 public: 
  interface();
  void provideInterface(char *argv);
  info parPack;
  std::string fileName;
 private:
  void readParFile(std::string const &fN);
};

#endif
