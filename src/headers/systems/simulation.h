#ifndef SIMULATION_CONTROL
#define SIMULATION_CONTROL

#include <complex>
#include <vector>
#include <string>
#include "parameters.h"
#include "problemOperators.h"
#include "templates/localMpo.h"
#include "Qsystem.h"

class simulation{
 public:
  simulation(problemParameters &pars, simulationParameters &simPars, double J, double g, double W, int pathPoints, double stepSize, double deltaPIn, std::string const &targetFile, int tSiteIn=-1, int jgScaleIn=1);
  void setMeasurement(mpo<std::complex<double> > &MPOperator, std::string &opName);
  void setLocalMeasurement(localMpo<std::complex<double> > &localMPOperator, std::string &opName);
  void setEntanglementMeasurement();
  void setEntanglementSpectrumMeasurement();
  std::vector<mpo<std::complex<double> > > measureTask;
  std::vector<localMpo<std::complex<double> > > localMeasureTask;
  std::vector<double> E0, dE;
  void run();
 private:
  //Order dependent, do not change
  problemParameters pars;
  simulationParameters simPars;
  Qsystem csystem;
  std::vector<double> convergedEigens;
  std::vector<std::string> operatorNames;
  std::vector<std::string> localOperatorNames;
  int pathLength;
  int tSite;
  int jgScale;
  int measureEE, measureES;
  double scaling, W, deltaP, targetDelta;
  std::complex<double> parDirection;
  std::string filePrefix;
  mpo<std::complex<double> > particleNumber;
  mpo<std::complex<double> > subChainParity;
  void singleRun();
};

#endif
