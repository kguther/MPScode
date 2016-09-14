#ifndef SIMULATION_CONTROL
#define SIMULATION_CONTROL

#include <vector>
#include <string>
#include "parameters.h"
#include "problemOperators.h"
#include "templates/localMpo.h"
#include "Qsystem.h"
#include "mpstype.h"

class simulation{
 public:
  simulation(problemParameters &pars, simulationParameters &simPars, double J, double g, double W, int pathPoints, double stepSize, double deltaPIn, std::string const &targetFile, int tSiteIn=-1, int jgScaleIn=1);
  void setMeasurement(mpo<mpsEntryType > &MPOperator, std::string &opName);
  void setLocalMeasurement(localMpo<mpsEntryType > &localMPOperator, std::string &opName);
  void setEntanglementMeasurement();
  void setEntanglementSpectrumMeasurement();
  std::vector<mpo<mpsEntryType > > measureTask;
  std::vector<localMpo<mpsEntryType > > localMeasureTask;
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
  mpo<mpsEntryType > particleNumber;
  mpo<mpsEntryType > subChainParity;
  void singleRun();
};

#endif
