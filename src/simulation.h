#ifndef SIMULATION_CONTROL
#define SIMULATION_CONTROL

#include <complex>
#include <vector>
#include <string>
#include "parameters.h"
#include "problemOperators.h"
#include "localMpo.h"
#include "mps.h"
#include "network.h"

class simulation{
 public:
  simulation();
  simulation(problemParameters &pars, simulationParameters &simPars, double const J, double const g, int const pathPoints, std::string const &targetFile);
  void generate(problemParameters &pars, simulationParameters &simPars, double const J, double const g, int const pathPoints, std::string &targetFile);
  void setMeasurement(mpo<lapack_complex_double> &MPOperator, std::string &opName);
  void setLocalMeasurement(localMpo<lapack_complex_double> &localMPOperator, std::string &opName);
  std::vector<mpo<std::complex<double> > > measureTask;
  std::vector<localMpo<std::complex<double> > > localMeasureTask;
  std::vector<double> E0, dE;
  void run();
 private:
  network TensorNetwork;
  simulationParameters simPars;
  problemParameters pars;
  std::vector<double> convergedEigens;
  std::vector<std::string> operatorNames;
  std::vector<std::string> localOperatorNames;
  int pathLength;
  std::complex<double> parDirection;
  std::string filePrefix;
  mpo<lapack_complex_double> particleNumber;
  mpo<lapack_complex_double> subChainParity;
  void singleRun();
  int measure(mpo<std::complex<double> > *const MPOperator, double &expectationValue, mps *const MPState=0);
  int measureLocal(localMpo<std::complex<double> > *const localMPOperator, std::vector<std::complex<double> > &result, mps *const MPState=0);
  void initialize(problemParameters &pars, simulationParameters &simPars, double const J, double const g, int const pathPoints, std::string const &targetFile);
};

#endif
