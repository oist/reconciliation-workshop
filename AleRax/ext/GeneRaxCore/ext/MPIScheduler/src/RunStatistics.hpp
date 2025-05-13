#pragma once

#include "Instance.hpp"
#include "Common.hpp"

#include <string>

namespace MPIScheduler {

class Logger;

class RunStatistics {
public:
  RunStatistics(const InstancesHistoric &historic, 
      Time begin, 
      Time end,
      unsigned int availableRanks,
      Logger &masterLogger);
  void printGeneralStatistics();
  void exportSVG(const std::string &svgfile);
private:
  const InstancesHistoric &_historic;
  Time _begin;
  Time _end;
  unsigned int _availableRanks;
  double _lbRatio;
  Logger &_masterLogger;
};

} // namespace MPIScheduler



