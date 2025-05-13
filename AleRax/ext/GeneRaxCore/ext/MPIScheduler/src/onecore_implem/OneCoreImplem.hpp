#pragma once

#include "../CommandsRunner.hpp"
#include <queue>


namespace MPIScheduler {

class OneCoreSlave {
public:
  OneCoreSlave():
    _masterRank(-1)
  {}
  ~OneCoreSlave();
  int main_core_slave(int argc, char **argv);

  int doWork(const CommandPtr command, 
    const std::string &outputDir) ;
private:
  void treatJobSlave();
  std::string _execPath;
  std::string _outputDir;
  int _masterRank;
  Timer _globalTimer;
  CommandsContainer _commands;
};

/*
 *  This allocator only returns single ranks
 */
class OneCoreRanksAllocator: public RanksAllocator {
public:
  // available threads must be a power of 2
  OneCoreRanksAllocator(unsigned int availableRanks, 
      const std::string &outoutDir);
  virtual ~OneCoreRanksAllocator() {}
  virtual bool ranksAvailable();
  virtual bool allRanksAvailable();
  virtual InstancePtr allocateRanks(unsigned int requestedRanks, 
      CommandPtr command);
  virtual void freeRanks(InstancePtr instance);
  virtual std::vector<InstancePtr> checkFinishedInstances();
  virtual void terminate();
private:
  unsigned int _cores;
  std::queue<int> _availableCores;
  std::string _outputDir;
  std::map<int, InstancePtr> _rankToInstances;

};

class OneCoreInstance: public Instance {
public:
  OneCoreInstance(const std::string &outputDir, 
      int rank, 
      CommandPtr command);

  virtual ~OneCoreInstance() {}
  virtual bool execute(InstancePtr self);
};

} // namespace MPIScheduler



