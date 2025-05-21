#pragma once

#include "../CommandsRunner.hpp"
#include <memory>
#include <unordered_set>
#include <queue>

namespace MPIScheduler {

class ForkInstance: public Instance {
public:
  ForkInstance(const std::string &outputDir, 
      const std::string &execPath,
      int coresOffset, 
      int cores, 
      CommandPtr command,
      const std::string &threadsArg);

  virtual ~ForkInstance() {}
  virtual bool execute(InstancePtr self);
  bool checkFinished();
  int getReturnValue() {return _returnValue;}
private:
  int executeChild(const CommandPtr command, 
    const std::string &outputDir); 
  int _pid;
  int _returnValue;
  std::string _execPath;
  Timer _timer;
  std::string _threadsArg;
};

class ForkRanksAllocator: public RanksAllocator {
public:
  // available threads must be a power of 2
  ForkRanksAllocator(unsigned int availableRanks, 
      const std::string &execPath,
      const std::string &outputDir,
      const std::string &threadsArg);
  virtual ~ForkRanksAllocator() {}
  virtual bool ranksAvailable();
  virtual bool allRanksAvailable();
  virtual InstancePtr allocateRanks(unsigned int requestedRanks, 
      CommandPtr command);
  virtual void freeRanks(InstancePtr instance);
  virtual std::vector<InstancePtr> checkFinishedInstances();
  virtual void terminate();
  struct Slot {
    Slot():
      startingRank(0),
      ranksNumber(0)
    {}
    Slot(int _startingRank, unsigned int _ranksNumber) : 
      startingRank(_startingRank),
      ranksNumber(_ranksNumber)
    {}
    int startingRank; 
    unsigned int ranksNumber;
  };
  std::queue<Slot> _slots;
private:
  std::unordered_set<std::shared_ptr<ForkInstance> > _runningInstances; 
  int _coresInUse;
  std::string _outputDir;
  std::string _execPath;
  std::string _threadsArg;
};


} // namespace MPIScheduler



