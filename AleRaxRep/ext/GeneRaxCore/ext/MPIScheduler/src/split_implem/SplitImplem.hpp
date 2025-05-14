#pragma once

#include "../CommandsRunner.hpp"
#include <queue>

namespace MPIScheduler {

int main_split_master(int argc, char **argv);
int main_split_slave(int argc, char **argv);
void terminate_slaves();




class SplitSlave {
public:
  SplitSlave(): 
    _localRank(-1),
    _globalRank(-1),
    _localMasterRank(-1),
    _globalMasterRank(-1)
  {}
  ~SplitSlave();
  int main_split_slave(int argc, char **argv);
private:
  int doWork(const CommandPtr command, MPI_Comm workersComm, const std::string &outputDir); 
  void splitSlave();
  void treatJobSlave();
  void terminateSlave();
  MPI_Comm _localComm; 
  int _localRank;
  int _globalRank;
  int _localMasterRank;
  int _globalMasterRank;
  CommandsContainer _commands;
  std::string _outputDir;
  Timer _globalTimer;
  std::string _libraryPath;
};

/*
 *  This allocator assumes that each request do
 *  not ask more ranks than the previous one and
 *  that the requests are powers of 2.
 */
class SplitRanksAllocator: public RanksAllocator {
public:
  // available threads must be a power of 2
  SplitRanksAllocator(unsigned int availableRanks, 
      const std::string &outoutDir);
  virtual ~SplitRanksAllocator() {}
  virtual bool ranksAvailable();
  virtual bool allRanksAvailable();
  virtual InstancePtr allocateRanks(unsigned int requestedRanks, 
      CommandPtr command);
  virtual void freeRanks(InstancePtr instance);
  virtual std::vector<InstancePtr> checkFinishedInstances();
  virtual void terminate();
  virtual void preprocessCommand(CommandPtr cmd);
  
  struct Slot {
    Slot():
      startingRank(0),
      ranksNumber(0)
    {}
    Slot(int _startingRank, unsigned int _ranksNumber) : 
      startingRank(_startingRank),
      ranksNumber(_ranksNumber)
    {}

    int startingRank; // relative to MPI_COMM_WORLD
    unsigned int ranksNumber;
  };

private:
  unsigned int _totalRanks;
  std::queue<Slot> _slots;
  int _ranksInUse;
  std::string _outputDir;
  std::map<int, InstancePtr> _rankToInstances;
};

class SplitInstance: public Instance {
public:
  SplitInstance(const std::string &outputDir, 
      int startingRank, 
      int ranksNumber,
      CommandPtr command);

  virtual ~SplitInstance() {}
  virtual bool execute(InstancePtr self);
};

} // namespace MPIScheduler


