#pragma once

#include "Command.hpp"
#include "Instance.hpp"

namespace MPIScheduler {

class RanksAllocator {
public:
  RanksAllocator() = default;
  virtual ~RanksAllocator() = default;
  virtual bool ranksAvailable() = 0;
  virtual bool allRanksAvailable() = 0;
  virtual InstancePtr allocateRanks(unsigned int requestedRanks, 
      CommandPtr command) = 0;
  virtual void freeRanks(InstancePtr instance) = 0;
  virtual std::vector<InstancePtr> checkFinishedInstances() = 0;
  virtual void terminate() {}
  virtual void preprocessCommand(CommandPtr /*command*/) {}
};

}
