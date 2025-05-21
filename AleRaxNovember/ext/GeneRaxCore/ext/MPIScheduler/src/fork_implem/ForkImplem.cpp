#include "ForkImplem.hpp"
#include <cassert>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <iostream> 
#include "../SchedulerArgumentsParser.hpp"

using namespace std;

namespace MPIScheduler {

ForkRanksAllocator::ForkRanksAllocator(unsigned int availableRanks, 
    const string &execPath,
    const string &outputDir,
    const string &threadsArg):
  _coresInUse(0),
  _outputDir(outputDir),
  _execPath(execPath),
  _threadsArg(threadsArg)
{
  Common::makedir(Common::joinPaths(outputDir, "per_job_logs"));
  Common::makedir(Common::joinPaths(outputDir, "running_jobs"));
  _slots.push(Slot(0, availableRanks));

}

bool ForkRanksAllocator::ranksAvailable()
{
  return !_slots.empty();
}

bool ForkRanksAllocator::allRanksAvailable()
{
  return !_coresInUse;
}

static void split(const ForkRanksAllocator::Slot &parent,
    ForkRanksAllocator::Slot &son1,
    ForkRanksAllocator::Slot &son2,
    unsigned int son1size)
{
  son1 = ForkRanksAllocator::Slot(parent.startingRank, son1size);
  assert(parent.ranksNumber > son1size);
  son2 = ForkRanksAllocator::Slot(parent.startingRank + int(son1size), parent.ranksNumber - son1size);
}
  
InstancePtr ForkRanksAllocator::allocateRanks(unsigned int requestedRanks, 
      CommandPtr command)
{
  Slot slot = _slots.front();
  _slots.pop();
  while (slot.ranksNumber > requestedRanks) {
    Slot slot1, slot2;
    split(slot, slot1, slot2, slot.ranksNumber / 2); 
    slot = slot1;
    _slots.push(slot2);
  }
  
  
  shared_ptr<ForkInstance> instance(new ForkInstance(_outputDir,
    _execPath,
    slot.startingRank,
    int(slot.ranksNumber),
    command,
    _threadsArg));
  _runningInstances.insert(instance);
  _coresInUse += slot.ranksNumber;
  return  instance;
}

void ForkRanksAllocator::freeRanks(InstancePtr instance) 
{
  _coresInUse -= instance->getRanksNumber();
  _slots.push(Slot(instance->getStartingRank(), 
        (unsigned int)instance->getRanksNumber()));
}

vector<InstancePtr> ForkRanksAllocator::checkFinishedInstances()
{
  vector<InstancePtr> res;
  auto instanceIt = _runningInstances.begin();
  while (instanceIt != _runningInstances.end()) {
    auto instance = *instanceIt;
    if (instance->checkFinished()) {
      res.push_back(instance);
      if (instance->getReturnValue()) {
        instance->onFailure(instance->getReturnValue());
      }
      instanceIt = _runningInstances.erase(instanceIt);
    } else {
      instanceIt++;
    }
  }
  return res;
}

void ForkRanksAllocator::terminate()
{

}
  
ForkInstance::ForkInstance(const string &outputDir, 
      const string &execPath,
      int coresOffset, 
      int cores, 
      CommandPtr command,
      const string &threadsArg):
  Instance(command, coresOffset, cores, outputDir),
  _pid(0),
  _returnValue(0),
  _execPath(execPath),
  _threadsArg(threadsArg)
{

}

bool ForkInstance::execute(InstancePtr)
{
  _timer.reset();
  pid_t pid = fork();
  assert(pid >= 0);
  if (pid == 0) {
    int res = executeChild(_command, getOutputDir());
    exit(res);
  } else if (pid > 0) {
    _pid = pid; 
  }
  return true;
}

int ForkInstance::executeChild(const CommandPtr command, 
    const string &outputDir) 
{
  string logsFile = Common::joinPaths(outputDir, "per_job_logs", command->getId() + "_out.txt");
  string runningFile = Common::joinPaths(outputDir, "running_jobs", command->getId());
  ofstream os(runningFile);
  os << logsFile << endl;
  os.close();
  const vector<string> &args  = command->getArgs();
  string systemCommand = _execPath;
  for (auto &arg: args) {
    systemCommand = systemCommand + " " + arg;
  }
  if (_threadsArg.size()) {
    systemCommand += " " + _threadsArg + " " + to_string(_ranksNumber);
  }
  int result = systemCall(systemCommand, logsFile, true);
  remove(runningFile.c_str());
  return result;
}

bool ForkInstance::checkFinished()
{
  int status;
  pid_t result = waitpid(_pid, &status, WNOHANG);
  assert(result != -1);
  if (result) {
    _returnValue = WEXITSTATUS(status);
    setElapsedMs((int)_timer.getElapsedMs());
  }
  return result != 0;

}

}
