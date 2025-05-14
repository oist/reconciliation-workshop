#include "CommandsRunner.hpp"

#include "Common.hpp"
#include "Logger.hpp"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <unistd.h>

using namespace std;

namespace MPIScheduler {

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
        }).base(), s.end());
}

// Read a line in a commands file and skip comments
// discards empty lines
static bool readNextLine(ifstream &is, string &os)
{
  while (getline(is, os)) {
    auto end = os.find("#");
    if (string::npos != end)
      os = os.substr(0, end);
    rtrim(os);
    if (os.size()) 
      return true;
  }
  return false;
}

CommandsContainer::CommandsContainer(const string &commandsFilename,
    bool addFakeExecutableName)
{
  ifstream reader(commandsFilename);
  if (!reader)
    throw MPISchedulerException("Cannot open commands file ", commandsFilename);
  
  string line;
  while (readNextLine(reader, line)) {
    
    string id;
    unsigned int ranks;
    long estimatedCost;
    
    istringstream iss(line);
    iss >> id;
    iss >> ranks;
    iss >> estimatedCost;
    
    if (ranks == 0) {
      cerr << "[mpi-scheduler warning] Found a job with 0 ranks... will assign 1 rank instead" << endl;
      ranks = 1;
    }

    vector<string> commandVector;
    if (addFakeExecutableName) {
      commandVector.push_back("mpi-scheduler");
    }
    while (!iss.eof()) {
      string argument;
      iss >> argument;
      commandVector.push_back(argument);
    }
    CommandPtr command(new Command(id, ranks, estimatedCost, commandVector));
    addCommand(command);
  }
}

void CommandsContainer::addCommand(CommandPtr command)
{
  _commands.push_back(command); 
  _dicoCommands[command->getId()] = command;
}

CommandPtr CommandsContainer::getCommand(string id) const
{
  auto res =  _dicoCommands.find(id);
  if (res == _dicoCommands.end())
    return CommandPtr(0);
  else 
    return res->second;
}

CommandsRunner::CommandsRunner(const CommandsContainer &commandsContainer,
      shared_ptr<RanksAllocator> allocator,
      const string &outputDir,
      bool jobFailureFatal,
      Logger &masterLogger):
  _outputDir(outputDir),
  _allocator(allocator),
  _checkpoint(outputDir),
  _finishedInstancesNumber(0),
  _verbose(true),
  _jobFailureFatal(jobFailureFatal),
  _masterLogger(masterLogger)
{
  Instance::_jobFailureFatal = _jobFailureFatal;
  _masterLogger.getCout() << "The master process runs on node " << Common::getHost() 
       << " and on pid " << Common::getPid() << endl;
  for (auto command: commandsContainer.getCommands()) {
    if (!_checkpoint.isDone(command->getId())) {
      _commandsVector.push_back(command);
    }
  }
  sort(_commandsVector.begin(), _commandsVector.end(), compareCommands);
  _commandIterator = _commandsVector.begin();
  _masterLogger.getCout() << "Remaining commands: " << _commandsVector.size() << endl;

}

void CommandsRunner::run(bool isMPI) 
{
  Timer globalTimer;
  Timer minuteTimer;
  while (!_allocator->allRanksAvailable() || !isCommandsEmpty()) {
    if (minuteTimer.getElapsedMs() > 1000 * 60) {
      _masterLogger.getCout() << "Runner is still alive after " << globalTimer.getElapsedMs() / 1000 << "s" << endl;
      minuteTimer.reset();
    } 
    if (!isCommandsEmpty()) {
      if (_allocator->ranksAvailable()) {
        executePendingCommand();
      }
    }
    vector<InstancePtr> finishedInstances = _allocator->checkFinishedInstances();
    for (auto instance: finishedInstances) {
      onFinishedInstance(instance);
    }
    if (!isMPI) {
      usleep(50); // ms
    }
  }
}
  


bool CommandsRunner::compareCommands(CommandPtr c1, CommandPtr c2)
{
  if (c2->getRanksNumber() == c1->getRanksNumber()) {
    return c2->getEstimatedCost() < c1->getEstimatedCost();
  }
  return c2->getRanksNumber() < c1->getRanksNumber();
}

  
bool CommandsRunner::executePendingCommand()
{
  auto command = getPendingCommand();
  InstancePtr instance = _allocator->allocateRanks(command->getRanksNumber(), command);
  if (!instance->execute(instance)) {
    _masterLogger.getCout() << "Failed to start " << command->getId() << ". Will retry later " << endl;
    _allocator->freeRanks(instance);
    return false;
  }
    
  if (_verbose) {
    _masterLogger.getCout() << "## Started " << command->getId() << " on [" 
      << instance->getStartingRank()  << ":"
      << instance->getStartingRank() + instance->getRanksNumber() - 1 
      << "] "
      << endl;
  }
  _historic.push_back(instance);
  ++_commandIterator;
  if (isCommandsEmpty()) {
    _masterLogger.getCout() << "All commands were launched" << endl;
  }
  return true;
}

void CommandsRunner::onFinishedInstance(InstancePtr instance)
{
  instance->onFinished();
  _checkpoint.markDone(instance->getId());
  if (_verbose) {
    _masterLogger.getCout() << "End of " << instance->getId() << " after " <<  instance->getElapsedMs() << "ms ";
    _masterLogger.getCout() << " (" << ++_finishedInstancesNumber << "/" << _commandsVector.size() << ")" << endl;
  }
 _allocator->freeRanks(instance);
}

}
