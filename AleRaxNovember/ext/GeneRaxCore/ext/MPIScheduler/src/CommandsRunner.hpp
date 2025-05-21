#pragma once

#include "Command.hpp"
#include "Instance.hpp"
#include "RanksAllocator.hpp"
#include "Logger.hpp"
#include <string>
#include <vector>

namespace MPIScheduler {

class SVGDrawer;
class Logger;

class CommandsContainer {
public:
  CommandsContainer() {}
  explicit CommandsContainer(const std::string &commandsFilename,
      bool addFakeExecutableName);

  CommandPtr getCommand(std::string id) const;
  std::vector<CommandPtr> &getCommands() {return _commands;}
  const std::vector<CommandPtr> &getCommands() const {return _commands;}
private:
  void addCommand(CommandPtr command);

  std::vector<CommandPtr> _commands;
  std::map<std::string, CommandPtr> _dicoCommands;
};

class CommandsRunner {
public:
  CommandsRunner(const CommandsContainer &commandsContainer,
      std::shared_ptr<RanksAllocator> allocator,
      const std::string &outputDir,
      bool jobFailureFatal,
      Logger &masterLogger);
  void run(bool isMPI);
  const InstancesHistoric &getHistoric() const {return _historic;} 
private:
  
  static bool compareCommands(CommandPtr c1, CommandPtr c2);
  CommandPtr getPendingCommand() {return *_commandIterator;}
  bool isCommandsEmpty() {return _commandIterator == _commandsVector.end();}
  
  bool executePendingCommand();
  void onFinishedInstance(InstancePtr instance);
  
  const std::string _outputDir;

  std::shared_ptr<RanksAllocator> _allocator;
  std::vector<CommandPtr> _commandsVector;
  std::vector<CommandPtr>::iterator _commandIterator;
  Checkpoint _checkpoint;
  InstancesHistoric _historic;
  int _finishedInstancesNumber;
  bool _verbose;
  bool _jobFailureFatal;
  Logger &_masterLogger;
};

}
