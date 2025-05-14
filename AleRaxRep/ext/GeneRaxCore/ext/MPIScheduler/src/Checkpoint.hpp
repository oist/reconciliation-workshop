#pragma once

#include <string>
#include <fstream>
#include <set>

namespace MPIScheduler {

class Checkpoint {
public:
  explicit Checkpoint(const std::string &outputDir); 
 
  bool isDone(const std::string &id);
  void markDone(const std::string &id);

private:
  static std::string getCheckpointCommandsFile(const std::string &outputDir);

private:
  std::set<std::string> _ids; 
  std::ofstream _os;
};

} // namespace MPIScheduler


