#pragma once

#include "Common.hpp"

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <stack>

#include "Checkpoint.hpp"

namespace MPIScheduler {


class Command {
public:
  Command(const std::string &id, 
      unsigned int ranks,
      long estimatedCost,
      const std::vector<std::string> &arguments);
  virtual ~Command();
 
  const std::string &getId() const {return _id;}
  long getEstimatedCost() const {return _estimatedCost;}
  void setRanksNumber(unsigned int ranks) {_ranksNumber = ranks;}
  unsigned int getRanksNumber() const {return _ranksNumber;}
  std::string toString() const;
  const std::vector<std::string> &getArgs() const {return _args;}
public:
  // initial information
  const std::string _id;
  const std::vector<std::string> _args;
  unsigned int _ranksNumber;
  long _estimatedCost;
};
using CommandPtr = std::shared_ptr<Command>;

}
