

#include "Command.hpp"
#include <cassert>

namespace MPIScheduler {
  

Command::Command(const std::string &id, 
    unsigned int ranks,
    long estimatedCost,
    const std::vector<std::string> &arguments):
  _id(id),
  _args(arguments),
  _ranksNumber(ranks),
  _estimatedCost(estimatedCost)
{
  assert(_ranksNumber > 0);
}

Command::~Command()
{

}

}
