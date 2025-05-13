#pragma once

#include <memory>
#include <string>

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace MPIScheduler {

class RanksAllocator;
class SchedulerArgumentsParser;

class ParallelImplementation {
  enum Impl {
    split,
    onecore,
    fork,
    invalid
  };
public:
  explicit ParallelImplementation(const std::string &implem);
  bool isValid() const;
  int getRank() const;
  unsigned int getRanksNumber() const;
  bool isMPI() const;
  void initParallelContext(int argc, char **argv, void *comm, int ranks);
  void closeParallelContext();
  bool slavesToStart() const;
  void startSlaves(int argc, char **argv);
  std::shared_ptr<RanksAllocator> getRanksAllocator(SchedulerArgumentsParser &arg,
                                    unsigned int ranksNumber);
  bool addFakeExecutableName() {return _impl == split;}
private:
  Impl _impl;
  int _rank;
  unsigned int _ranksNumber;
  bool _ownMPIContext;
#ifdef WITH_MPI
  MPI_Comm _comm;
#endif
};

}

