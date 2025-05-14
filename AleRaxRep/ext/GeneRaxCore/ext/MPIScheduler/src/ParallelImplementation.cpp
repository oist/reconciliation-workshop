#include "ParallelImplementation.hpp"

#include <assert.h>

#ifdef WITH_MPI
#include "split_implem/SplitImplem.hpp"
#include "onecore_implem/OneCoreImplem.hpp"
#include <mpi.h>
#endif
#include "fork_implem/ForkImplem.hpp"

#include "Command.hpp"
#include "Common.hpp"
#include "SchedulerArgumentsParser.hpp"

using namespace std;

namespace MPIScheduler {

ParallelImplementation::ParallelImplementation(const string &implem): _rank(0), _ranksNumber(0), _ownMPIContext(false) {
  if (implem == "--split-scheduler")
    _impl = split;
  else if (implem == "--onecore-scheduler") 
    _impl = onecore;
  else if (implem == "--fork-scheduler")
    _impl = fork;
  else
    _impl = invalid;
}

bool ParallelImplementation::isValid() const {
#ifndef WITH_MPI
  if (isMPI()) { 
    cerr << "Error: trying to use an MPI implementation that was not compiled" << endl;
    return false;
  }
#endif
  return _impl != invalid;
}

int ParallelImplementation::getRank() const {return _rank;}
unsigned int ParallelImplementation::getRanksNumber() const {return _ranksNumber;}

bool ParallelImplementation::isMPI() const {
  return (_impl == split) || (_impl == onecore); 
}

void ParallelImplementation::initParallelContext(int argc, char **argv, void *comm, int ranks) {
  if (isMPI()) {
#ifdef WITH_MPI
    _ownMPIContext = (comm == 0);
    if (_ownMPIContext) {
      Common::check(MPI_Init(&argc, &argv));
      _comm = MPI_COMM_WORLD;
    } else {
      _comm = *((MPI_Comm*)comm);
    }
    
    MPI_Comm_rank(_comm, &_rank);
    int temp;
    MPI_Comm_size(_comm, &temp);
    _ranksNumber = static_cast<unsigned int>(temp);
    assert(static_cast<unsigned int>(ranks) == _ranksNumber);
#else
    assert(0);
#endif
  } else {
    _rank = 0;
    _ranksNumber = ranks;
  }
}

void ParallelImplementation::closeParallelContext() {
  if (isMPI()) {
#ifdef WITH_MPI
    if (_ownMPIContext) {
      MPI_Finalize();
    }
#else
    assert(0);
#endif
  }
}

bool ParallelImplementation::slavesToStart() const {
  return isMPI();
}

void ParallelImplementation::startSlaves(int argc, char **argv) {
  if (getRank() != int(getRanksNumber()) - 1) {
#ifdef WITH_MPI
    if (_impl == split) {
      SplitSlave slave;
      slave.main_split_slave(argc, argv);
    } else if (_impl == onecore) {
      OneCoreSlave slave;
      slave.main_core_slave(argc, argv);
    }
#else
    (void*)argc;
    (void*)argv;
#endif
  }
}

shared_ptr<RanksAllocator> ParallelImplementation::getRanksAllocator(SchedulerArgumentsParser &arg,
                                  unsigned int ranksNumber) {
#ifdef WITH_MPI
  if (_impl == split) {
    return shared_ptr<RanksAllocator>(new SplitRanksAllocator(ranksNumber,
        arg.outputDir));
  } else if (_impl == onecore) {
    return shared_ptr<RanksAllocator>(new OneCoreRanksAllocator(ranksNumber,
        arg.outputDir));
  }
#endif
  if (_impl == fork) {
    return shared_ptr<RanksAllocator>(new ForkRanksAllocator(ranksNumber,  arg.library, arg.outputDir, arg.threadsArg));
  }
  assert(0);
  return shared_ptr<RanksAllocator>(0);
}

}


