#include "SplitImplem.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "../SchedulerArgumentsParser.hpp"
#include "../DynamicLibrary.hpp"
#include <cassert>

using namespace std;

namespace MPIScheduler {

const int SIGNAL_SPLIT = 1;
const int SIGNAL_JOB = 2;
const int SIGNAL_TERMINATE = 3;

const int TAG_END_JOB = 1;
const int TAG_START_JOB = 2;
const int TAG_SPLIT = 3;
const int TAG_MASTER_SIGNAL = 4;

const int MSG_SIZE_END_JOB = 3;


static int getRank(MPI_Comm comm) {
  int rank = 0;
  Common::check(MPI_Comm_rank(comm, &rank));
  return rank;
}

static int getSize(MPI_Comm comm) {
  int size = 0;
  Common::check(MPI_Comm_size(comm, &size));
  return size;
}


static int mpiSend(void* data,
    int count,
    MPI_Datatype datatype,
    int destination,
    int tag,
    MPI_Comm communicator)
{
  return MPI_Send(data, count, datatype, destination, tag, communicator);
}

static int mpiRecv(        void* data,
    int count,
    MPI_Datatype datatype,
    int source,
    int tag,
    MPI_Comm communicator,
    MPI_Status* status)
{
  return MPI_Recv(data, count, datatype, source, tag, communicator, status);
}

static int mpiBcast( void *buffer, int count, MPI_Datatype datatype, int root, 
                   MPI_Comm comm)
{
  return MPI_Bcast(buffer, count, datatype, root, comm);
}

static int mpiIprobe(int source, int tag, MPI_Comm comm, int *flag,
        MPI_Status *status)
{
  return MPI_Iprobe(source, tag, comm, flag, status);
}

static int mpiSplit(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
  return MPI_Comm_split(comm, color, key, newcomm);
}

SplitSlave::~SplitSlave()
{
}

static void move_file(const string &from, const string &to)
{
  if (std::rename(from.c_str(), to.c_str()) != 0) {
    cout << "failed to move " << from << " to " << to << endl;
    perror( "perror:" );
  }
}

int SplitSlave::doWork(const CommandPtr command, 
    MPI_Comm workersComm,
    const string &outputDir) 
{
  bool isMaster = !getRank(workersComm);
  string logsFile = Common::joinPaths(outputDir, 
      "running_jobs", command->getId() + "_out.txt");
  string errFile = Common::joinPaths(outputDir, 
      "running_jobs", command->getId() + "_err.txt");
  auto library = DynamicLibrary::getLibrary(_libraryPath);
  MPI_Comm raxmlComm;
  Common::check(MPI_Comm_dup(workersComm, &raxmlComm));
  int res = library->run(logsFile, errFile, command->getArgs(), (void*)raxmlComm);
  Common::check(MPI_Barrier(raxmlComm));
  MPI_Comm_free(&raxmlComm);
  if (isMaster) {
    move_file(logsFile,  Common::joinPaths(outputDir, 
          "per_job_logs", command->getId() + "_out.txt"));
    move_file(errFile, Common::joinPaths(outputDir, 
          "per_job_logs", command->getId() + "_err.txt"));
  }
  return res;
}


void SplitSlave::splitSlave() 
{
  MPI_Status status;
  int splitSize;
  if (_localMasterRank == _localRank) {
    Common::check(mpiRecv(&splitSize, 1, MPI_INT, _globalMasterRank, TAG_SPLIT, MPI_COMM_WORLD, &status));
    int signal = SIGNAL_SPLIT;
    Common::check(mpiBcast(&signal, 1, MPI_INT, _localMasterRank, _localComm));
  }
  Common::check(mpiBcast(&splitSize, 1, MPI_INT, _localMasterRank, _localComm));
  bool inFirstSplit = (_localRank < splitSize);
  int newRank = inFirstSplit ? _localRank : (_localRank - splitSize);
  MPI_Comm newComm;
  Common::check(mpiSplit(_localComm, !inFirstSplit , newRank, &newComm));
  _localComm = newComm;
  _localRank = newRank;
}

void SplitSlave::treatJobSlave()
{
  bool alone = getSize(_localComm) == 1;
  MPI_Status status;
  const int maxCommandSize = 200;
  char command[maxCommandSize];
  if (_localMasterRank == _localRank) {
    Common::check(mpiRecv(command, maxCommandSize, MPI_CHAR, _globalMasterRank, TAG_START_JOB, MPI_COMM_WORLD, &status));
    if (!alone) {
      int signal = SIGNAL_JOB;
      Common::check(mpiBcast(&signal, 1, MPI_INT, _localMasterRank, _localComm));
    }
  }

  if (!alone) {
    Common::check(mpiBcast(command, maxCommandSize, MPI_CHAR, _localMasterRank, _localComm));
  }
  Timer timer;
  auto startingTime = _globalTimer.getElapsedMs();
  auto jobResult = doWork(_commands.getCommand(string(command)), _localComm, _outputDir);
  auto elapsedMS = timer.getElapsedMs();
  if (_localMasterRank == _localRank) {
    int endJobMsg[MSG_SIZE_END_JOB];
    endJobMsg[0] = jobResult;
    endJobMsg[1] = int(startingTime);
    endJobMsg[2] = int(elapsedMS);
    Common::check(mpiSend(endJobMsg, MSG_SIZE_END_JOB, MPI_INT, _globalMasterRank, TAG_END_JOB, MPI_COMM_WORLD));
  }
}

void SplitSlave::terminateSlave()
{
  if (_localMasterRank == _localRank) {
    int signal = SIGNAL_TERMINATE;
    Common::check(mpiBcast(&signal, 1, MPI_INT, _localMasterRank, _localComm));
  }
}

int SplitSlave::main_split_slave(int argc, char **argv)
{
  SchedulerArgumentsParser arg(argc, argv);
  _libraryPath = arg.library;
  _outputDir = arg.outputDir;
  _commands = CommandsContainer(arg.commandsFilename, true);
  _globalRank = getRank(MPI_COMM_WORLD);
  _localMasterRank = 0;
  _globalMasterRank = getSize(MPI_COMM_WORLD) - 1;
  _localRank = _globalRank;
  // split master and slave
  Common::check(mpiSplit(MPI_COMM_WORLD, 1, _localRank, &_localComm));
  while (true) {
    _localRank = getRank(_localComm);
    if (!_localRank) {
      int signal;
      MPI_Status status;
      Common::check(mpiRecv(&signal, 1, MPI_INT, _globalMasterRank, TAG_MASTER_SIGNAL, MPI_COMM_WORLD, &status));
      if (SIGNAL_SPLIT == signal) {
        splitSlave();
      } else if (SIGNAL_JOB == signal) {
        treatJobSlave();         
      } else if (SIGNAL_TERMINATE == signal) {
        terminateSlave();
        break;
      }
    } else { 
      int signal;
      Common::check(mpiBcast(&signal, 1, MPI_INT, _localMasterRank, _localComm));
      if (SIGNAL_SPLIT == signal) {
        splitSlave();
      } else if (SIGNAL_JOB == signal) {
        treatJobSlave();         
      } else if (SIGNAL_TERMINATE == signal) {
        terminateSlave();
        break;
      }
    }
  }
  return 0;
}


SplitRanksAllocator::SplitRanksAllocator(unsigned int availableRanks,
    const string &outputDir):
  _totalRanks(availableRanks),
  _ranksInUse(0),
  _outputDir(outputDir)
{
  Common::makedir(Common::joinPaths(outputDir, "per_job_logs"));
  Common::makedir(Common::joinPaths(outputDir, "running_jobs"));
  MPI_Comm fakeComm;
  // split master and slave
  Common::check(mpiSplit(MPI_COMM_WORLD, 0, 0, &fakeComm));
  _slots.push(Slot(0, availableRanks - 1));
}


bool SplitRanksAllocator::ranksAvailable()
{
  return !_slots.empty();
}
 
bool SplitRanksAllocator::allRanksAvailable()
{
  return _ranksInUse == 0;
}
  
void SplitRanksAllocator::terminate()
{
  while (_slots.size()) {
    int signal = SIGNAL_TERMINATE;
    Common::check(mpiSend(&signal, 1, MPI_INT, _slots.front().startingRank, TAG_MASTER_SIGNAL, MPI_COMM_WORLD));
    _slots.pop();
  }
}

static void split(const SplitRanksAllocator::Slot &parent,
    SplitRanksAllocator::Slot &son1,
    SplitRanksAllocator::Slot &son2,
    int son1size)
{
  // send signal
  int signal = SIGNAL_SPLIT;
  Common::check(mpiSend(&signal, 1, MPI_INT, parent.startingRank, TAG_MASTER_SIGNAL, MPI_COMM_WORLD));
  Common::check(mpiSend(&son1size, 1, MPI_INT, parent.startingRank, TAG_SPLIT, MPI_COMM_WORLD));
  son1 = SplitRanksAllocator::Slot(parent.startingRank, (unsigned int)son1size);
  assert(int(parent.ranksNumber) >= son1size);
  son2 = SplitRanksAllocator::Slot(parent.startingRank + son1size, (unsigned int)(parent.ranksNumber - (unsigned int)son1size));
}


InstancePtr SplitRanksAllocator::allocateRanks(unsigned int requestedRanks, 
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
  _ranksInUse += slot.ranksNumber;
  
  InstancePtr instance(new SplitInstance(_outputDir,
    slot.startingRank,
    int(slot.ranksNumber),
    command));
  _rankToInstances[instance->getStartingRank()] = instance;
  return instance;
}
  
void SplitRanksAllocator::freeRanks(InstancePtr instance)
{
  _ranksInUse -= instance->getRanksNumber();
  auto splitInstance = static_pointer_cast<SplitInstance>(instance);
  _slots.push(Slot(instance->getStartingRank(), 
        (unsigned int)instance->getRanksNumber()));
}

vector<InstancePtr> SplitRanksAllocator::checkFinishedInstances()
{
  vector<InstancePtr> finished;
  Timer t;
  while (true) {
    MPI_Status status;
    int flag;
    Common::check(mpiIprobe(MPI_ANY_SOURCE, TAG_END_JOB, MPI_COMM_WORLD, &flag, &status)); 
    if (!flag) {
      break;
    }
    int source = status.MPI_SOURCE;
    int endJobMsg[MSG_SIZE_END_JOB];
    Common::check(mpiRecv(endJobMsg, MSG_SIZE_END_JOB, MPI_INT, source,  
        TAG_END_JOB, MPI_COMM_WORLD, &status));
    auto foundInstance = _rankToInstances.find(source);
    if (foundInstance == _rankToInstances.end()) {
      cerr << "Error: did not find the instance started on rank " << source << endl;
      break;
    }
    auto instance = static_pointer_cast<SplitInstance>(foundInstance->second);
    instance->setStartingElapsedMS(endJobMsg[1]);
    instance->setElapsedMs(endJobMsg[2]);
    if (endJobMsg[0]) {
      instance->onFailure(endJobMsg[0]);
    }
    finished.push_back(instance);
    //_rankToInstances.erase(source);

  }
  return finished;
}

void SplitRanksAllocator::preprocessCommand(CommandPtr cmd)
{
  auto ranksNumber = cmd->getRanksNumber();
  unsigned int newNumber = _totalRanks;
  while (newNumber > ranksNumber) {
    newNumber /= 2;
  }
  cmd->setRanksNumber(newNumber);
}

SplitInstance::SplitInstance(const string &outputDir, 
  int startingRank, 
  int ranksNumber,
  CommandPtr command):
  Instance(command, startingRank, ranksNumber, outputDir)
{
}


bool SplitInstance::execute(InstancePtr self)
{
  if (_ranksNumber == 0) {
    throw MPISchedulerException("Error in SplitInstance::execute: invalid number of ranks ", to_string(_ranksNumber));
  }
  int signal = SIGNAL_JOB;
  Common::check(mpiSend(&signal, 1, MPI_INT, self->getStartingRank(), TAG_MASTER_SIGNAL, MPI_COMM_WORLD));
  Common::check(mpiSend((char *)self->getId().c_str(), int(self->getId().size() + 1), MPI_CHAR, self->getStartingRank(), TAG_START_JOB, MPI_COMM_WORLD));
  return true;
}
  
} // namespace MPIScheduler

