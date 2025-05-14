
#include <iostream>
#include <string>
#include <chrono>
#include <ctime>  
#include <cassert>
#include "RanksAllocator.hpp"
#include "ParallelImplementation.hpp"
#include "Command.hpp"
#include "Common.hpp"
#include "CommandsRunner.hpp"
#include "RunStatistics.hpp"
#include "SchedulerArgumentsParser.hpp"

using namespace std;

namespace MPIScheduler {

void printStart(int argc, char **argv)
{
  auto start = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(start);
  std::cout << "Program started at " << std::ctime(&start_time) << std::endl;
  std::cout << "MPIScheduler was called as follow:" << std::endl;
  for (int a = 0; a < argc; ++a) {
    std::cout << argv[a] << " ";
  }
  std::cout << std::endl;
}

  
static int main_scheduler(int argc, char **argv, void* comm)
{
  // Init
  SchedulerArgumentsParser arg(argc, argv);
  ParallelImplementation implem(arg.implem);
  if (!implem.isValid()) {
    cerr << "Invalid scheduler implementation: " << arg.implem << endl;
    return 1;
  }
  implem.initParallelContext(argc, argv, comm, arg.coresNumber);
  if (!implem.getRank()) {
    printStart(argc, argv);
  }
  if (implem.slavesToStart()) {
    implem.startSlaves(argc, argv);
    if (implem.getRank() != int(implem.getRanksNumber()) - 1) {
      implem.closeParallelContext();
      return 0;
    }
  }
  Logger masterLogger;
  if (arg.outputLogs.size() != 0) {
    masterLogger.redirectLogs(arg.outputLogs);
  }
  Time begin = Common::getTime();
  CommandsContainer commands(arg.commandsFilename,
      implem.addFakeExecutableName());
  auto allocator = implem.getRanksAllocator(arg, implem.getRanksNumber());
  for (auto command: commands.getCommands()) {
    allocator->preprocessCommand(command);
  }
  
  // Run 
  CommandsRunner runner(commands, allocator, arg.outputDir, arg.jobFailureFatal, masterLogger);
  runner.run(implem.isMPI());
  masterLogger.getCout() << "end of run" << endl;
  // End
  Time end = Common::getTime();
  assert(implem.getRanksNumber() > 0);
  RunStatistics statistics(runner.getHistoric(), begin, end, (unsigned int)(implem.getRanksNumber()), masterLogger);
  statistics.printGeneralStatistics();
  if (runner.getHistoric().size()) {
    statistics.exportSVG(Common::getIncrementalLogFile(arg.outputDir, "statistics", "svg"));
  }
  allocator->terminate();
  implem.closeParallelContext();
  masterLogger.getCout() << "End of Multiraxml run" << endl;
  return 0;
}

} // namespace MPIScheduler



#ifdef MPISCHEDULER_BUILD_AS_LIBRARY

extern "C" int mpi_scheduler_main(int argc, char** argv, void* comm)
{
  assert(argv);
  assert(argv);
  int res =  MPIScheduler::main_scheduler(argc, argv, comm);
#ifdef WITH_MPI
  if (comm) { 
    MPI_Barrier(*((MPI_Comm*)comm)); 
  }
#endif
  return res;
}

#else

int main(int argc, char** argv) 
{
  exit(MPIScheduler::main_scheduler(argc, argv, 0));
}


#endif
