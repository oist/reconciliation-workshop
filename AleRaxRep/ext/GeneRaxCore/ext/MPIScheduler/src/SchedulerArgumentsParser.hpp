#pragma once
#include <string>


namespace MPIScheduler {

class SchedulerArgumentsParser {
public:
  SchedulerArgumentsParser(int argc, char** argv):
    commandsFilename()
  {
    if (argc < 6) {
      std::cerr << "INVALID SYNTAX" << std::endl;
      print_help();
      exit(1);
    }
    int i = 1;
    implem = std::string(argv[i++]);
    coresNumber = atoi(argv[i++]);
    library = std::string(argv[i++]);
    commandsFilename = std::string(argv[i++]);
    outputDir = std::string(argv[i++]);
    jobFailureFatal = 0;
    for (;i < argc; ++i) {
      auto currentArg = std::string(argv[i]);
      if (currentArg == std::string("--jobs-failure-fatal")) {
        jobFailureFatal = true;
      } else if (currentArg == std::string("--threads-arg")) {
        threadsArg = std::string(argv[++i]);
      } else if (currentArg == std::string("--logs")) {
        outputLogs = std::string(argv[++i]);
      }
    }
  }
  
  void print_help() 
  {
    std::cout << "Syntax:" << std::endl;
    std::cout << "mpi-scheduler implementation cores_number library_or_exec command_file output_dir [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --jobs-failure-fatal\t Stop scheduler when at least one job fails" << std::endl;
    std::cout << "  --threads-arg <string>\t Only makes sense with fork implem. Set the command argument to use to set the number of threads for the scheduled executable" << std::endl;
    std::cout << "  --logs <string>\t Path to the file in which the scheduler should write logs" << std::endl;
  }

  std::string implem;
  std::string library;
  std::string commandsFilename;
  std::string outputDir;
  int coresNumber;
  int jobFailureFatal;
  std::string threadsArg;
  std::string outputLogs;
};

} 
