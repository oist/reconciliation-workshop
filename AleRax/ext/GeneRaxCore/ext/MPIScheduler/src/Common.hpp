
#pragma once

#ifdef WITH_MPI
#include <cstdio>
#include <mpi.h>
#endif

#include <fstream>
#include <chrono>
#include <thread>
#include <cstdio>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <signal.h>  
#include <iostream>
#include <string>
#include <cassert>

namespace MPIScheduler {

class MPISchedulerException: public std::exception {
public:
  explicit MPISchedulerException(const std::string &s): msg_(s) {}
  MPISchedulerException(const std::string &s1, 
      const std::string &s2): msg_(s1 + s2) {}
  const char* what() const noexcept {return msg_.c_str();}
private:
  std::string msg_;
};

using Time = std::chrono::time_point<std::chrono::system_clock>;


class Common {
public:
  // todobenoit not portable
  // todobenoit not portable
  static void makedir(const std::string &name) {
    mkdir(name.c_str(), 0755);
  }
  
  static std::string joinPaths(const std::string &path1, const std::string &path2) {
    return path1 + "/" + path2;
  }

  static std::string joinPaths(const std::string &path1, 
      const std::string &path2,
      const std::string &path3) {
    return joinPaths(joinPaths(path1, path2), path3);
  }

  static Time getTime() {
    return std::chrono::system_clock::now();
  }

  static std::string getIncrementalLogFile(const std::string &path, 
      const std::string &name,
      const std::string &extension);

  static long getElapsedMs(Time begin, Time end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>
      (end-begin).count();
  }

#ifdef WITH_MPI
  static void check(int mpiError) {
    if (mpiError != MPI_SUCCESS) {
      std::cout << "MPI ERROR !!!!" << std::endl;
      throw MPISchedulerException("MPI error !");
    }
  }
#endif

  static int getPid() {
    return getpid();
  }
  
  static std::string getHost() {
#ifdef HOST_NAME_MAX // not defined on OSX
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    return std::string(hostname);
#else
    return std::string();
#endif
  }


};

int systemCall(const std::string &command, const std::string &outputFile, bool threadSafe = false);


class Timer {
public:
  Timer() {
    reset();
  }

  long getElapsedMs() const {
    auto end = Common::getTime();
    return Common::getElapsedMs(_start, end);
  }

  void reset() {
    _start = std::chrono::system_clock::now();
  }
private:
  Time _start;
};

} // namespace MPIScheduler

