#pragma once

#include <iostream>
#include <fstream>

namespace MPIScheduler {

class Logger {
public:
  Logger(): _os(0) {

  }

  ~Logger() {
    delete _os;
  }

  void redirectLogs(const std::string &file) {
    getCout() << "Redirecting logs to " << file << std::endl;
    delete _os;
    _os = new std::ofstream(file);
  }

  std::ostream &getCout() {
    if (_os) {
      return *_os;
    }
    return std::cout;
  }

private:
  std::ofstream *_os;
};










}
