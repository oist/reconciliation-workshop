#pragma once

#include <memory>
#include <string>
#include <vector>

#ifdef MPISCHEDULER_STATIC_SCHEDULED_MAIN
extern "C" int static_scheduled_main(int argc, char** argv, void* comm);
#endif

typedef int (*mainFct)(int,char**,void*);  

class DynamicLibrary {
public:
  DynamicLibrary() = default;
  ~DynamicLibrary();
  static std::shared_ptr<DynamicLibrary> getLibrary(const std::string &libraryPath);
  int run(const std::string &logsFile,
      const std::string &errFile,
      const std::vector<std::string> &args,
      void *comm); 
private:
  bool setPath(const std::string &libraryPath);
  void *_handle;
  mainFct _raxmlMain;
};




