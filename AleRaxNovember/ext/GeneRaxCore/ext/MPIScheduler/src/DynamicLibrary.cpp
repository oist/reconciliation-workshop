#include "DynamicLibrary.hpp"
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <cassert>
using namespace std;

  
shared_ptr<DynamicLibrary> DynamicLibrary::getLibrary(const string &libraryPath)
{
  auto library = make_shared<DynamicLibrary>();
  if (!library->setPath(libraryPath)) {
    return 0;
  }
  return library;
}
  
int DynamicLibrary::run(const string &logsFile,
      const string &errFile,
      const vector<string> &args,
      void *comm) 
{
  std::ofstream out(logsFile);
  std::streambuf *coutbuf = std::cout.rdbuf(); 
  std::cout.rdbuf(out.rdbuf()); 
  std::ofstream err(errFile);
  std::streambuf *cerrbuf = std::cerr.rdbuf(); 
  std::cerr.rdbuf(err.rdbuf());
  auto argc = args.size(); 
  auto argv = new char*[argc];
  for (unsigned int i = 0; i < argc; ++i) {
    argv[i] = const_cast<char *>(args[i].c_str());
    cout << argv[i] << " ";
  }
  cout << endl;
  int res = -1;
  try {
    res =_raxmlMain((int)argc, argv, (void*)&comm);
  } catch (exception &e) {
    cerr << "Catched exception: " << e.what() << endl;
  }                    
  delete[] argv;
  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  out.close();
  err.close();
  return res;
}

bool DynamicLibrary::setPath(const string &libraryPath)
{
  if (libraryPath == "--static_scheduled_main") {
#ifdef MPISCHEDULER_STATIC_SCHEDULED_MAIN
    _raxmlMain = (mainFct) static_scheduled_main;
    assert(_raxmlMain);
    return true;
#else
    cerr << "To use static_scheduled_main, please compile with MPISCHEDULER_STATIC_SCHEDULED_MAIN cmake variable on" << endl;
    return false;;
#endif
  }
  _handle = dlopen(libraryPath.c_str(), RTLD_LAZY);
  if (!_handle) {
    cerr << "Cannot open shared library " << libraryPath << endl;
    cerr << "Error: " << dlerror() << endl;
    return false;
  }
  _raxmlMain = (mainFct) dlsym(_handle, "dll_main");
  const char *dlsym_error = dlerror();
  if (dlsym_error) {
    cerr << "Cannot load symbole dll_main " << dlsym_error << endl;
    dlclose(_handle);
    _handle = nullptr;
    return false;
  }
  return true;
}

DynamicLibrary::~DynamicLibrary()
{
  if (_handle)
    dlclose(_handle);
}



