#include "Common.hpp"
#include<iostream>

#include <stdio.h>
#include <ftw.h>
#include <unistd.h>

using namespace std;

namespace MPIScheduler {

string Common::getIncrementalLogFile(const string &path, 
      const string &name,
      const string &extension)
{
  string file = joinPaths(path, name + "." + extension);
  int index = 1;
  while (ifstream(file).good()) {
    file = joinPaths(path, name + "_" + to_string(index) + "." + extension);
    index++;
  }
  return file;
}

int systemCall(const string &command, const string &outputFile,
    bool threadSafe)
{

  if (threadSafe) {
    return system((command + " > " + outputFile).c_str());
  } else {
    int result = 0;
    FILE *ptr, *file;
    file = fopen(outputFile.c_str(), "w");
    if (!file) {
      cerr << "[MPIScheduler error] Cannot open output file " << outputFile << endl;
      return 0;
    }
    //cout << command.c_str() << endl;
    if ((ptr = popen(command.c_str(), "r")) != NULL) {
      char buf[BUFSIZ];
      while (fgets(buf, BUFSIZ, ptr) != NULL) {
        fprintf(file, "%s", buf);
      }
      result = pclose(ptr);
    }
    fclose(file);
    return WEXITSTATUS(result);
  }
}

} // namespace MPIScheduler

