#include "Instance.hpp"

#include "SVGDrawer.hpp"
#include "Common.hpp"

using namespace std;

namespace MPIScheduler {

bool Instance::_jobFailureFatal = false;

Instance::Instance(CommandPtr command,
    int startingRank,
    int ranksNumber,
    const string &outputDir):
  _startingElapsedMS(0),
  _outputDir(outputDir),
  _command(command),
  _startingRank(startingRank),
  _ranksNumber(ranksNumber),
  _endTime(Common::getTime()),
  _elapsed(0)
{
  assert(_ranksNumber > 0);
}

void Instance::writeSVGStatistics(SVGDrawer &drawer, const Time &)
{
  drawer.writeSquare(getStartingRank(),
    _startingElapsedMS,
    getRanksNumber(),
    getElapsedMs());

}
void Instance::onFinished()
{
  _endTime = Common::getTime();
}

void Instance::onFailure(int errorCode)
{
  string outputFile = Common::joinPaths(_outputDir, "failed_commands.txt");
  bool exists = ifstream(outputFile.c_str()).good();
  ofstream os(outputFile, fstream::out | fstream::app);
  if (exists)
    os << endl;
  os << _command->getId();
  cerr << "Warning, command " << _command->getId() << " failed ";
  if (errorCode) {
    cerr << " with code " << errorCode;
  }
  cerr << endl;
  if (_jobFailureFatal) {
    cerr << "Job failures are fatal, aborting" << endl;
    exit(242);
  }
}

}
