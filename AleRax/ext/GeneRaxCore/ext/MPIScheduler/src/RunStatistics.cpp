
#include "RunStatistics.hpp"

#include "SVGDrawer.hpp"
#include "Logger.hpp"

using namespace std;

namespace MPIScheduler {

RunStatistics::RunStatistics(const InstancesHistoric &historic,
    Time begin,
    Time end,
    unsigned int availableRanks,
    Logger &masterLogger):
  _historic(historic),
  _begin(begin),
  _end(end),
  _availableRanks(availableRanks),
  _lbRatio(1.0),
  _masterLogger(masterLogger)
{
  assert(_availableRanks > 0);
}

void RunStatistics::printGeneralStatistics()
{
  long totalElapsedTime = Common::getElapsedMs(_begin, _end);
  long cumulatedTime = 0;
  for (auto instance: _historic) {
    cumulatedTime += instance->getElapsedMs() * instance->getRanksNumber();
  }
  _lbRatio = double(cumulatedTime) / double(_availableRanks * totalElapsedTime);
  
  _masterLogger.getCout() << "Finished running commands. Total elasped time: ";
  _masterLogger.getCout() << totalElapsedTime / 1000  << "s" << endl;
  _masterLogger.getCout() << "Load balance ratio: " << _lbRatio << endl;
}


void RunStatistics::exportSVG(const string &svgfile) 
{
  Timer t;
  _masterLogger.getCout() << "Saving svg output in " << svgfile << endl;
  auto totalWidth = _availableRanks + 1;
  _masterLogger.getCout() << "total width " << totalWidth << endl;
  auto totalHeight = Common::getElapsedMs(_begin, _end);
  string caption = "t = " + to_string(totalHeight / 1000) + "s";
  caption += ", lb = " + to_string(_lbRatio);
  SVGDrawer svg(svgfile, double(totalWidth), double(totalHeight));
  svg.writeHeader(caption);
  for (auto instance: _historic) {
    instance->writeSVGStatistics(svg, _begin);
  }
  svg.writeFooter();
  _masterLogger.getCout() << "Time spent writting svg: " << t.getElapsedMs() / 1000 << "s" << endl;
}

} // namespace MPIScheduler



