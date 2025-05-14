#pragma once

#include <string>
#include <fstream>

namespace MPIScheduler {

class SVGDrawer {
public:
  SVGDrawer(const std::string &filepath,
      double maxXValue,
      double maxYValue);
  ~SVGDrawer();
  void writeSquare(double x, double y, double w, double h, const char *color = 0);
  void writeSquareAbsolute(double x, double y, double w, double h, const char *color = 0);
  static std::string getRandomHex();
  void writeHeader(const std::string &caption);
  void writeFooter();

private:
  void writeCaption(const std::string &text);
  std::ofstream _os;
  double _width;
  double _height;
  double _ratioWidth;
  double _ratioHeight;
};

} // namespace MPIScheduler

