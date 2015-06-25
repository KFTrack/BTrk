#ifndef GNUPLOT_HH
#define GNUPLOT_HH

#include <fstream>
#include "BTrk/BbrGeom/HepPoint.h"
#include <vector>
//
//  Structure for GNU plotting, used in GnuPlot functions.
//  This is used with std vectors
//
struct GnuPlot{
  std::ofstream* output;
  bool active;
  std::vector<HepPoint>& stlPoints;
  GnuPlot(std::ofstream& stream,bool act,
	  std::vector<HepPoint>& pnts):
    output(&stream),active(act),stlPoints(pnts)
  {;}
};

#endif
