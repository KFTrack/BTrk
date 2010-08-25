//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
// $Id: MinMax.cc,v 1.1 2010/06/22 16:05:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:05:18 $
//
// Original author Rob Kutschke

#include "BaBar/BaBar.hh"
#include "mu2eFast/MinMax.hh"

#include <limits>
#include <cmath>

using namespace std;


MinMax::MinMax():
  _min(   std::numeric_limits<double>::max()),
  _max(  -std::numeric_limits<double>::max()),
  _small( std::numeric_limits<double>::max())
{}

MinMax::MinMax(double x):
  _min(   std::numeric_limits<double>::max()),
  _max(  -std::numeric_limits<double>::max()),
  _small( std::numeric_limits<double>::max()){
  compare(x);
}


void MinMax::compare(double x){
  _min   = (x < _min ) ? x : _min;
  _max   = (x > _max ) ? x : _max;
  _small = ( abs(x) < _small ) ? abs(x) : _small;
}
