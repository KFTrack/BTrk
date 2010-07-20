#ifndef MinMax_HH
#define MinMax_HH

//
// Find minimum, maximum and smallest in magnitude of a set of numbers
// presented one at a time.
//
//
// Original author Rob Kutschke

#include <limits>
#include <cmath>

class MinMax{
  
public:

  MinMax();
  MinMax(double x);
  void compare(double x);

  // Return limiting values.
  double min() const {return _min;}
  double max() const {return _max;}
  double smallest() const {return _small;}

private:

  // Limiting values of the numbers compared so far.
  double _min;
  double _max;
  double _small;


};

#endif
