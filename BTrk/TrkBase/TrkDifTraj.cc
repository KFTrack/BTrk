//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifTraj.cc,v 1.8 2000/02/04 22:05:03 jalbert Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/difAlgebra/DifVector.hh"

//Constructor
TrkDifTraj::TrkDifTraj(const double lowlim,const double hilim) :
  Trajectory(lowlim, hilim) {
}
 
TrkDifTraj::~TrkDifTraj() {}

void  
TrkDifTraj::getDFInfo2(double fltLen, DifPoint& pos, DifVector& direction) const {
  // Slow default implementation.  Override in subclasses where speed matters
  DifVector dummy;
  getDFInfo(fltLen, pos, direction, dummy);
}




