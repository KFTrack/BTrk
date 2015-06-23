//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFit.cc,v 1.9 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkFit.hh"
#include <iostream>
using std::ostream;

//------------------------------------------------------------------------
TrkFit::~TrkFit() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkFit::TrkFit() {
//------------------------------------------------------------------------
}

void 
TrkFit::printType(ostream& ostr) const 
{
  ostr << "Particle type: " << particleType().name();
}

// default implementation of validFlightLength
bool
TrkFit::validFlightLength(double fltl,double tolerance) const {
  return fltl+tolerance >= startValidRange() &&
    fltl-tolerance <= endValidRange();
}
