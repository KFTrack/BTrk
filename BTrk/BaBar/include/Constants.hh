//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Constants.hh 475 2009-12-04 13:34:21Z stroili $
//
// Description:
//	The Constants class contains various constant values as
//      static members.
//
//      Other constants can be found in other classes.  For example,
//      AbsTrack/BField.hh contains various constants for converting
//      curvature to momentum, etc.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen 
//
// Copyright Information:
//	Copyright (C) 1995, 1996
//
//------------------------------------------------------------------------

#ifndef CONSTANTS_HH
#define CONSTANTS_HH

//-----------------
// BaBar Headers --
//-----------------
#include "BaBar/BaBar.hh"

class Constants {

  public:

  // Initialization of these is done in the .cc file. Although this
  // could be done here, it is not supported by the DEC C++ compiler
  // and is a recent change to the language reference.
  static const double pi;
  static const double twoPi;
  static const double radToDegrees;

  static const double c;

  // small value, replacing the BABAR_EPSILON from BaBar.hh
  static const double epsilon;

};

#endif

