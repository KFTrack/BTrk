//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:$
//
// Description:
//	The Constants class contains various constant values as
//      static members.
//
// Environment:
//      This version is a stripped down version for the port to Mu2e.
//
// Author List:
//      Bob Jacobsen 
//      Rob Kutschke, Dave Brown (port to Mu2e).
//
// Copyright Information:
//	Copyright (C) 1995, 1996
//
//------------------------------------------------------------------------
#include "BaBar/Constants.hh"
#include "CLHEP/Units/PhysicalConstants.h"

const double Constants::pi           = CLHEP::pi;
const double Constants::twoPi        = CLHEP::twopi;
const double Constants::radToDegrees = 180./Constants::pi;

const double Constants::c  = CLHEP::c_light/CLHEP::centimeter*CLHEP::second;

