//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BbrAngle.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//	Class BbrAngle
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Luca Lista	          	Original Author
//      Bob Jacobsen                    Added some implementation
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "BbrGeom/BbrAngle.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

const double BbrAngle::pi        = Constants::pi;
const double BbrAngle::twoPi     = Constants::twoPi;
const double BbrAngle::toDegrees = Constants::radToDegrees;

// The followings characters are used in DegString() 
// for printout in degrees
#ifdef HPUX
const HepString BbrAngle::degChar = "°";
#else
const HepString BbrAngle::degChar = "^";
#endif

const HepString BbrAngle::deg1Char =  "'";
const HepString BbrAngle::deg2Char = "\"";

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

HepString BbrAngle::degString() const
{
  float tmp;
  int deg, deg_, deg__;
  HepString sign = "";
  
  if ((tmp = this->deg()) < 0) 
  { sign = "-"; tmp = -tmp; };
  deg   = int(tmp);
  deg_  = int(tmp = 60*(tmp - deg));
  deg__ = int(60*(tmp - deg_));
  
  return 
  (
    sign + 
    HepString(deg)+degChar+
    HepString(deg_)+deg1Char+
    HepString(deg__)+deg2Char
  );
}


