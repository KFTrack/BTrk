//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchFWire.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchFWire
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchFWire.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <assert.h>
using std::endl;
using std::ostream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] =
    "$Id: DchFWire.cc 123 2010-04-29 14:41:45Z stroili $";

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchFWire::DchFWire(HepPoint& rP, HepPoint& fP, double sag) :
  _traj(rP, fP), _rear(rP), _forward(fP), _sag(sag), _zOffset(_rear.z())
{
}

//--------------
// Destructor --
//--------------

DchFWire::~DchFWire()
{
  ;
}
//-------------
// Operators --
//-------------
int
DchFWire::operator ==(const DchFWire& fw) const
{
  if (&_traj != fw.getTraj()) return 0;
  return 1;
}

//-------------
// Selectors --
//-------------
HepPoint
DchFWire::getPoint(double z)
{
  // this is valid in the assumption that the wire trajectory is a line
  // will have to change this when it will be a parabola
  double zLoc = z - _zOffset;
  register double zdir = _traj.direction(zLoc).unit().z();
  assert(zdir != 0.);
  register double flt = zLoc / zdir;

  return _traj.position(flt);
}

void
DchFWire::print(ostream& o) const
{
  o << " DchFWire rear point: " << *getRearPoint() << endl;
}

void
DchFWire::printAll(ostream& o) const
{
  o << " DchFWire rear point: " << *getRearPoint() << " forw point: "
      << *getForwPoint() << endl;
}

ostream&
operator <<(ostream& o, const DchFWire& w)
{
  w.print(o);
  return o;
}
