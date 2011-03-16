//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchFWire.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchFWire.
//      for the moment inherits completely from DchSWire
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

#ifndef DCHFWIRE_HH
#define DCHFWIRE_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <iostream>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/HepPoint.h"
#include "TrajGeom/TrkLineTraj.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchFWire {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  //  Constructors
  DchFWire(HepPoint& rP, HepPoint& fP, double sag=0);

  //  Destructor
  virtual ~DchFWire( );

  //  Operators
  virtual int operator==( const DchFWire& ) const;

  //  Selectors (const)
  const TrkLineTraj* getTraj(void)   const { return &_traj; }
  const HepPoint* getRearPoint(void) const { return &_rear; }
  const HepPoint* getForwPoint(void) const { return &_forward; }

  //  return point of wire at z (in the BaBar reference frame)
  double gotoLocal(double z) { return z-_zOffset; }
  double gotoGlob(double zLoc) { return zLoc+_zOffset; }
  HepPoint getPoint(double z);
  Hep3Vector getDirection(double z)
                        { return getTraj()->direction(gotoLocal(z)).unit(); }
  double rEnd(void) const { return sqrt(_rear.x()*_rear.x()+
					_rear.y()*_rear.y()); }
  double dPhiz(void) const { return (_forward.phi() - _rear.phi())*0.5; }
  double zLength(void) const { return getForwPoint()->z()-
				       getRearPoint()->z(); }
  double stereo(void) const
                        { return 2.*rEnd()*sin(dPhiz())/zLength(); }

  void print(std::ostream& o=std::cout) const ;
  void printAll(std::ostream& o=std::cout) const ;

private:

  //  Data members
  //  Data members: THESE ARE IN THE GLOBAL (BaBar) REFERENCE SYSTEM
  TrkLineTraj     _traj;  // wire trajectory
  HepPoint        _rear;  // wire position on rear end-plate
  HepPoint     _forward;  // wire position on forward end-plate
  double           _sag;  // wire sagitta
  double       _zOffset;  // z offset

  //  copy constructor & != operator
  DchFWire( const DchFWire& );
  DchFWire& operator=( const DchFWire& );
  int operator!=( const DchFWire& ) const;

};

std::ostream& operator << (std::ostream& o, const DchFWire&);

#endif // DCHFWIRE_HH
