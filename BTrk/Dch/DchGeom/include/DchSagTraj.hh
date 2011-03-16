//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSagTraj.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchSagTraj.
//      the flight length at the moment is the flight length of a linear 
//      trajectory, the sag is approximated with a parabola (for the
//      moment, maybe at a later time we can replace with a more accurate
//      description, cosh for example)
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1998	INFN & Padova University
//
//------------------------------------------------------------------------

#ifndef DCHSAGTRAJ_HH
#define DCHSAGTRAJ_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "BbrGeom/TrkGeomTraj.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class Trajectory;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchSagTraj : public TrkGeomTraj {

public:

  // Constructors
  DchSagTraj(const double sag, const HepPoint& point1, const HepPoint& point2);
  DchSagTraj(const DchSagTraj& traj);
 
  DchSagTraj* clone() const;
  // Destructor
  virtual ~DchSagTraj( );

  // Operators
  DchSagTraj& operator=(const DchSagTraj&);

  double sag( void ) const { return _sag; }

// needed implementations for intersection with a Surface
  HepPoint position( double ) const;
  Hep3Vector direction( double ) const;
  double curvature( double f = 0. ) const;
  Hep3Vector delDirect( double ) const;
  const Hep3Vector& rawDirection( void ) const { return _direction; }
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction) const;
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction, 
			Hep3Vector& delDirect) const;

  virtual double distTo1stError(double s, double tol, int pathDir) const;
  virtual double distTo2ndError(double s, double tol, int pathDir) const;
  double deltaY(double fltLen) const { return (_a*fltLen+_b)*fltLen; }
  // Support Visitor pattern (see TrkGeomTraj.hh)
  void accept(TrkGeomTrajVisitor& visitor) const;

//    virtual int operator==( const DchSagTraj& ) const;
//            int operator!=( const DchSagTraj& ) const;

private:

  // Data members
  double _sag;
  double _a;
  double _b;
  double _length;
  HepPoint _start; // where the trajectory starts
  HepPoint _stop; // where the trajectory stops
  Hep3Vector _direction; // direction (unit) vector for 
                         // null sag (straight line)

  Hep3Vector deviation(double) const; // displacement from the line trajectory
                                      // at a given flightlength
  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
//   DchSagTraj( const DchSagTraj& );       // Copy Constructor
//   DchSagTraj&       operator= ( const DchSagTraj& );  // Assignment op

};

#endif // DCHSAGTRAJ_HH
