//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSagTraj.cc 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchSagTraj
//      
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
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchSagTraj.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
#include "BbrGeom/Trajectory.hh"

//----------------
// Constructors --
//----------------
DchSagTraj::DchSagTraj( const double sag, const HepPoint& start, 
			const HepPoint& stop )
  : TrkGeomTraj(0.0,start.distanceTo(stop)), _sag(sag), _start(start), _stop(stop)
{
  _direction = stop - start;
  _length = _direction.mag();
  assert( _length != 0 );

  _direction.setMag(1.0);
  _a = _sag*4./(_length*_length);
  _b = -_a*_length;

}

DchSagTraj::DchSagTraj( const DchSagTraj& other )
  : TrkGeomTraj(0.0,other._start.distanceTo(other._stop)), 
    _sag(other._sag), _a(other._a), _b(other._b), _length(other._length),
    _start(other._start), _stop(other._stop), 
    _direction(other._direction)
{}

DchSagTraj::~DchSagTraj(){;}
//-------------
// Methods   --
//-------------
DchSagTraj*
DchSagTraj::clone() const 
{
  return new DchSagTraj(*this);
}
    
//-------------
// Operators --
//-------------
//-----------------------------------------------------------------------------
DchSagTraj& 
DchSagTraj::operator = (const DchSagTraj& other) 
{
  if(&other != this){
    for(int iend=0;iend<2;iend++)
      flightrange[iend] = other.flightrange[iend];
    _start = other._start;
    _stop = other._stop;
    _sag = other._sag;
    _a = other._a;
    _b = other._b;
    _length = other._length;
    _direction = other._direction;
  }
  return *this;
}

Hep3Vector 
DchSagTraj::deviation( double flightlen ) const 
{
  // only correction in y
//   Hep3Vector deviation(0., (_a*flightlen+_b)*flightlen, 0.);
//   Hep3Vector deviaH(0.,cosh((flightlen-_length/2.)/_a_H)+_b_H,0.);
//   return deviation;
  return Hep3Vector(0., (_a*flightlen+_b)*flightlen, 0.);
}

HepPoint
DchSagTraj::position(double flightlen) const 
{
  static HepPoint tmppos;
  tmppos = _start;
  tmppos += _direction*flightlen;
  tmppos.setY(tmppos.y()+(_a*flightlen+_b)*flightlen);
  return tmppos;
}

Hep3Vector 
DchSagTraj::direction( double flightlen ) const {
  if ( flightlen <= 0. ) return _direction;
//    Hep3Vector dir =  _direction*flightlen + delDirect(flightlen);
  static Hep3Vector tmpdir;
  tmpdir =  _direction*flightlen;
//   register double newy = tmpdir.y() + 2.*_a*flightlen+_b;
//   tmpdir.setY(newy);
  tmpdir.setY(tmpdir.y() + 2.*_a*flightlen+_b);
//   tmpdir += delDirect(flightlen);
  tmpdir.setMag(1.0);
  return tmpdir;
}
  
Hep3Vector 
DchSagTraj::delDirect( double flightlen ) const 
{
  return Hep3Vector(0., 2.*_a, 0.);
}
 
void  
DchSagTraj::getInfo(double flightlen, HepPoint& pos, Hep3Vector& dir) const 
{
  // Written using +=, etc to avoid temporaries
  pos = _start;
  pos += _direction*flightlen;
  
  dir = _direction;
  if ( flightlen > 0. ) {
    pos.setY(pos.y() + deltaY(flightlen));
    dir.setY(dir.y() + 2.*_a*flightlen+_b);
    dir.setMag(1.0);
  }
}

void  
DchSagTraj::getInfo( double flightlen, HepPoint& pos, Hep3Vector& dir, 
		     Hep3Vector& delDir ) const 
{
  pos = _start;
  pos += _direction*flightlen ;
  pos.setY( pos.y() + (_a*flightlen+_b)*flightlen) ;

  dir = _direction;
  dir.setY( dir.y() + 2*_a*flightlen+_b) ;
  // Note: `dir' is on purpose not normalized (WDH, Jan 2003)

  delDir.setX(0.);
  delDir.setY(2.*_a);
  delDir.setZ(0.);
}

double
DchSagTraj::curvature( double f ) const
{
  return _sag;
}

double 
DchSagTraj::distTo1stError(double flightlen, double tol, int pathDir) const 
{
  double dtmp = pathDir*2.*_a*flightlen+_b;

  return dtmp==0. ? 9999.e4 : fabs(tol/dtmp);
}
  
double 
DchSagTraj::distTo2ndError(double s, double tol, int pathDir) const
{ 
  return 999.e4 ;
  //return _a==0. ? 999.e4 : tol/(2.*_a); 
}

// Support Visitor pattern (see TrkGeomTraj.hh)
void 
DchSagTraj::accept(TrkGeomTrajVisitor& visitor) const 
{
  ErrMsg(error)<<"accept visitor NOT implemented yet"<<endmsg;
}

