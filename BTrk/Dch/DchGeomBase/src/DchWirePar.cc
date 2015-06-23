//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchWirePar.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchWirePar
//      class used to store basic wire information in the database
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
#include "DchGeomBase/DchWirePar.hh"

//---------------
// C++ Headers --
//---------------
#include <math.h>
using std::endl;
using std::ostream;

//----------------
// Constructors --
//----------------
DchWirePar::DchWirePar()
  : _ID(-1), _radRear(0.), _phiRear(0.), _zRear(0.), _radForw(0.),
    _phiForw(0.), _zForw(0.), _sag(0.), _diameter(0.)
{
}

DchWirePar::DchWirePar( int id )
  : _ID(id), _radRear(0.), _phiRear(0.), _zRear(0.), _radForw(0.),
    _phiForw(0.), _zForw(0.), _sag(0.), _diameter(0.)
{
}

DchWirePar::DchWirePar( int id, double radR, double phiR, double zR, 
			double radF, double phiF, double zF, double sag, 
			double diameter ) 
  : _ID(id), _radRear(radR), _phiRear(phiR), _zRear(zR), _radForw(radF),
    _phiForw(phiF), _zForw(zF), _sag(sag), _diameter(diameter)
{
}

DchWirePar::DchWirePar( int id, double radR, double phiR, double zR, 
			double length, double twist, double sag, 
			double diameter )
  : _ID(id), _radRear(radR), _phiRear(phiR), _zRear(zR), _radForw(radR),
    _sag(sag) , _diameter(diameter)
{
  _phiForw = _phiRear+twist;
  _zForw = _zRear+length;
}

DchWirePar::DchWirePar( const DchWirePar& w ) 
{
  _ID = w.getID();
  _radRear = w.getRadRear();
  _phiRear = w.getPhiRear();
  _zRear = w.getZRear();
  _radForw = w.getRadForw();
  _phiForw = w.getPhiForw();
  _zForw = w.getZForw();
  _sag = w.getSag();
  _diameter = w._diameter;
}

//--------------
// Destructor --
//--------------
DchWirePar::~DchWirePar() 
{}

//-------------
// Operators --
//-------------
DchWirePar& DchWirePar::operator= ( const DchWirePar& w ) 
{
  if (&w == this) return *this;

  _ID = w.getID();
  _radRear = w.getRadRear();
  _phiRear = w.getPhiRear();
  _zRear = w.getZRear();
  _radForw = w.getRadForw();
  _phiForw = w.getPhiForw();
  _zForw = w.getZForw();
  _sag = w.getSag();
  _diameter = w.getDiameter();

  return *this;
}

int DchWirePar::operator== ( const DchWirePar& w ) const
{ 
  // very simple equality test

  if ( _ID == w.getID() ) return 1; 
 
  return 0;
} 
 
    
//-------------
// Selectors --
//-------------
double DchWirePar::stereo() const
{
  return 2.*getRadRear()*sin(twist()*0.5)/zLength(); 
}

void DchWirePar::printAll(ostream& o) const {
  o << _ID << " -> " << _radRear << " " << _phiRear << " " << _zRear 
    << " " << _radForw << " " << _phiForw << " " << _zForw
    << "\n\t" << _sag << " " << _diameter << endl;
}

ostream&  operator<<(ostream& o, DchWirePar& cell) {
  cell.printAll(o);
  return o;
}
