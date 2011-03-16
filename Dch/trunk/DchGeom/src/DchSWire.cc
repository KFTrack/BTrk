//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSWire.cc 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchSWire
//      Do not use this for DchSWired class (foo<T>).  use DchSWireDchSWire.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchSWire.hh"

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BbrGeom/BbrAngle.hh"
#include "CLHEP/Geometry/Transformation.h"
#include "DchGeom/DchLayer.hh"
#include "DchGeomBase/DchWirePar.hh"
#include "DetectorModel/DetAlignElem.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ostream;


//----------------
// Constructors --
//----------------
DchSWire::DchSWire(const HepPoint& rP, const HepPoint& fP, double sag)
  :  _traj(sag,rP,fP), _rear(rP), _forward(fP), _sag(sag),
     _phiend(rP.phi()), _id(0)
{
  _rend = sqrt(xRear()*xRear() + yRear()*yRear());
  BbrAngle fPphi(fP.phi());
  BbrAngle rPphi(rP.phi());
  _twist = BbrAngle(fP.phi() - rP.phi()) * 0.5;
}

//--------------
// Destructor --
//--------------
DchSWire::~DchSWire( ) 
{;}

//-------------
// Modifiers --
//-------------
void
DchSWire::wireAlign( const DchWirePar& corr ) 
{
  double newr = _rend+corr.getRadRear();
  double newphi = _phiend+corr.getPhiRear();
  double xcorr = newr*cos(newphi);
  double ycorr = newr*sin(newphi); 
  double zcorr = corr.getZRear();

  HepPoint newRear(xcorr, ycorr, zRear()+zcorr);

  _rear = newRear;

  newr = sqrt(xForw()*xForw()+yForw()*yForw())+corr.getRadForw();
  newphi = _forward.phi()+corr.getPhiForw();
  xcorr = newr*cos(newphi); 
  ycorr = newr*sin(newphi);  
  zcorr = corr.getZForw(); 
  HepPoint newForw(xcorr, ycorr, zForw()+zcorr); 
 
  _forward = newForw; 

  _rend = sqrt(xRear()*xRear() + yRear()*yRear()); 
  _phiend = _rear.phi(); 
  BbrAngle fPphi(_forward.phi()); 
  BbrAngle rPphi(_rear.phi()); 
  _twist = (fPphi - rPphi) * 0.5; 
 
  _sag += corr.getSag();

  _traj = DchSagTraj(_sag,_rear,_forward); 
}

void 
DchSWire::wireAlign( const DetAlignElem& align )
{
  wireAlign(align.transform());
}

void 
DchSWire::wireAlign( const HepTransformation& glob )
{

  HepPoint newRear = _rear.transform(glob);
 
  _rear = newRear; 
 
  HepPoint newForw = _forward.transform(glob); 
  
  _forward = newForw;  
 
  _rend = sqrt(xRear()*xRear() + yRear()*yRear());  
  _phiend = _rear.phi();  
  BbrAngle fPphi(_forward.phi());  
  BbrAngle rPphi(_rear.phi());  
  _twist = (fPphi - rPphi) * 0.5;  

  _traj = DchSagTraj(_sag,_rear,_forward); 
}

void
DchSWire::wireAlign( const HepTransformation& reartransf,   
		     const HepTransformation& forwtransf )
{
  HepPoint newRear = _rear.transform(reartransf); 
  _rear = newRear;  
  
  HepPoint newForw = _forward.transform(forwtransf);  
  _forward = newForw;   
  
  _rend = sqrt(xRear()*xRear() + yRear()*yRear());   
  _phiend = _rear.phi();   
  BbrAngle fPphi(_forward.phi());   
  BbrAngle rPphi(_rear.phi());   
  _twist = (fPphi - rPphi) * 0.5;   
  _traj = DchSagTraj(_sag,_rear,_forward);
}

Hep3Vector
DchSWire::yAxis( double z )
{
  Hep3Vector vec( xWireDC(z), yWireDC(z), z);
  return vec.unit();
}

void DchSWire::print(ostream &o) const 
{
  o << "rear end-plate coordinate:    " << *getRearPoint() <<"\n"
    << "forward end-plate coordinate: " << *getForwPoint() <<"\n"
    << "sagitta:                      " << getSag() <<"\n"
    << "radius on rear end-plate:     " << rEnd() <<"\n"
    << "radius at mid chamber:        " << rMid() <<"\n"
    << "z length:                     " << zLength() <<"\n"
    << "twist:                        " << dPhiz() <<"\n"
    << "stereo:                       " << stereo() <<"\n"
    << "wire id:                      " << Id() <<"\n"
    << "layer number:                 " << layer()->layNum() <<"\n"
    << "x-y rear:                     " << xRear() << " - " << yRear() <<"\n"
    << "x-y mid:                      " << xMid() << " - " << yMid() << endl;
}

void DchSWire::printInfo(ostream &o) const 
{
  o << "rear " << getRearPoint()->x() << " " << getRearPoint()->y() << " " 
    << getRearPoint()->z() <<"\n"
    << "forward " << getForwPoint()->x() <<" " << getForwPoint()->y() 
    << " " << getForwPoint()->z() << "\n"
    << "sagitta: " << getSag()
    << " z length: " << zLength()
    << " twist: " << dPhiz() 
    << " stereo: " << stereo() <<"\n"
    << "x-y mid: " << xMid() << " - " << yMid() << endmsg;
}

ostream&  operator << (ostream& o, const DchSWire& w) 
{ 
  w.print(o);
  return o;
}

