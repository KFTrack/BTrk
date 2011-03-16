//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHit.cc 127 2010-09-13 08:34:00Z stroili $
//
// Description:
//	Class Implementation |DchHit|
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      S. Schaffner - Original
//	A. Snyder - Mods
//      S. Sen                  Feb 17 1998     New digi()->MCPtr() references
//
// Copyright Information:
//	Copyright (C) 1996	SLAC
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"

#include "DchData/DchHit.hh"


#include "DchData/DchDigi.hh"
#include "DchData/DchHitData.hh"
#include "DchData/DchDigiWF.hh"
#include "DchData/DchDigiMC.hh"
#include "DchData/DchDigiStatus.hh"
#include "DchGeom/DchDetector.hh"
#include "DchGeom/DchLayer.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "ErrLogger/ErrLog.hh"
#include "DchCalib/DchTimeToDist.hh"
using std::endl;
using std::ostream;


DchHit::DchHit( const DchDigi& aDigi,
                const DchDetector& det,
                const DchTimeToDistList &t2d,
                unsigned iTdc)
  : DchHitBase(aDigi.layernumber(), aDigi.wirenumber(),aDigi.TdcTime(iTdc), iTdc,aDigi.charge(),det,t2d),
    _digiPtr(&aDigi),
    _phi(_layerPtr->phiWire(_wire)),_cosphi(cos(_phi)),_sinphi(sin(_phi)),
    _rmid(_layerPtr->rMid()),_zlen(_layerPtr->zLength()),
    _status(aDigi.status())
{
  // check WF status
  unsigned status = aDigi.status() & 0xe00;
  unsigned badWF     =  status & DchDigiStatus::badWaveform;
//  unsigned badPed    =  status & DchDigiStatus::badPedestal;
//  unsigned clippedWF =  status & DchDigiStatus::clippedWaveform;
  
  if ( badWF !=0  ) _charge = 0.;
  
  if ( status!=0 ) {
    const DchDigiWF *wv = aDigi.waveform();
  // if you have the waveform try to do something...
    if (wv!=0) _charge = wv->recomputeCharge(aDigi.status(), _charge);
  }
}

DchHit::DchHit( const DchHitData& hitData,
                const DchDetector& det,
                const DchTimeToDistList &t2d)
  : DchHitBase(hitData.layernumber(), hitData.wirenumber(),hitData.tdcTime(),hitData.tdcIndex(),hitData.charge(),det,t2d),
    _digiPtr(0),
    _phi(_layerPtr->phiWire(_wire)),_cosphi(cos(_phi)),_sinphi(sin(_phi)),
    _rmid(_layerPtr->rMid()),_zlen(_layerPtr->zLength()),
    _status(hitData.status())
{
}

DchHit::DchHit(const DchHit& other) :
  DchHitBase(other),_digiPtr(other._digiPtr),_phi(other._phi),
  _cosphi(other._cosphi),_sinphi(other._sinphi),_rmid(other._rmid),
  _zlen(other._zlen),_status(other._status)
{
}

DchHit&
DchHit::operator = (const DchHit& other) {
  if(&other != this){
    DchHitBase::operator =(other);
    _digiPtr     = other._digiPtr;
    _phi         = other._phi;
    _cosphi      = other._cosphi;
    _sinphi      = other._sinphi;
    _rmid        = other._rmid;
    _zlen        = other._zlen;
    _status      = other._status;
  }
  return *this;
}


//Destructor
DchHit::~DchHit() 
{
  // This is ugly and inefficient.  This, along with the rest of 
  //  the hitList mess, should be cleaned up by handling the 
  //  association in an external map

  // Not written as a loop because removeHit() modifies TrkFundHit::_hitList
  short count = 0;
  while (nUsedHits() > count) {
    bool removed = _hitList[count]->parentTrack()->hits()->removeHit(this);
    if (!removed) count++;
  }
}

bool
DchHit::operator==( const DchHit& rhs ) const 
{
  return (this == &rhs);
}

void 
DchHit::print( ostream& o ) const 
{
  o << "layer:" <<  _layer 
    << "\nwire:" << _wire 
    << "\nrawTime:" << rawTime()  << " charge:" << charge() 
    << "\nphi:" << _phi << " cos(phi):" << _cosphi << " sin(phi):" 
    << _sinphi 
    << "\nrmid:" << _rmid << " zlen:" << _zlen 
    << endl;
}

double 
DchHit::x( double z ) const 
{
  return geom()->xWire(layernumber(),wire(),z);
}

double 
DchHit::y( double z ) const 
{
  return geom()->yWire(layernumber(),wire(),z);
}

extern ostream& 
operator<<( ostream &o, const DchHit& hit ) 
{
  hit.print(o); return o;
} 

const GTrack* 
DchHit::getGTrack() const 
{
  const DchDigiMC *d=digi()->MCPtr();
  return d==0?0:d->getGTrack(0);
}

TrkDetElemId 
DchHit::elemId() const 
{
  return TrkDetElemId(DchCellAddr::cellIs(_wire,_layer),TrkDetElemId::dch);
}

unsigned
DchHit::status() const
{
  return _status;
}

// FIXME: cannot inline because DDL compiler cannot handle DchCalib/DchTimeToDist.hh...
double 
DchHit::driftDist(double arrivalTime, int ambig, double entranceAngle,
                   double dipAngle, double z) const
{ 
  return _t2d->timeToDist(driftTime(arrivalTime),ambig, entranceAngle,dipAngle,z,charge()); 
}
