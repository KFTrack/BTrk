//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: WireHit.cc 127 2010-09-13 08:34:00Z stroili $
//
// Description:
//	Class Implementation |WireHit|
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

#include "DchData/WireHit.hh"

#include "DchGeom/DchDetector.hh"
#include "DchGeom/DchLayer.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "ErrLogger/ErrLog.hh"
#include "DchCalib/DchTimeToDist.hh"
using std::endl;
using std::ostream;

WireHit::WireHit(double dist,double sigma){
  _dist=dist;
  _sigma=sigma;
}

//Destructor
WireHit::~WireHit() 
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

void 
WireHit::print( ostream& o ) const 
{
  o<<" dist "<<_dist<<" sigma "<<_sigma<<endl;
}

extern ostream& 
operator<<( ostream &o, const WireHit& hit ) 
{
  hit.print(o); return o;
} 

const GTrack* 
WireHit::getGTrack() const 
{
  return 0;
}

TrkDetElemId 
WireHit::elemId() const 
{
  return TrkDetElemId(DchCellAddr::cellIs(_wire,_layer),TrkDetElemId::null);
}

unsigned
WireHit::status() const
{
  return _status;
}

