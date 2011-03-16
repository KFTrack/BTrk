//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHitBase.cc 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//	Class Implementation |DchHitBase|
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
#include "DchData/DchHitBase.hh"

#include "DchGeom/DchDetector.hh"
#include "DchCalib/DchTimeToDistList.hh"
#include "DchCalib/DchTimeToDist.hh"
#include "DchGeom/DchLayer.hh"
#include "ErrLogger/ErrLog.hh"



//Declare static members
double DchHitBase::_addFindingSig = 200.e-4;
double DchHitBase::_fittingSigMult = 1.0;


DchHitBase::DchHitBase( unsigned layer, unsigned wire, double rawTime, unsigned iTdc, double charge,
                        const DchDetector& det,
                        const DchTimeToDistList& t2d)
  : _rawTime(rawTime), _charge(charge), _geomPtr(&det),  _layerPtr(det.getDchLayer(layer)), _t2d(t2d.getTimeToDist(layer,wire)), _layer(layer), _wire(wire), _iTdc(iTdc)
{
  assert( _layerPtr!=0 && _layerPtr->exist() );
  assert( 0 <= wire && wire < _layerPtr->nWires() );
}

DchHitBase::~DchHitBase()
{
}

const Trajectory*
DchHitBase::hitTraj() const
{
   return layer()->makeHitTrajInGlobalCoords(wire(),0.0);
}

int 
DchHitBase::whichView() const
{
   return _layerPtr->view(); 
}
