//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalHitCache.cc,v 1.6 2002/07/15 20:31:49 echarles Exp $
//
// Description:
//      class KalHitCache.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown10/20/98
//------------------------------------------------------------------------

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalHitCache.hh"
#include "BTrk/KalmanTrack/KalHit.hh"


KalHitCache::KalHitCache(const TrkHitUse& hituse,
			 KalHit* kalhit, double chisq) :
  _kalhit(kalhit),
  _hituse(&hituse),
  _chisquared(chisq),
  _chi2Valid(true)
{;}

KalHitCache::KalHitCache() : _kalhit(0), _hituse(0), _chisquared(-1.0), _chi2Valid(false)
{;}

KalHitCache::KalHitCache(const KalHitCache& other) :
  _kalhit(other._kalhit),
  _hituse(other._hituse),
  _chisquared(other._chisquared),
  _chi2Valid(other._chi2Valid)
{;}

KalHitCache::~KalHitCache()
{;}

KalHitCache::KalHitCache(const KalHitCache& other,const KalRep* rep) :
  _kalhit((KalHit*)(other._kalhit->clone(rep))),
  _hituse(other._hituse),
  _chisquared(other._chisquared),
  _chi2Valid(false)
{;}

KalHitCache&
KalHitCache::operator = (const KalHitCache& other) {
  if(&other != this){
    _kalhit = other._kalhit;
    _hituse = other._hituse;
    _chisquared = other._chisquared;
    _chi2Valid = other._chi2Valid;
  }
  return *this;
}

bool
KalHitCache::operator == (const TrkHitUse& hituse) const {
// use hituse equivalence operator eventually
//  return *_hituse == hituse;
// for now, use pointer comparison
  return _hituse == &hituse;
}

bool
KalHitCache::operator == (const KalHitCache& other) const {
  return _hituse == other._hituse &&
    _kalhit == other._kalhit;
}

void
KalHitCache::deleteAll() {
  if(_kalhit != 0)
    _kalhit->deleteHOT();
  delete _kalhit;
}

