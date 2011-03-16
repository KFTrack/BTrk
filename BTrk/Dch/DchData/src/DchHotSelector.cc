//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchHotSelector.cc 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//   See header file
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	UC, San Diego
//
// Author List:
//      Gerhard Raven 7/19/01
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "DchData/DchHotSelector.hh"
#include "DchData/DchHitOnTrack.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include <math.h>

DchHotSelector::DchHotSelector(bool activeOnly, double minDoca, double maxDoca) :
  _ignoreActive(!activeOnly),_minDoca(minDoca),_maxDoca(maxDoca)
{}

DchHotSelector::~DchHotSelector()
{}

bool
DchHotSelector::useHot(const TrkHitOnTrk& hot) const 
{
  const DchHitOnTrack *h = hot.dchHitOnTrack();
  bool use =  (h != 0) && (_ignoreActive || hot.isActive());
  double doca(0);
  if (use && (_minDoca>0 || _maxDoca>0)) {
    // note: if dcaWire is too large, the hit SHOULD be skipped;
    //       the uncertainty in the t->d relation is large beyond 0.9 cm
    //       if it's too small, the intrinsic resolution isn't that
    //       great, and the near-wire calibration isn't that great...
          doca = fabs(h->dcaToWire());
          use = ((_minDoca>0?doca>_minDoca:true)
              || (_maxDoca>0?doca<_minDoca:true));
  }
  return use;
}

