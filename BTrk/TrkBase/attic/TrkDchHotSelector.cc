//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDchHotSelector.cc,v 1.1 2001/07/18 00:42:42 brownd Exp $
//
// Description:
//   class TrkDchHotSelector.  Trivial implementation of TrkHotSelector for
//   selecting Dch hots.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 7/17/01
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkDchHotSelector.hh"
#include "TrkBase/TrkHitOnTrk.hh"

TrkDchHotSelector::TrkDchHotSelector(bool activeOnly) :
  _ignoreActive(!activeOnly)
{}

TrkDchHotSelector::~TrkDchHotSelector()
{}

bool
TrkDchHotSelector::useHot(const TrkHitOnTrk& hot) const {
  return (hot.dchHitOnTrack() != 0) &&
    (_ignoreActive || hot.isActive());
}

