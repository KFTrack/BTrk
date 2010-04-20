//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkSvtHotSelector.cc,v 1.1 2001/07/18 00:42:43 brownd Exp $
//
// Description:
//   class TrkSvtHotSelector.  Trivial implementation of TrkHotSelector for
//   selecting Svt hots.
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
#include "TrkBase/TrkSvtHotSelector.hh"
#include "TrkBase/TrkHitOnTrk.hh"

TrkSvtHotSelector::TrkSvtHotSelector(bool activeOnly) :
  _ignoreActive(!activeOnly)
{}

TrkSvtHotSelector::~TrkSvtHotSelector()
{}

bool
TrkSvtHotSelector::useHot(const TrkHitOnTrk& hot) const {
  return (hot.svtHitOnTrack() != 0) &&
    (_ignoreActive || hot.isActive());
}

