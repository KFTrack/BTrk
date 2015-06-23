//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkActiveHotSelector.cc,v 1.1 2001/07/18 00:42:42 brownd Exp $
//
// Description:
//   class TrkActiveHotSelector.  Trivial implementation of TrkHotSelector for
//   selecting active hots.
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
#include "TrkBase/TrkActiveHotSelector.hh"
#include "TrkBase/TrkHitOnTrk.hh"

TrkActiveHotSelector::TrkActiveHotSelector(){}

TrkActiveHotSelector::~TrkActiveHotSelector(){}

bool
TrkActiveHotSelector::useHot(const TrkHitOnTrk& hot) const {
  return hot.isActive();
}
