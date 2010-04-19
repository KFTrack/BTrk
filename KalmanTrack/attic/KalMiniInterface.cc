//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniInterface.cc,v 1.2 2000/12/05 19:56:17 brownd Exp $
//
// Description:
//     class KalMiniInterface.
//     
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown, 11/14/00
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "ProxyDict/IfdIntKey.hh"

KalMiniInterface::KalMiniInterface()
{}

KalMiniInterface::~KalMiniInterface()
{}

const IfdKey& 
KalMiniInterface::myKey() const {
  static IfdIntKey _theKey(2154);
  return _theKey;
}

const KalMiniRep*
KalMiniInterface::kalmanMiniRep() const{
  return (const KalMiniRep*) myConstRep();
}

KalMiniRep*
KalMiniInterface::kalMiniRep() {
  return (KalMiniRep*) myRep();
}
