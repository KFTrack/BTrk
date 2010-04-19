//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairMaker.cc,v 1.11 2000/07/17 18:27:03 roberts Exp $
//
// Description:
//   Creates tracks with KalPairReps.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 3/21/98
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalPairMaker.hh"

#include "ErrLogger/ErrLog.hh"
#include "KalmanTrack/KalInterface.hh"
#include "KalmanTrack/KalPairRep.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkFitter/TrkContextCustom.hh"

KalPairMaker::KalPairMaker() {;}
KalPairMaker::~KalPairMaker() {;}

TrkRecoTrk*
KalPairMaker::makePairTrack(TrkRecoTrk& oldtrk, KalPairConduit* conduit) {
  // oldtrk must have a KalRep

  KalInterface kalInter;
  TrkRecoTrk* newTrk(0);
  PdtPid::PidType ptype = oldtrk.defaultType();
  if( oldtrk.attach(kalInter, ptype) ) {
    // Steal the context from the old track
    TrkContextCustom context(oldtrk);
    // create the new track
    newTrk = createTrack(ptype, context, oldtrk.trackT0());
    // construct the pair rep from the KalRep
    KalPairRep* pairRep =
      new KalPairRep(newTrk, kalInter.kalmanRep(), conduit);
    //  Install this rep in the track
    setRep(*newTrk,pairRep);
  }
  else
    ErrMsg(error) << " No KalRep in track" << endmsg;
  
  return newTrk;
}
