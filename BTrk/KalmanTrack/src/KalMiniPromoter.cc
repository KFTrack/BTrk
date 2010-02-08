//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniPromoter.cc,v 1.4 2007/07/09 21:54:59 brownd Exp $
//
// Description:
//   Promote KalMiniRep-based track fits to refit mode
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2004	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, Brian Aagaard
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalMiniPromoter.hh"
#include "TrkBase/TrkHotSelector.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalInterface.hh"
#include "ErrLogger/ErrLog.hh"

TrkErrCode 
KalMiniPromoter::promoteToHots(TrkRecoTrk* trk,
			       PdtPid::PidType hypo,
			       const TrkHotSelector* selector,
			       bool restorefail) {
  TrkErrCode retval(TrkErrCode::fail,101,"Can't attach KalMiniRep");
  static KalMiniInterface kmiface;
  static KalInterface kiface;
// 'null' hypo is keyword for default.  I hate secret codes!!!!
  if(hypo == PdtPid::null)
    hypo = trk->defaultType();
  if ( trk->attach(kmiface, hypo)) {
    KalMiniRep* minirep = kmiface.kalMiniRep();
    KalMiniRep::miniState oldstate = minirep->currentState();
    if(oldstate == KalMiniRep::hots){
      if(minirep->fitCurrent())
	retval = TrkErrCode(TrkErrCode::succeed,102,"Fit already in hots state");
      else
// attempt to fit
	retval = minirep->fit();
    } else {
      retval = minirep->changeState(KalMiniRep::hots);
      if( retval.success() ) {
// if requested, select hots
	if(selector != 0)selectHots(minirep,selector);
// perform the fit
	retval = minirep->fit();
	if (retval.failure() && restorefail ) {
// if the refit failed, try restoring the previous state
	  TrkErrCode fitResult = minirep->changeState(oldstate);
	  if( fitResult.success() )
	    fitResult = minirep->fit();
	  if (fitResult.failure()) {
// if the restore fails downgrade to seed
	    fitResult = minirep->changeState(KalMiniRep::seed);
	    fitResult = minirep->fit();
	  }
	}
      }
    }
// try the  regular KalInterface
  } else if(trk->attach(kiface, hypo))
// track has a regular Kalman rep; that's good enough
    retval = TrkErrCode(TrkErrCode::succeed,103,"Regular KalRep");
  return retval;
}

void
KalMiniPromoter::selectHots(KalMiniRep* minirep,const TrkHotSelector* selector) {
  TrkHotList* trkhotlist = minirep->hotList();
  if(trkhotlist != 0 && trkhotlist->hitCapable()){
    TrkHotList::nc_hot_iterator it=trkhotlist->begin();
    while(it!=trkhotlist->end()) {
      TrkHitOnTrk* hot = it.get();
      if (!selector->useHot(*hot))
	minirep->deactivateHot(hot);
      it++;
    }
  } else
    ErrMsg(error) << "MiniRep has no hit-capable hot list!" << endmsg;
}

KalRep* 
KalMiniPromoter::findKalRep(TrkRecoTrk* trk,
                            PdtPid::PidType hypo) {
  KalRep* retval(0);
  static KalMiniInterface kmiface;
  static KalInterface kiface;
// 'null' hypo is keyword for default.  I hate secret codes!!!!
  if(hypo == PdtPid::null)
    hypo = trk->defaultType();
  if ( trk->attach(kmiface, hypo)) {
    KalMiniRep* minirep = kmiface.kalMiniRep();
    KalMiniRep::miniState oldstate = minirep->currentState();
    if(oldstate == KalMiniRep::hots){
      retval = minirep->kalRep();
    } else {
      TrkErrCode change =  minirep->changeState(KalMiniRep::hots);
      if( change.success() ) {
        retval = minirep->kalRep();
      }
    }
// try the  regular KalInterface
  } else if(trk->attach(kiface, hypo))
// track has a regular Kalman rep; that's good enough
    retval = kiface.kalRep();
  return retval;
}
