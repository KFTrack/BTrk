//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMaker.cc,v 1.51 2004/05/06 11:18:14 raven Exp $
//
// Description:
//   Creates tracks with KalReps.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 4/16/97
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalMaker.hh"
#include "KalmanTrack/KalCodes.hh"
#include "KalmanTrack/KalStub.hh"
#include "KalmanTrack/KalInterface.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "CLHEP/Alist/AList.h"
#include "TrkBase/TrkErrCode.hh"
#include "KalmanTrack/KalInterface.hh"
#include "ProxyDict/Ifd.hh"

//  Default constructor
//
KalMaker::KalMaker(const KalContext& context) : _kalcon(context){;}
KalMaker::~KalMaker(){;}
//
// Note that all KalMaker track/rep creation functions use the KalContext
// object to define what PID to assign by default.
//
TrkRecoTrk*
KalMaker::makeTrack(const TrkExchangePar& helix,
		    const TrkContext& con, double t0,
		    PdtPid::PidType hypo) const {
// check if the default hypo has been overridden
  if(hypo == PdtPid::null)
    hypo = _kalcon.defaultType();
//  Build a new track
  TrkRecoTrk* newtrack = createTrack(hypo, con, t0);
  if(newtrack != 0){
    KalRep* krep = new KalRep(helix, newtrack, _kalcon,hypo);
//  Install this rep in the track
    if(krep != 0)
      setRep(*newtrack,krep);
  }
//  Return it
  return newtrack;
}

TrkRecoTrk*
KalMaker::makeTrack(const TrkSimpTraj& seedtraj,
		    TrkHotList* hots,
		    const TrkContext& con, double t0,
		    PdtPid::PidType hypo,
		    double cfltlen,bool* cparams,
		    bool stealhots) const {
// check if the default hypo has been overridden
  if(hypo == PdtPid::null)
    hypo = _kalcon.defaultType();
//  Build a new track
  TrkRecoTrk* newtrack = createTrack(hypo, con, t0);
  if(newtrack != 0){
// build the KalRep
    KalRep* krep = new KalRep(seedtraj, hots,
			      newtrack, _kalcon, hypo,cfltlen,cparams,stealhots);
//  Install this rep in the track
    if(krep != 0)
      setRep(*newtrack,krep);
    else {
// we failed
      delete newtrack;
      newtrack = 0;
    }
  }
//  Return it
  return newtrack;
}


// create a new track from a KalStub
TrkRecoTrk* 
KalMaker::makeTrack(KalStub& stub) const {
  TrkRecoTrk* newtrack(0);
  if(stub.isPrimary()){
// get the KalRep from the stub
    const KalRep& srep = stub.kalRep();
// get the TrkRecoTrk from the rep
    const TrkRecoTrk* parent = srep.parentTrack();
// work with the tracks current default type
    PdtPid::PidType hypo = parent->defaultType();
// create a new track using this context
    newtrack = new TrkRecoTrk(*parent);
    if(newtrack != 0){
      TrkErrCode sappend(TrkErrCode::fail);
// copy
      KalInterface kinter;
      if(newtrack->attach(kinter,hypo)) {
	KalRep* krep = kinter.kalRep();
// append the stub on the new rep
	if(krep != 0)
	  sappend = krep->append(stub);
      }
// we failed; cleanup
      if(sappend.failure()){
	delete newtrack;
	newtrack = 0;
      }
    }
  }
//  Return it
  return newtrack;
}


TrkErrCode
KalMaker::changeFit(TrkRecoTrk& theTrack) const {
// work with the tracks current default type
  PdtPid::PidType dhypo = theTrack.defaultType();
  PdtPid::PidType hypo = _kalcon.defaultType();
// get helix parameters from the current default fit result (if possible)
  TrkRep* defrep = getRep(theTrack,dhypo);
  if(defrep != 0){
    TrkHotList *hots = defrep->hotList();
    assert(hots!=0);
    KalRep* krep = new KalRep(defrep->helix(0),*hots,&theTrack,_kalcon,hypo,true);
    if(krep != 0){
// copy over the history for the existing rep
      krep->addHistory(theTrack.status()->history());
      setRep(theTrack, krep);   // theTrack will delete all old Reps
      if(dhypo != hypo) changeDefault(theTrack,hypo);
    } else
      return TrkErrCode(TrkErrCode::fail,KalCodes::makerep,
			"Cannot create KalRep object");
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::fitresult,
		      "Cannot find initial track fit");
  return TrkErrCode(TrkErrCode::succeed);
}

TrkErrCode 
KalMaker::changeFit(TrkRecoTrk& theTrack,
		    const TrkSimpTraj& ctraj,
		    double cfltlen,
		    bool* cparams) const {
// work with the tracks current default type
  PdtPid::PidType dhypo = theTrack.defaultType();
  PdtPid::PidType hypo = _kalcon.defaultType();

  TrkRep* defrep = getRep(theTrack,hypo);
  if(defrep != 0){
    KalRep* krep = new KalRep(ctraj, defrep->hotList(), &theTrack,_kalcon,hypo,cfltlen,cparams,false);
//  Install this rep in the track
    if(krep != 0){
// copy over the history for the existing rep
      krep->addHistory(theTrack.status()->history());
      setRep(theTrack, krep);   // theTrack will delete all old Reps
      if(dhypo != hypo) changeDefault(theTrack,hypo);
    } else
      return TrkErrCode(TrkErrCode::fail,KalCodes::makerep,
			"Cannot create KalRep object");
  } else 
    return TrkErrCode(TrkErrCode::fail,KalCodes::fitresult,
		      "Cannot find initial track fit");
  return TrkErrCode(TrkErrCode::succeed);
}

// create the rep for this track
KalRep*
KalMaker::makeRep(TrkRecoTrk& theTrack,
                  const TrkExchangePar& helix,
                  const TrkHotList* hotlist) const {
  PdtPid::PidType hypo = theTrack.defaultType();
// This sets the rep PID to be that specified by the existing track
  KalRep *rep(0);
  if (hotlist!=0)
    rep = new KalRep(helix,*hotlist,&theTrack,_kalcon,hypo);
  return rep;
}


TrkErrCode 
KalMaker::addHypo(TrkRecoTrk& theTrack,
		  PdtPid::PidType hypo,
		  bool makeDefault) const {
  TrkErrCode retval(TrkErrCode::succeed);
// check the default hypo first
  PdtPid::PidType dhypo = theTrack.defaultType();
  KalInterface kinter;
  if(theTrack.attach(kinter,dhypo)) {
    if(dhypo == hypo)
      retval = TrkErrCode(TrkErrCode::succeed,KalCodes::hypoexists,
			  "Requested mass hypo is already a KalRep");
    else {
// get the KalRep from the interface
      KalRep* krep = kinter.kalRep();
// use it to clone the new rep
      TrkRep* newrep = (krep->cloneNewHypo(hypo));
      if(newrep != 0) {
// add this hypo to the track
	addHypoTo(theTrack,newrep,hypo);
// if requested, make this the new default PID
	if(makeDefault)
	  changeDefault(theTrack,hypo);
      } else
	retval = TrkErrCode(TrkErrCode::fail,KalCodes::makerep,
			    "Cannot create KalRep object");
    }
  } else
// can't do this if the existing default rep isn't a KalRep
    retval  = TrkErrCode(TrkErrCode::fail,KalCodes::notkalrep,
			 "Existing default rep isn't a KalRep");
  return retval;
}

    
TrkErrCode
KalMaker::fitHypo(TrkRecoTrk& theTrack,
		  PdtPid::PidType hypo) const {
// which hypo to fit
  if(hypo == PdtPid::null)
    hypo = theTrack.defaultType();
// attach the interface
  KalInterface kinter;
  if(theTrack.attach(kinter,hypo)) {
    KalRep* krep = kinter.kalRep();
    if(krep != 0)
// fit it
      return krep->fit();
    else
      return TrkErrCode(TrkErrCode::fail,KalCodes::norep,
			"Can't find specified KalRep");
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::notkalrep,
		      "Specified rep isn't a KalRep");
}
